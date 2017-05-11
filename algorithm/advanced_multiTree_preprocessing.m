%%
% !!compared to the original method (multiTree_preprocessing), a clever way of avoiding the enumeration of all trees!!
% we do the same preprocessing for the subtrees (calc leaves, diffpoints, cells), but avoid
% enumerating all combinatorics. instead we do "precalc_4_prime" to speed up the calculation of the
% pprime routine
% 
function [transformedOverall ]=  advanced_multiTree_preprocessing(trees,WL,timeframe,tree_startingTimes,commonTree_struct)

%% ======some sanity checks on the inputs==========================

%make sure both optional args are there or non
    assert((~exist('tree_startingTimes','var') & ~exist('commonTree_struct','var'))| (exist('tree_startingTimes','var') & exist('commonTree_struct','var')),'the last two args have both to be specified or non of them')

    if exist('tree_startingTimes','var')
        %check by filename that the subtree->commonUpTree relation is ok
        check_subtree_uptree_relation(trees,commonTree_struct)
    end
    
    %optinal args in case of subtrees
    if ~exist('tree_startingTimes','var')
        tree_startingTimes = zeros(size(trees));
        disp('assuming 0 starting time for all trees: they are centred')
    end

    %optinal args in case of subtrees
    if ~exist('commonTree_struct','var')
        commonTree_struct.trees = {[]};
        commonTree_struct.index2commonUpperTree = ones(1,length(trees));
        disp('assuming 0 starting time for all trees: they are centred')
    end
    assert(length(tree_startingTimes)==length(trees),'starting times array must have same length as the tree array')
    assert(length(commonTree_struct.index2commonUpperTree)==length(trees),'commonTree_struct.index2commonUpperTree must have same length as trees array')

%% =============================================


%=======================
% CHANGED 28.10
% to keep trakc of the commonUPTrees, here all used cells are remembered (initialize with [])
commonTree_struct.usedCells = cell(1,length(commonTree_struct.trees));
for i = 1:length(commonTree_struct.trees)
    commonTree_struct.usedCells{i} = [];
end

%*********** pruning to onset
pC_time = @(t)t.(WL)==1;
g = @(mT)tUtil_pruneTree_TIME(mT,pC_time);
mTrees = cellfun(g, trees,'UniformOutput',0);

assert(all(cellfun(@(t)t.(timeframe)(1)==0,mTrees)),'All trees must start with timepoint 0, just rescale them!')

transformed = cell(1,length(mTrees));

%this is the loop over all observed trees
textprogressbar('tranforming trees ')
for i =1:length(mTrees)
    textprogressbar(i*100/length(mTrees))
    sanityCheck(mTrees{i},WL)
    
    %===============PART 1=================
    all_dp = sort(unique(mTrees{i}.cellNr));
    subtrees = cell(1,length(all_dp));
    %iterate over all subtrees (number of subtrees = number of cells)
    for j=1:length(all_dp)
        currentDP = all_dp(j); %decision in cell currentDP
        
        % den entsprechenden subtree holen
        c_subtree = tUtil_getSubtree(mTrees{i},currentDP);
        
        %let the rootnode be cell nr 1
        assert(currentDP == min(c_subtree.cellNr), 'the subtree root is not the diff point, but it should be'); %added when exchangeing inversetransformCellNrs() for tUtil_renumberTree(). makes sure the subtree is only starting the the diff point (this is enforce by inverse... but not be renumber)
        
        oldRoot = currentDP; % for speed give the old root to tUtil_renumberTree
        c_subtree = tUtil_renumberTree(c_subtree,1, oldRoot); % relabel the subtree so the root is cell nr 1

        %also create all variants of this subtree that have different root
        %length, needed for the integral thing
        [trees_cutAtdiffPoint, diffPoint] = createTrees_varyingRootLength(c_subtree);
        
        subtreeStruct= struct();
        subtreeStruct.diffPoints = diffPoint;
        % subtreeStruct.leaves = tUtil_getLeaves(c_subtree); % this one is way slower than below
        subtreeStruct.leaves = tUtil_getLeaves(c_subtree,1);
        subtreeStruct.allCells = sort(unique(c_subtree.cellNr));

        % GET RID OF ORIGINAL SUBTREE: instead put in some fields required later
        subtreeStruct.first_absoluteTime = c_subtree.absoluteTime(1);
        subtreeStruct.(['hasOnset_' WL]) = arrayfun(@(c)any(c_subtree.(WL)(c_subtree.cellNr==c)), subtreeStruct.allCells); % for each cell determine if it has an onset
        %sanity check: all onsets must be in leaves
        assert(all(ismember(subtreeStruct.allCells(subtreeStruct.(['hasOnset_' WL])), subtreeStruct.leaves)))
        birth = arrayfun(@(c)min(c_subtree.absoluteTime(c_subtree.cellNr==c)), subtreeStruct.allCells);
        death = arrayfun(@(c)max(c_subtree.absoluteTime(c_subtree.cellNr==c)), subtreeStruct.allCells);
        subtreeStruct.lifetime = death-birth;
        

        subtrees{j} = subtreeStruct;
    end
    
    %===============PART 2=================
    % do the precalculation for the pprime
    pprime_Structure = precalc_4_Prime(mTrees{i},WL);
    
    %==========PART 3=============
    %assemble the precalculated data into a struct:
    % 1. the subtrees of all cells (including their leaves, diffpoints etc)
    % 2. the structure for the pprime calculation
    
    %for some stupid reason struct( ...) tries to create an array! using
    %the old way of creating a struct
    transformed{i}.subtreeRoots = all_dp;
    transformed{i}.subtrees = subtrees;
    transformed{i}.pprime_Structure = pprime_Structure;
    
%====================
% CHANGED 28.10
% IN CASE WE CUT SUBTREES from the full original tree,
%
% we add information when the tree actually started, 
    transformed{i}.startingTime = tree_startingTimes(i);
% additionally we account for the time before, considering all cells obsered previously,
% but also take care that we dont use the same "prev." cells twice
% -> see projects/Trees/code/MATLAB/MarkerDelay/linearDiff/testUncenteredTrees/README
    upTree_ix = commonTree_struct.index2commonUpperTree(i);
    associatedUpTree = commonTree_struct.trees{upTree_ix};
    if ~isempty(associatedUpTree) % is emtpy when we actually dont have subtrees at all OR if the particular tree is really the root of the large tree
        %just take the rootpath
        RP = rootpath(mTrees{i}.originalRoot);
        RP = setdiff(RP,mTrees{i}.originalRoot); %remove the cell itself
                
        % check whioch of those cells have already been used
        currentUsedCells = commonTree_struct.usedCells{upTree_ix};
        unusedCells = setdiff(RP,currentUsedCells); % contained in the rootpath but not yet used
        
        %if all cells are used up, nothing has to be accounted in this cell
        % eg if subtrees 2,3 are selected, when the method arives at subtree 3, cell 1 is already
        % used up -> nothing left!
        if isempty(unusedCells)
            %use this 0 length interval (in thelikelihood this should only become a factor of 1)
            transformed{i}.commonUpTreeCorrection_start_stop = [ tree_startingTimes(i) tree_startingTimes(i)];
            continue;
        end
        % the unused cells form a branch that is undifferentiated
        % we account for this by including a correction, spanning the time from beginning of the
        % branch to its end
        ix_cells = ismember(associatedUpTree.cellNr,unusedCells);
        branchTime = associatedUpTree.(timeframe)(ix_cells);
        transformed{i}.commonUpTreeCorrection_start_stop = [min(branchTime) max(branchTime)];
        
        %mark them as used
        commonTree_struct.usedCells{upTree_ix} = union(commonTree_struct.usedCells{upTree_ix},unusedCells );
        
    else % otherwise use this 0 length interval (in thelikelihood this should only become a factor of 1)
        transformed{i}.commonUpTreeCorrection_start_stop = [ tree_startingTimes(i) tree_startingTimes(i)];
    end
end
textprogressbar(' done ')

%% figure out the different timepoints where we need to evaluate the CME later
disp('calculating timerange')

timeRange_firstProcess = gatherTimepoints_acrossTrees(trees);

% this one here takes loads of time, but returns just the same as the above
% NO IDEA WHY I DID IT LIKE THIS
% timerange is essentiallu just all unique timepoints over all trees
%timeRange_secondProcess = secondProcess_gatherTimepoints_acrossTrees(transformed,WL);
timeRange_secondProcess = timeRange_firstProcess;

transformedOverall.transformedTrees = transformed;
transformedOverall.timeRange_firstProcess = timeRange_firstProcess;
transformedOverall.timeRange_secondProcess = timeRange_secondProcess;

    
end

function timeRange = gatherTimepoints_acrossTrees(trees)

    timeRangeArray = cell(1,length(trees));
    for i= 1:length(trees)
        cT = trees{i};
        cells = unique(cT.cellNr);
        
        t_perCell = [];
        for j = 1:length(cells)
            t = cT.absoluteTime(cT.cellNr==cells(j));
            t_perCell= [t_perCell t-t(1)];
        end
        timeRangeArray{i} = sort(unique(t_perCell));
    end
    
    timeRange = sort(unique([timeRangeArray{:}]));
end

%% ===================================
% figures out at which timepoints we have to solve the CME of the
% differentiaion process!!
% 
% we require it at: the end of each cellcycle -> figure out all cellcycles
% but we also require it at all timepoints of the deciding cells
function timeRange = firstProcess_gatherTimepoints_acrossTrees(transformedTrees)
    error('deprecated. relies on originalSubtree, which is unnecessary. use (gatherTimepoints_acrossTrees)')
    timeRangeArray = cell(1,length(transformedTrees));
    for i= 1:length(transformedTrees)
        

        %% as we calc the transition matrix for the first process in this function,
        %but also need it later to assign probs to the roots of the delay tree (contain still some undiff. timepoints)
        % we also have to feed in the timepoints of the deciding cells to get
        % transition matrix also at these timepoints
        decindingCellTPs = [];
        for j = 1:length(transformedTrees{i}.subtrees)
           cT =  transformedTrees{i}.subtrees{j}.originalSubtree;
           time = cT.absoluteTime(cT.cellNr==1);
           
           %make time relative to birth
           time_rel = time-time(1);
           decindingCellTPs = [decindingCellTPs time_rel];
        end
        decindingCellTPs = unique(decindingCellTPs)    ;
        
        %% now find the different cell cycle length!!
        % kind of stupid here, one could just take the observed tree and calc its cellcyle
        % but i just leave it as it was before...
        cellCycleLengths = [];
        for j = 1:length(transformedTrees{i}.subtrees)
            %checkIntegrity(undiff_trees{i},WL,timescale,'absoluteTime')
            cellCycles = findcellcyle(transformedTrees{i}.subtrees{j}.originalSubtree);
            cellCycleLengths = [cellCycleLengths cellCycles];
        end
        
        %evaluate the process at all these times + the timepoints of deciding cells
        %(last variable of the function call)
        timeRangeArray{i} = sort(unique([0 cellCycleLengths decindingCellTPs]));
    end
    
    %throw it all together
    timeRange = sort(unique([timeRangeArray{:}]));

end

%% ================================
%figures out at what timepoints we need to solve the CME for the delay
%process
%
% we require it at all cell cycle length of the delay trees
% + the cell cycle of the root cell can vary from 0 to full length
function timeRange = secondProcess_gatherTimepoints_acrossTrees(transformedTrees,WL)


    timeRangeArray = cell(1,length(transformedTrees));
    for outerCounter= 1:length(transformedTrees)

        subtrees = transformedTrees{outerCounter}.subtrees;
        
        %% lets crate all the subtrees, by cutting the root of the subtrees at any timepoint
        % note that for each subtree we have to consider different points of onset
        % within the rootcell, -> j-loop
        subtrees_diffPoints= cell(1,length(subtrees)); % here we memorize for each subtree the possible diff poitns, relative to its birth
        subtrees_shortened= cell(1,length(subtrees)); % here are the corresponding shortened trees
        subtrees_Index = cell(1,length(subtrees)); %this will help to get the linearized version (eg of the pvalues) back into the array form
        for i =1:length(subtrees)

            [trees_cutAtdiffPoint, diffPoint] = createTrees_varyingRootLength(subtrees{i}.originalSubtree);
            subtrees_diffPoints{i} = diffPoint; %relative to the first timepoint!!
            subtrees_shortened{i} = trees_cutAtdiffPoint;
            subtrees_Index{i} = i.*ones(1,length(trees_cutAtdiffPoint));

        end
        
        % having create all subtrees with different rootlength, thrwo them
        % all together and figure out the cell cycle times
        trees = [subtrees_shortened{:}];
        
        cellCycleLengths = [];
        for i = 1:length(trees)
            checkTruncationOfDelayTrees(trees{i},WL)
            %checkIntegrity(undiff_trees{i},WL,timescale,'absoluteTime')
            cellCycles = findcellcyle(trees{i});
            cellCycleLengths = [cellCycleLengths cellCycles];
        end
        
        timeRangeArray{outerCounter}=cellCycleLengths;
    end
    
    timeRange = sort(unique([0 timeRangeArray{:}]));
end


%% ======================================
% given on tree, create varaitions of this one by changing the root cells
% length
% return these trees (their first timepoint is now 0!) and also return the
% timepoint (in the original tree) where we cut it
function [trees_cutAtdiffPoint diffPoint] = createTrees_varyingRootLength(tree)
    currentST = tree;
    currentRoot_time = currentST.absoluteTime(currentST.cellNr==1);
    %consider each timpeoint as possible diffpoint
    %however skipt the last point, causes trouble ATM
    OLD =0;
    if OLD ==1
        warning('kiping last tp as diff Point due to some bugs')
        currentRoot_time = currentRoot_time(1:end-1);
    else
        currentRoot_time = currentRoot_time;
    end
    trees_cutAtdiffPoint = cell(1,length(currentRoot_time)); %here we store the subtree, cut a diff points
    for j=1:length(currentRoot_time)
        currentDiffPoint = currentRoot_time(j);
        currentUndiffIx = currentST.absoluteTime<currentDiffPoint & currentST.cellNr==1;
        %kick out all these timepoint to get a tree with shorter root
        tmp  = tUtil_treeIndexing(currentST,~currentUndiffIx);
        tmp.absoluteTime = tmp.absoluteTime-currentDiffPoint; %make the time in this tree relative to the first TP
        trees_cutAtdiffPoint{j} = tmp;
    end
    diffPoint = currentRoot_time-currentRoot_time(1);
end

%% ======================================
function checkTruncationOfDelayTrees(tree,wavelength)
cells = unique(tree.cellNr);
for i = 1:length(cells)
    
   ix=tree.cellNr ==  cells(i);
   %wenn der marker an ist, dann nur zum letzten zeitpunk
   markerData = tree.(wavelength)(ix);
   blub = find(markerData==1,1,'first');
   
   %entweder gar nicht oder als letztes
   assert(isempty(blub)||blub==length(markerData),'make sure the tree is already truncated, that is once the marker goes on, no more timepoints afterwards')
end
end

%% ======================================
%returns all cells in the tree that are not contained in diffpoints and diffpoints subtrees
function undiffCells = getUndiffCells(tree,diffpoints)
    overallCells = unique(tree.cellNr);
    subtreeCells = [];
    for i =1:length(diffpoints)
       subtree = tUtil_getSubtree(tree,diffpoints(i)) ;
       subtreeCells = [subtreeCells unique(subtree.cellNr)];
    end
    
    undiffCells = setdiff(overallCells,subtreeCells);
end


function sanityCheck(michiTree,WL)
    markerCells = unique(michiTree.cellNr(michiTree.(WL)==1));
    
    for i =1:length(markerCells)
       
        currentCell = markerCells(i);
        
        %es darf keine tochterzelle geben
        tmp_l = any(michiTree.cellNr == currentCell.*2); 
        tmp_r = any(michiTree.cellNr == currentCell.*2+1); 
        
        assert(~tmp_l && ~tmp_r)
    end  
end