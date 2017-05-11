% prunes the tree exactly at the timepoints where the condition is fullfilled the first time
% (the timepoint where the condition is fullfilled the first time stays in the tree)
function [prunedTree, leaves_ofPruned ]= tUtil_pruneTree_TIME(currentTree,pruningCond, startingCell_dirtyHack)

% tmp = precalcCellIndexing(currentTree); %very slow


   leaves_ofPruned = [];
   
   if ~exist('startingCell_dirtyHack','var')
        nodes = [1];
        assert(any(currentTree.cellNr==1),'pruning function doesnt work for cells without a root')
   else
       st=dbstack;
       if ~strcmp(st(2).name, 'spliceTree') %this funtioncality is only supoosed to be used from within this spliceTree fcuntion, if called from another function issue this warning
        warning('this is some dirty hack to make pruneing work for trees that dont have root=1. not really tested, only used in treeSplicing')
       end
       nodes = startingCell_dirtyHack;
   end
   
   retainedNodes = [];
   
   
   retained_indices = zeros(size(currentTree.cellNr)); %which datapoints are retained after pruning

   pC = pruningCond(currentTree);
   assert(length(pC)==length(currentTree.cellNr)); % the prunding condition must return one value per timepoint for a given cell
   
   allCells = unique(currentTree.cellNr);

   %for speed precalc the cellNr == cC operatoins
%    [precalc_cellix, precalc_cc ]= precalcCellIndexing_faster(currentTree);
%    sparse_precalc_cellix = sparse(precalc_cellix);
   while (~isempty(nodes))
       nextLevel = [];
       for i=1:length(nodes)
           currentNode = nodes(i);
           
           assert( currentNode ~= 0)
           
           retainedNodes= [retainedNodes currentNode];
           
           cellix = currentTree.cellNr == currentNode;
%            cellix2 = precalc_cellix(precalc_cc==currentNode,:);
%            dummyix = unidrnd(length(precalc_cc));
%            cellix3 = precalc_cellix(dummyix,:);
%            cellix4 = sparse_precalc_cellix(dummyix,:);
           
           
           pruneThisCell = any(pC & cellix); %pruning cond fulfilled in this cell?
           % if not, process children
           if ~pruneThisCell 
               
               % all timepoints of the current cell are retained
               retained_indices = retained_indices+cellix;
               
               %left
               if any(allCells == 2.*currentNode)
                   nextLevel = [nextLevel 2.*currentNode ];
               end
               %right
               if any(allCells== 2.*currentNode+1)
                   nextLevel = [nextLevel 2.*currentNode+1 ];
               end
               
           else
               % only those datapoints are retained up to (and including) the time when the
               % condition is fulfilled
               firstPoint_ofCell = find(cellix,1,'first');
               first_pC_point = find(pC & cellix,1,'first');
               retained_indices(firstPoint_ofCell:first_pC_point)=1;
               
               %remember the leave
               leaves_ofPruned = [leaves_ofPruned currentNode];
           end
       end
          
       nodes = nextLevel;
   end
   assert(all(retained_indices<2))
   
   %% put together the pruned tree from the retained indices   
   prunedTree = tUtil_treeIndexing(currentTree,retained_indices==1);
end


function [cellIndices, cells] = precalcCellIndexing(tree)
error('dont use; too slow')
cells = unique(tree.cellNr);

cellIndices = zeros(length(cells),length(tree.cellNr)); % #uniqueCells x #datapoints
[~,ix] = ismember(tree.cellNr,cells); % calculates for each datapoint which index in cells it corresponds to: if we are at cellNr(i) we write to cellIndices(ix(i))
for i =1:length(tree.cellNr)

    cellIndices(ix(i),i) = 1;
end
    
end

function [cellIndices, cells] = precalcCellIndexing_faster(tree)

    cells = unique(tree.cellNr);

    cellIndices = false(length(cells),length(tree.cellNr)); % #uniqueCells x #datapoints
    [~,ix] = ismember(tree.cellNr,cells); % calculates for each datapoint which index in cells it corresponds to: if we are at cellNr(i) we write to cellIndices(ix(i))

    ixlinear = sub2ind(size(cellIndices),ix,1:length(tree.cellNr)); % calculates to positions where the matrix should have 1s it replaces : "for i: cellIndices(ix(i),i) = 1;end"
    cellIndices(ixlinear)=true;
end