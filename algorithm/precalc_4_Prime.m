% in order to speed uup things we precalculate some quantities for the calculation of PPRime
% that do not change at all for a given tree!
% this is 
%   * the cell lifetimes
%   * the cells contributing their hazard to PPrime(i)

function precalcStructure = precalc_4_Prime(tree,markerField)

precalcStructure = [];
% first calc the birth and death times of all cells in the trees
% to speed up the hazard calculation per cell

cells = sort(unique(tree.cellNr));

%birth and divions times of cells RELATIVE TO TREEROOT (nedd to add starting time when used with noncentered trees)
birthtimes = zeros(1,length(cells));
deathtimes = zeros(1,length(cells));
for i = 1:length(cells)
    cC = cells(i);
    t = tree.absoluteTime(tree.cellNr==cC);
    deathtimes(i) = max(t);
    birthtimes(i) = min(t);
end

leaves = tUtil_getLeaves(tree);
%% ========================
% now that we have it for each cell, combinethem into PPrime(i)
% which is the prob to be undiff until (and including )i
% however accounting only for left branches
% that is, calc the prb on the rootpath until the current cell is the right sister somwhere
% (however, we also have to check the the left sister indeed exists!)

% one entry in PPrime cosists of the product of hazard rates of upstream cell (so that nothing is overlapping)
% determine which cells contribute to each entry 
relevantCells_perEntry = cell(1,length(cells));
relevantCells_perEntry_IX = cell(1,length(cells)); % even more usefull: determine the indices in cells 

isRootFlag = []; % does this entry contain the root of the tree? if so we also include the startingTime thing later in the calculation
for i = 1:length(cells)
    tmp = intersect(rootpath(cells(i)),cells); % only take the rootpath up to the root (which is not nec 1)
    
    %first condition is being a right sister, second is that its left sister does exist
    isRightSister_hasLeftSister= mod(tmp,2)~=0  & ismember(tmp-1,cells); %rootpath(i)-1 = left sisters of rootpath(i) if rootpath i is a right sister
    
    ix = find(isRightSister_hasLeftSister,1,'last'); % the first(from the leave) cell that is a right sister and does have a left sister
    
    if isempty(ix) % no cell could be found on the rootpath that has a left sister (without being theleft sister itself). hence we want all cells from leave to root
        ix =1;
        isRootFlag(i)=true; %set a flag for later
    else
        isRootFlag(i)=false;
    end
    
    
    % the path from the first right sister to the elave, including the leave!
    % in case the current cells is already the right sister, its just empty
    relevantCells = tmp(ix:end);
    relevantCells_perEntry{i} = relevantCells;
    [dummy, ix] = ismember(relevantCells, cells); assert(all(dummy));
    relevantCells_perEntry_IX{i} = ix;
end

%% SOME STUFF USED IN THE COMBINE_S_P
generations = floor(log2(cells));  %which generation are the cell in
markerPositive = arrayfun(@(c)any(tree.(markerField)(tree.cellNr==c)==1),cells); % which cells are marker postive


% put all necc info together
precalcStructure.cells = cells;
precalcStructure.leaves = leaves;
precalcStructure.birthtimes = birthtimes;
precalcStructure.deathtimes = deathtimes;
precalcStructure.relevantCells_perEntry = relevantCells_perEntry;
precalcStructure.relevantCells_perEntry_IX = relevantCells_perEntry_IX;
precalcStructure.isRootFlag = isRootFlag;
precalcStructure.generations=generations;
precalcStructure.markerPositive = markerPositive;