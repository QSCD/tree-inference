% looks for the leaves of the tree (or potentially a subtree, if specified) and returns the cellNrs of the leaves in
% the tree
function leaves= tUtil_getLeaves(mTree,subtreeRoot)


if exist('subtreeRoot','var')
    assert(length(subtreeRoot)==1)
    nodes = subtreeRoot;
else
    
%else figure out the root of the tree (doesnt have to be the 1, but often it is)
% if its not the 1, make sure its realy a single conntexted tree 
    cells = unique(mTree.cellNr);
    
    %in case of empty trees, just return empty list
    if isempty(cells)
        leaves = [];
        return
    end
    
    
    oldestAncestor = zeros(1,length(cells));%the earliest ancestro of that cell that ecists in the tree
    for i = 1:length(cells)
        rp = rootpath(cells(i));
        existingAncestors = rp(ismember(rp,cells));
        oldestAncestor(i) = min(existingAncestors);
    end
    assert(all(oldestAncestor==oldestAncestor(1)),'treee seems to be disconnected');

    nodes = oldestAncestor(1);
    
    if nodes~=1
        warning(['tree root is not cell 1 but ' num2str(nodes)])
    end
end

assert(any(mTree.cellNr==nodes)) %make sure we have the "root"

leaves = [];
while ~isempty(nodes)
    nextlevel = [];
    for i =1:length(nodes)
        currentNode = nodes(i);
        % ist i die mutter zweier blätter?
        lD = 2.*currentNode;
        rD = 2.*currentNode+1;
        if any(mTree.cellNr==lD) | any(mTree.cellNr==rD)
            if any(mTree.cellNr==lD)
                nextlevel = [nextlevel lD];
            end
            if any(mTree.cellNr==rD)
                nextlevel = [nextlevel rD];
            end
        else
            leaves = [leaves currentNode];
        end
    end
    nodes = nextlevel;
    
end


end