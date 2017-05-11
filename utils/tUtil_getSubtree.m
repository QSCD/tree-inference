% warning: the root of the subtree is not 1 afterwards but still has the original cellnumber!
function subtree = tUtil_getSubtree(mTree,subtreeroot)

% warning('should add check if tree is connected!')
% checkConnected(mTree) %<- bad runtime

    %check wich nodes are under the subtreeroot
    nodes = subtreeroot;
    subtreeNodes = [subtreeroot];
    
    cells = unique(mTree.cellNr);
    while ~isempty(nodes)
        nextlevel = [];
        for i=1:length(nodes)
            
            ld = nodes(i).*2;
            rd = nodes(i).*2+1;
            
            if any(cells==ld)
                nextlevel = [nextlevel ld];
                subtreeNodes = [subtreeNodes ld];
            end

            if any(cells==rd)
                nextlevel = [nextlevel rd];
                subtreeNodes = [subtreeNodes rd];
            end            
            
        end
        nodes = nextlevel;
    end

    
    % what cells belong to the subtree
    sub_ix = ismember(mTree.cellNr,subtreeNodes);
    
    %do the indexing to get the subtree
    subtree = tUtil_treeIndexing(mTree,sub_ix);
    
end

function checkConnected(tree)

cells = unique(tree.cellNr);
reqCells = [];

%collect the rootpath of all cells
for i = 1:length(cells);
   rp = rootpath(cells(i));
   reqCells = [reqCells rp];
end
reqCells = unique(reqCells);


%if the tree is fully connected, all cells in the rootpaths are in the tree

flag_cellsPresent = ismember(reqCells,cells);
assert(all(flag_cellsPresent),'tree seems to be diconntected, so subtrees wont work proberly')

end