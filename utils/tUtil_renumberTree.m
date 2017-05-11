% for a given input tree,
% change the cellNr of the root towards "newRootnumber" and also change all
% its children correspondingly

function [tree, cellNr_mapping] =tUtil_renumberTree(tree,newRootnumber ,oldRootnumber)
    
    % one can either specify the old root of the tree (fast)
    % or determine it in here (slower)
    if ~exist('oldRootnumber','var')
        %make sure its a genuine tree (no loose ends, all have on common ancestor)
        cells = unique(tree.cellNr);
        oldestAncestor = zeros(1,length(cells));%the earliest ancestro of that cell that ecists in the tree
        for i = 1:length(cells)
            rp = rootpath(cells(i));
            existingAncestors = rp(ismember(rp,cells));
            oldestAncestor(i) = min(existingAncestors);
        end
        assert(all(oldestAncestor==oldestAncestor(1)),'treee seems to be disconnected');
      
        % the root of the tree 
        oldRootnumber = oldestAncestor(1);
    end
    
    % calc the new cellnumbers for the addTree
    cellNr_mapping.oldCellNr = [];
    cellNr_mapping.newCellNr = [];
    nodes =[oldRootnumber];
    newNodes = [newRootnumber];
    
    
    while(~isempty(nodes))
        nextLevel= [];
        nextLevel_new= [];
       for i = 1:length(nodes)
          
           cN = nodes(i);
           
           cN_new = newNodes(i);
           cellNr_mapping.oldCellNr = [cellNr_mapping.oldCellNr cN];
           cellNr_mapping.newCellNr = [cellNr_mapping.newCellNr cN_new];
           
           d1 = cN.*2;
           d2 = cN.*2+1;
           if any(tree.cellNr==d1)
               nextLevel = [nextLevel d1];
               nextLevel_new = [nextLevel_new 2.*cN_new];
           end
           if any(tree.cellNr==d2)
               nextLevel = [nextLevel d2];
               nextLevel_new = [nextLevel_new 2.*cN_new+1];
           end           
       end
       nodes = nextLevel;
       newNodes = nextLevel_new;

    end
    
    %% now renumber the tree
    % now put the new number onto addTree
    newCellNr = zeros(1,length(tree.cellNr));
    for i = 1:length(cellNr_mapping.oldCellNr)
        ix = tree.cellNr==cellNr_mapping.oldCellNr(i);
        
        newCellNr(ix) = repmat(cellNr_mapping.newCellNr(i),1,sum(ix));
    end
    
    tree.cellNr = newCellNr;
end