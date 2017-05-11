function tree=colorPredTrees_withTime(tree, predCells,predTimepoints_ix,markerWL)

    trueState = zeros(size(tree.cellNr));
    for i =1:length(predCells)
        cC =  predCells(i);
        subtree = tUtil_getSubtree(tree,cC  );
        subtreeCells = unique(subtree.cellNr);

        % color everything in these cells to 1
        [ix_sub, ~]= ismember(tree.cellNr,subtreeCells);
        trueState(ix_sub)=1;

        % color the first half of the decided cell to 0 again making it
        % half hjalf

        ix_CC = tree.cellNr==cC;
        time = tree.absoluteTime(ix_CC);  
        diffTP_absolute = time(predTimepoints_ix(i));
        trueState(ix_CC & tree.absoluteTime<diffTP_absolute) =0 ;
    end

    %finally, the timepoints that have marker are colored 2
    trueState(tree.(markerWL)==1)=2;
    tree.trueState = trueState;
end