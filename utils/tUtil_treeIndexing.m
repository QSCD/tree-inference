% applies the indexing to all fields of the tree, except those that are not
% real tree datapoint but e.g. the treeId, position etc

function thinnedTree = tUtil_treeIndexing(mTree,ix)

singletonFields = {'movieID','tttFile','quantHash'}; %fields that dont have to be subindexed anyways

thinnedTree = struct();
fields_of_struct = fieldnames(mTree);
for j = 1:length(fields_of_struct)
    currentfield = fields_of_struct{j};
    %manche felder gibts nicht für jede zelle, nur einmal kopieren
    if size(mTree.(currentfield),2)~= length(mTree.cellNr) | any(strcmp(singletonFields,currentfield))
        thinnedTree.(currentfield) = mTree.(currentfield);
    else
        %bei 1d attributes, its easy
        if size(mTree.(currentfield),1)==1 |size(mTree.(currentfield),2)==1
            thinnedTree.(currentfield) = mTree.(currentfield)(ix);
        else %for 2d we have to know where to index
            assert(size(mTree.(currentfield),2)==length(mTree.cellNr),'problem indexing a matrix, time must be along the columns!!')
            thinnedTree.(currentfield) = mTree.(currentfield)(:,ix);
        end
    end
end

end