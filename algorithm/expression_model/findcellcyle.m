% how long are the  cells in the tree
function cellCycles = findcellcyle(tree)

cells = unique(tree.cellNr);

cellCycles = zeros(1,length(cells));
for i = 1:length(cells)
    ix = tree.cellNr==cells(i);
    time = tree.absoluteTime(ix);
    cellCycles(i) = time(end)-time(1);
end

end