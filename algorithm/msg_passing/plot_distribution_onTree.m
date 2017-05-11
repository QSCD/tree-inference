function plot_distribution_onTree(cells,birthDist,endDist,xlimits)

assert(length(cells)==length(birthDist) & length(cells) == length(endDist));

generations = floor(log2(cells));
maxGen = max(generations);

hSubplots = 2.^maxGen+1; 
vSubplots = maxGen+1;


for i = 1:length(cells)
    currentCell = cells(i);
    currentGen = generations(i);
    maxCells_inthisGen = 2.^currentGen;
    
    %each cell gets equal space here, depending on how many there are in
    %this gen
    number_of_subplots_occupied = hSubplots./maxCells_inthisGen;
    ix = 1:number_of_subplots_occupied:hSubplots;
    ix2 = number_of_subplots_occupied:number_of_subplots_occupied:hSubplots;
    
    ix_withinGeneration = 1+currentCell-2.^currentGen; 
    
    space_vect = [ix(ix_withinGeneration) ix2(ix_withinGeneration)]; %thats the "subplots" this one will occupy
    
    subplotix = currentGen.*hSubplots + space_vect;
    subplot(vSubplots,hSubplots,subplotix)
    
    assert(size(birthDist{i},1)==1 | size(birthDist{i},2)==1,'the birth distribution is 2dimensianol!!' )
    assert(size(endDist{i},1)==1 | size(endDist{i},2)==1,'the birth distribution is 2dimensianol!!' )
    
    plot(birthDist{i});
    hold on
    plot(endDist{i},'g');
    xlim(xlimits)
end