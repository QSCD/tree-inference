%applies perfect halfing to the density of the statespace (e.g proteins)

function densityAfterDivision =  applyCellDivision(density,statespace)

    assert(length(density)==length(statespace)+1); %+1 for the absoring state
    

 %% new stuff with some uggly indxing, however its fast
% state x(i) gets contributions from:   1* x(i*2) + 0.5 x(i*2 +1)  + 0.5 x(i*2 -1)

% in my indexing: index1 = protein0, hence shift i by j:=i-1
% i = 1:end/2
% 
% proteinAmount = 1:(length(density)-1)./2-1  
oldDens = density;
density = density(1:end-1); %get rid of the sink state, this one is not passed on anyway

proteinAmount = 0:length(density)-1; %the amount of protein in each state
maxProteinAmount = ceil(proteinAmount(end)./2); %soviel hat man maxc nach der halbierung
proteinAmount = 0:maxProteinAmount;

ix_double = proteinAmount*2 +1;
ix_doublePlus1 = proteinAmount*2 -1 +1 ;
ix_doubleMinus1 = proteinAmount*2 +1 +1;

dummyState_ix = length(density)+1;

ix_double(ix_double<1) = dummyState_ix;
ix_doublePlus1(ix_doublePlus1<1) = dummyState_ix;
ix_doubleMinus1(ix_doubleMinus1<1)=dummyState_ix;
ix_double(ix_double>length(density)) = dummyState_ix;
ix_doublePlus1(ix_doublePlus1>length(density)) = dummyState_ix;
ix_doubleMinus1(ix_doubleMinus1>length(density))=dummyState_ix;
density = [density 0];

x = density(ix_double) + 0.5.*density(ix_doublePlus1) + 0.5.*density(ix_doubleMinus1 );
densityAfterDivision = [x zeros(1,length(oldDens)-length(x))];
end

