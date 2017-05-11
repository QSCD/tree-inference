% creates a structure that helps to efifciently evaluate the likelihood of a tree (sum over all posibilities)
% without this huge sum over all enumerations
% 
% tree          : the tree structure itself
% S             : precalculated array of probabilites of delayTrees (one per cell in the tree)
% PDiff_Prime   : probability for a particular cell (and its upstream cells) to be undifferentiated (just a product of lambdas), but
%                 in such a way that it does not count ancestors twice for related cells ("takes only left sisters")
%                  for a better explanation, see my notes!
% P_cells       : the cells corresponding to the entries in pdiffprime
% 
% PPrime_rawHazardperCell: the raw hazards of a cell differenitating (needed to correct pprime sometimes)
%
% startingFactor: if the cell is not at the start of thew movie (because cutting) theres a certain prob that the cell
%                  even made it that far.
%
% D: called kappa in the paper!
%  
function [probOfObservedTree, D,combinationString, log_probOfObservedTree, logD]= combine_S_and_PDiff_Prime(precalcStructure,S,S_cells,PDiff_Prime,P_cells,PPrime_rawHazardperCell)

leaves = precalcStructure.leaves;
cells = precalcStructure.cells;
generations = precalcStructure.generations;
markerPositive = precalcStructure.markerPositive; % which cells are marker postive

assert(all(S_cells==P_cells),'make sure those two s_cell and p_cell have the same ordering')
% starting from the leaves, create the D
D = NaN(1,length(cells));

logD = NaN(1,length(cells)); %trying to deal with numerics here

%starting from the leaeves go upwards
combinationString = cell(1,length(cells)); %to visually keep track of the mathj thats going on
for g= max(generations):-1:0

    cellsOfGen = cells(generations==g);
    % iteraeting over all cells int the given generation
    for i = 1:length(cellsOfGen);
        cC=cellsOfGen(i);
        ixS = find(S_cells==cC);  % find the cells in the two arrays
        ixP = find(P_cells==cC);        
        
        pprimeCorrected = PDiff_Prime(ixP)./PPrime_rawHazardperCell(ixP); % prob to observed ancestors of the current cell undiff. (pprime is the prob of ancestors undiff * currentCell undiff, hence the division)
        if PDiff_Prime(ixP) == 0 && PPrime_rawHazardperCell(ixP)==0
            pprimeCorrected = 0;
        end
        
        if any(leaves==cC) % if the cell is a leave
            %if it has a marker onset, its the usual formula
            if markerPositive(cells==cC)
                D(ixS) = S(ixS).* pprimeCorrected; %prob of the subtree (S) multiplied by the prob of upstream(mother and ancestors) undif
%                 combinationString{ixS} = sprintf('S[%.2d]*C[%.2d]',cC,cC);
                logD(ixS) = log(S(ixS)) + log(pprimeCorrected);
            else
                D(ixS) = S(ixS).* pprimeCorrected + PDiff_Prime(ixP); %either starts a subtree here(1term), or it doesnt (2nd term, just the prob to be undiff)
%                 combinationString{ixS} = sprintf('S[%.2d]*C[%.2d] + p[%.2d]',cC,cC,cC);

                logD(ixS) = log(S(ixS).* pprimeCorrected + PDiff_Prime(ixP));
            end
        else %otherwise its a function of its children
            childrenIx1 = P_cells == cC.*2;
            childrenIx2 = P_cells == cC.*2+1;
            
            %in case a child is missing just omitt by< setting to one
            if all(childrenIx1==0)
                dL = 1;
                logDL = 0;
            else
                dL = D(childrenIx1);
                logDL = logD(childrenIx1);
            end
            
            if all(childrenIx2==0)
                dR = 1;
                logDR = 0;
            else
                dR = D(childrenIx2);
                logDR = logD(childrenIx2);
            end            
            %neat, we no longer need a separate case for the root, it self corrects via pprimeCorrected
            % (before, we had to explicitly mulitply the root prob by the startingFactor)
            D(ixS) = S(ixS).* pprimeCorrected + dL*dR;
            combinationString{ixS} = sprintf('S[%.2d]*C[%.2d] + (%s) * (%s)',cC,cC,combinationString{childrenIx1},combinationString{childrenIx2});
            
            % calc via logsumexp
            X1 = log(S(ixS)) + log(pprimeCorrected); 
            X2 = logDL + logDR;
            logD(ixS) = logsumexp([X1, X2]);
        end
        
%         if isnan(logD(ixS))
%            'd' 
%         end
    end
end
% the probability of observing any hidden tree that gave rise to the observed tree is now contained
% in the entry of the root
probOfObservedTree = D(1); 
log_probOfObservedTree = logD(1);
end


function overall =  logsumexp(Xi)
% calculcates log[ exp(X1) +...+ exp(Xn)]

if all(Xi==-Inf)  %edge case, as factoring out -Inf doesnt work: -Inf- (-Inf) = NAN
    overall = -Inf;
else
    maxX = max(Xi);
    overall = maxX + log( sum(exp(Xi-maxX)) );
end
end
