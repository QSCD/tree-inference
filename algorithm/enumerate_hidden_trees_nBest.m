% function to figure out the most likely hidden tree given this D_array
% which is used to calculate the likelihood of the sum over all hidden trees
%
%   - D_array: 1 x n vector, where n is the number of cells in the tree
%   - cells:   1 x n vector; cells(i) is the cellnumber correspoding to D_array(i)
%   - nBest:    find the nBest combinations
%
% NOTE:
% better use the log-version (log_enumerate...) for numeric stability!!!

function [prunings,bestP] = enumerate_hidden_trees_nBest(tree,D_array,logD_array,cells,pprime,nBest,MARKER)

currentCell = 1;
marginals = D2marginal(D_array,cells);  % S(i) in the paper

globalStruct.D = D_array;
globalStruct.logD = logD_array;
globalStruct.S = marginals;
globalStruct.cells = cells;
globalStruct.PPrime = pprime;
globalStruct.positive = ismember(cells,unique(tree.cellNr(tree.(MARKER)==1)));

bestP = zeros(1,nBest);
log_bestP = -Inf(1,nBest);
[prunings,bestP,logbestP_dummy] = enum_helper(currentCell,bestP,log_bestP,globalStruct,nBest);

%do a final sorting which is missed sometimes if the S(i) case is best
[bestP,sortix] = sort(bestP,'descend');
prunings = prunings(sortix);
end

%hlper for recursion
% D_array,marginals, cells are just global variables needed in each recursion
%
function [pruning, bestP, logbestP] = enum_helper(currentCellNr,bestP_stack,logbestP_stack, globalVarsStruct,nBest)

    cells = globalVarsStruct.cells;
    currentS = globalVarsStruct.S(currentCellNr==cells);
    currentPrime = globalVarsStruct.PPrime(currentCellNr==cells);
    
    lD = currentCellNr.*2;
    rD = lD + 1;
    isLeave = ~any(cells == lD) & ~any(cells == rD);
    isMarkerPositive = globalVarsStruct.positive(cells == currentCellNr);
    %this is the D(2i)*D(2i+1) : by definition its D(i) - marginal(i)
    %careful here, depends on wether cell is a leave and if it haas the marker
    % 1. D(i) = mar(i)                  if leave and positive
    % 2.      = mar(i) + PDiff_Prime(i) if leavev and not positive
    % 3.      = mar(i) + D(2i)D(2i+1)   if not leave    
    
    %==============================
    %1. if its a positive leave
    if isLeave & isMarkerPositive
        upperBound = currentS; % because its the only posibility
        if upperBound < min(bestP_stack)
            pruning = {};
            bestP = bestP_stack;
            logbestP = logbestP_stack;
            return 
        else
            pruning = {currentCellNr};
            bestP = currentS;
            logbestP = log(currentS);
            return
        end
    %==============================    
    %2. if its a negative leave: 2 posiblities
    elseif isLeave & ~isMarkerPositive
        % this is actually two terms: diff + undiff alternatves
        first_P = currentS-currentPrime;
        sec_P = currentPrime;
        upperBound = max(first_P,sec_P); 
        if upperBound < min(bestP_stack)
            pruning = {};
            bestP = bestP_stack;
            logbestP = logbestP_stack;
            return 
        else        
            pruning = {currentCellNr, []};
            bestP = [first_P, sec_P];
            logbestP = log(bestP);
            return
        end
        
    %==============================
    % 3. has some children
    elseif ~isLeave & ~isMarkerPositive
        d2i_x_d2iplus1 = globalVarsStruct.D(currentCellNr==cells) - currentS;
        
        x1 = globalVarsStruct.logD(currentCellNr==cells);
        x2 = log(currentS);
        log_d2i_x_d2iplus1 =  LogDiffExp(x1,x2) ;
        
        upperBound = max(currentS,d2i_x_d2iplus1);
        % if it is worse than the worst (best)
        if upperBound < min(bestP_stack)
            pruning = {};
            bestP = bestP_stack;
            logbestP = logbestP_stack;
            return
        else
            % need the recursion
            
            if any(cells == lD)
                [pruning_left, bestP_left, log_bestP_left] = enum_helper(lD,bestP_stack,logbestP_stack, globalVarsStruct,nBest);
            else
                pruning_left = {};
                bestP_left = 0;
                log_bestP_left = -Inf;
            end
            if any(cells == rD)
                [pruning_right, bestP_right, log_bestP_right] = enum_helper(rD,bestP_stack,logbestP_stack, globalVarsStruct,nBest);
            else
                pruning_right = {};
                bestP_right = 0;
                log_bestP_right =-Inf;
            end
            
            % put the results together
            [combinations, bestP, logbestP] = combine_helper(pruning_left,bestP_left, pruning_right, bestP_right, nBest, log_bestP_left, log_bestP_right);
            %     alternativ kann man auch einfach die root prunen, falls noch
            %     keinerder beiden subtrees über das pruning limit hinaus ist!!
            combinations{end+1}= currentCellNr; %(**)
            bestP(end+1) = globalVarsStruct.S(globalVarsStruct.cells==currentCellNr);
            logbestP(end+1) = log(globalVarsStruct.S(globalVarsStruct.cells==currentCellNr));
            pruning = combinations;
        end
    else
        error('not possible')
    end

end


function [combinations, ps, log_ps] = combine_helper(lD_prunings,best_left, rD_prunings, best_right, nBest, log_best_left, log_best_right  )

ASSERTING = false;
    combinations = {};
    
    % only right daughter
    if isempty(lD_prunings)& ~isempty(rD_prunings)
        %only take the nBest
         [best_right, sortIx] = sort(best_right,'descend');
         
         [log_best_right, sortIx2] = sort(log_best_right,'descend');
         if ASSERTING
            assert(all(sortIx2==sortIx))
         end
         combinations = rD_prunings(sortIx);
        if length(combinations) > nBest
            combinations = combinations(1:nBest);
            best_right = best_right(1:nBest);
            log_best_right = log_best_right(1:nBest);
        end
        
        ps = best_right;
        log_ps = log_best_right;
        counter = length(combinations)+1;
        
    elseif isempty(rD_prunings)& ~isempty(lD_prunings)
        [best_left, sortIx] = sort(best_left,'descend');
        
        [log_best_left, sortIx2] = sort(log_best_left,'descend');
        if ASSERTING
            assert(all(sortIx2==sortIx))
        end
        combinations = lD_prunings(sortIx);
        if length(combinations) > nBest
            combinations = combinations(1:nBest);
            best_left = best_left(1:nBest);
            log_best_left = log_best_left(1:nBest);
        end  
        counter = length(combinations)+1;
        ps = best_left;
        log_ps = log_best_left;
    %if both daughters are present, combine them
    elseif ~isempty(rD_prunings)& ~isempty(lD_prunings)
        l1 = length(lD_prunings); l2 = length(rD_prunings); % should actually be <nBest
        combinations = cell(1, l1.*l2);  %prealloc for some speed
        counter = 1;
        prob = zeros(1, l1.*l2);
        logprob = zeros(1, l1.*l2);
        for i = 1:l1
            pr_1_i=lD_prunings{i}; %prealloc for speed;
            prob_left = best_left(i);
            log_prob_left = log_best_left(i);
            for j=1:l2
                c = [pr_1_i rD_prunings{j}];
                prob(counter) = prob_left.* best_right(j);
                logprob(counter) = log_prob_left+  log_best_right(j);
                
                combinations{counter} = c;
                counter = counter+1;
            end
        end

        %only take the best
        [best_combined, sortIx] = sort(prob,'descend');
        [log_best_combined, sortIx2] = sort(logprob,'descend');
        if ASSERTING
            assert(all(sortIx==sortIx2))
        end
        
        combinations = combinations(sortIx);
        ps = best_combined;
        log_ps = log_best_combined;
        if length(combinations)>nBest
            combinations = combinations(1:nBest);
            ps = ps(1:nBest);
            log_ps = log_ps(1:nBest);
        end
    else
        combinations = {};
        ps = [];
        log_ps = [];
%         error('both are empty, should be caught be leave case')
    end 
    if ASSERTING
     assert(all(log(ps)- log_ps < eps));
    end
end


function overall = LogDiffExp(X1,X2)
maxX = max(X1,X2);

overall = maxX + log( exp(X1-maxX) - exp(X2-maxX) );
end
