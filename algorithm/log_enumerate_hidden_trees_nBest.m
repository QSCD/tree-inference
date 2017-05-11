% function to figure out the most likely hidden tree given this D_array
% which is used to calculate the likelihood of the sum over all hidden trees
% 
% works in logspace !!
% => see enumerate_hidden_trees_nBest(). 
%
%   - D_array: 1 x n vector, where n is the number of cells in the tree
%   - cells:   1 x n vector; cells(i) is the cellnumber correspoding to D_array(i)
%   - nBest:    find the nBest combinations
%
function [prunings,logbestP] = log_enumerate_hidden_trees_nBest(tree,logD_array,cells,pprime,nBest,MARKER)

currentCell = 1;
logmarginals = log_D2marginal(logD_array,cells);

globalStruct.logD = logD_array;
globalStruct.logS = logmarginals;
globalStruct.cells = cells;
globalStruct.PPrime = pprime;
globalStruct.positive = ismember(cells,unique(tree.cellNr(tree.(MARKER)==1)));

log_bestP = -Inf(1,nBest);

[prunings, logbestP] = log_enum_helper(currentCell,log_bestP, globalStruct,nBest);

%do a final sorting which is missed sometimes if the S(i) case is best
[logbestP,sortix] = sort(logbestP,'descend');
prunings = prunings(sortix);
end

% version working in logspace of probabilities
function [pruning, logbestP] = log_enum_helper(currentCellNr,logbestP_stack, globalVarsStruct,nBest)

    cells = globalVarsStruct.cells;
    log_currentS = globalVarsStruct.logS(currentCellNr==cells);
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
        logupperBound = log_currentS; % because its the only posibility
        if logupperBound < min(logbestP_stack)
            pruning = {};
            logbestP = logbestP_stack;
            return 
        else
            pruning = {currentCellNr};
            logbestP = log_currentS;
            return
        end
    %==============================    
    %2. if its a negative leave: 2 posiblities
    elseif isLeave & ~isMarkerPositive
        % this is actually two terms: diff + undiff alternatves
        first_P = exp(log_currentS)-currentPrime;
        if first_P<0
            first_P=0;
            warning('first_P <0');
        end
        
        sec_P = currentPrime;
        logupperBound = max(log(first_P),log(sec_P)); 
        if logupperBound < min(logbestP_stack)
            pruning = {};
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
         
        % we just want to calc globalVarsStruct.D(currentCellNr==cells) -
        % currentS, but in logspacve
        x1 = globalVarsStruct.logD(currentCellNr==cells);
        x2 = log_currentS;
        log_d2i_x_d2iplus1 =  LogDiffExp(x1,x2) ;
        
        logupperBound = max(log_currentS,log_d2i_x_d2iplus1);
        % if it is worse than the worst (best)
        if logupperBound < min(logbestP_stack)
            pruning = {};
            logbestP = logbestP_stack;
            return
        else
            % need the recursion
            
            if any(cells == lD)
                [pruning_left, log_bestP_left] = log_enum_helper(lD,logbestP_stack, globalVarsStruct,nBest);
            else
                pruning_left = {};
                log_bestP_left = -Inf;
            end
            if any(cells == rD)
                [pruning_right, log_bestP_right] = log_enum_helper(rD,logbestP_stack, globalVarsStruct,nBest);
            else
                pruning_right = {};
                log_bestP_right =-Inf;
            end
            
            % put the results together
            [combinations, logbestP] = log_combine_helper(pruning_left,log_bestP_left, pruning_right, log_bestP_right, nBest);
            %     alternativ kann man auch einfach die root prunen, falls noch
            %     keinerder beiden subtrees über das pruning limit hinaus ist!!
            combinations{end+1}= currentCellNr; %(**)
            logbestP(end+1) = log_currentS;
            pruning = combinations;
        end
    else
        error('not possible')
    end

end


% version working in logspace of probabilities
function [combinations, log_ps] = log_combine_helper(lD_prunings,log_best_left, rD_prunings, log_best_right, nBest )


    combinations = {};
    
    % only right daughter
    if isempty(lD_prunings)& ~isempty(rD_prunings)
        %only take the nBest         
         [log_best_right, sortIx] = sort(log_best_right,'descend');

         combinations = rD_prunings(sortIx);
        if length(combinations) > nBest
            combinations = combinations(1:nBest);
            log_best_right = log_best_right(1:nBest);
        end
       
        log_ps = log_best_right;
        counter = length(combinations)+1;
        
    elseif isempty(rD_prunings)& ~isempty(lD_prunings)

        [log_best_left, sortIx] = sort(log_best_left,'descend');

        combinations = lD_prunings(sortIx);
        if length(combinations) > nBest
            combinations = combinations(1:nBest);
            log_best_left = log_best_left(1:nBest);
        end  
        counter = length(combinations)+1;
        log_ps = log_best_left;
    %if both daughters are present, combine them
    elseif ~isempty(rD_prunings)& ~isempty(lD_prunings)
        l1 = length(lD_prunings); l2 = length(rD_prunings); % should actually be <nBest
        combinations = cell(1, l1.*l2);  %prealloc for some speed
        counter = 1;
        logprob = -Inf(1, l1.*l2);
        for i = 1:l1
            pr_1_i=lD_prunings{i}; %prealloc for speed;
            log_prob_left = log_best_left(i);
            for j=1:l2
                c = [pr_1_i rD_prunings{j}];
                logprob(counter) = log_prob_left+  log_best_right(j);
                combinations{counter} = c;
                counter = counter+1;
            end
        end
        
        %only take the best
        [log_best_combined, sortIx] = sort(logprob,'descend');
        
        combinations = combinations(sortIx);
        log_ps = log_best_combined;
        if length(combinations)>nBest
            combinations = combinations(1:nBest);
            log_ps = log_ps(1:nBest);
        end
    else
        combinations = {};
        log_ps = [];
%         error('both are empty, should be caught be leave case')
    end 
end


% logExpSum trick but for two varialbes and difference instead of sum
function overall = LogDiffExp(X1,X2)
maxX = max(X1,X2);

overall = maxX + log( exp(X1-maxX) - exp(X2-maxX) );
end