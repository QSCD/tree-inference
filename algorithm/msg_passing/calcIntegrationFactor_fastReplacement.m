% calculates the integration factor of the graphical model at the root of the subtree
% this version is a speed up version of "calcIntegrationFactor_slowButDefinitlyCorrect.m"
% THIS FASTER VERSION IS DEFINITLY CORRECT, we just keep the slower one because its easier to understand!
%
% 
% the for i= ... loop ion the slow version can be replaced via vectorisation
% essentially we just do a matrix product in each loop, the matrix changes
% from iter to iter

% the results are:
% p_delayTrees(i) is the prob to see the subtree, given differentiation happened
% at diffPoints(i)!
function p_delayTrees = calcIntegrationFactor_fastReplacement(end_dists,subtree,processStruct,WL)

    assert(~iscell(end_dists) ); %<- enforces using only a single element of the original array
    
    maxT = max(subtree.diffPoints);
    tmp_time = maxT-subtree.diffPoints; %the time that the delay process has to cover 
    
    %figure out the matrix indices belonging to the timepoints
    [a_tmp, time_ix]=ismember(tmp_time,processStruct.timescale);assert(all(a_tmp))
    
    % is it a single cell (else its a treee)    
    if length(subtree.allCells)==1
        
        %does the cell even have an onset?
        hasOnset = subtree.(['hasOnset_' WL])(subtree.allCells==1);
        if hasOnset
            %stupid exception for t=1
            exept = processStruct.transitionDensities(1,end,1);
            
            % ok this one is weird:
            % we do a diff in the time direction of the matrix (the matrix still has forward time: t1,t2,t3,...)
            % however, we have to resort acc to time_ix (which is backward)
            diff_vect = squeeze(diff(processStruct.transitionDensities(1,end,:)));% do the difference in time dimension, and the start and end state
            p_delayTrees = [diff_vect(time_ix(2:end)); exept]';
            
        else
            p_delayTrees =squeeze( 1-processStruct.transitionDensities(1,end,time_ix))';
        end
    else % is it a tree
        
        %note that we dont even need to do full multiplication ABOVE, as we are
        %only interested in the first component!!
        % hence just mulitply the first row of the matrix into
        % end_dists{1}'
        tmp2 = mtimesx(processStruct.transitionDensities(1,:,time_ix),end_dists');
        p_delayTrees = squeeze(tmp2)';
        
        
    end
end