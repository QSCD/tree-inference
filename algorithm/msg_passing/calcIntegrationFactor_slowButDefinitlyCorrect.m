
% calculates the integration factor of the graphical model at the root of the subtree
% this is the older, slower, but tested and correct version of "calcIntegrationFactor_fastReplacement.m"
% 
% we just keep it, because its maybe easier to understand (just a foor loop no vectorisation)
% THE FASTER VERSION IS DEFINITLY ALSO CORRECT!
%
% the results are:
% p_delayTrees(i) is the prob to see the subtree, given differentiation happened
% at diffPoints(i)!
function p_delaytrees4 =calcIntegrationFactor_slowButDefinitlyCorrect(mine_dist_end_fpt,subtree,delayProcess_struct,WL)

    maxT = max(subtree.diffPoints);
    p_delaytrees4 = zeros(1,length(subtree.diffPoints));
    
%     error('This here does not work if cells dont have makrer onset')
    for i = 1:length(p_delaytrees4)

        tmp_time = maxT-subtree.diffPoints(i); %how much is the time between the end_dist and the diff timepoint
        % careful, if its only a branch/single cell we have to work with the fpt
        if length(subtree.allCells)==1
            
            %use the FPT
            
            time_ix = find(delayProcess_struct.timescale==tmp_time);
            transMat= delayProcess_struct.transitionDensities(:,:,time_ix);
            hasOnset = subtree.(['hasOnset_' WL])(subtree.allCells==1);

            if hasOnset %whats the prob to see the onset at the observed timepoint
                assert(sum(mine_dist_end_fpt)==1,'the end of the leaves must be a delta dist')
                if time_ix==1 %stuid exception at instant onset
                    p_delaytrees4(i)=transMat(1,end);
                else
                    transMat_before= delayProcess_struct.transitionDensities(:,:,time_ix-1);
                    p_delaytrees4(i) = transMat(1,end)-transMat_before(1,end);
                end                
                
            else 
               %if we didnt see an onset, we calc the probability
               %that the onset didnt occur until now
               % THATS just 1-prob(x->target) (the prob to not end up in the absorbing state)
               p_delaytrees4(i)= 1-transMat(1,end)';  
            end

        else %just propagate the distribution
            transMat= delayProcess_struct.transitionDensities(:,:,delayProcess_struct.timescale==tmp_time);
            tmp = transMat*mine_dist_end_fpt';
            p_delaytrees4(i) = tmp(1);
        end
    end
end