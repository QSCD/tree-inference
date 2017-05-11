%uses the backwards approach to calc the probability of a delay tree (calculating actually the prob for all possible timepoints in the root)
% MODE Sspecifies how the prob is evaluated: 
% MODE = 0 -> at leaves, it just propagates backwards the delta peak
% MODE = 1 -> it uses the FPT at leaves
% MODE = 2 -> we dont exactly know where it turned on, because the
%             timepoints before the anotated onset are just not quantified
%             therefore we have to integrate over all timepoint in the
%             leaves from between the insptected timepoint in the beginning
%             to the insptectd tps in the end

function [p_delaytrees] =  backwards_precalculatedCME_calcProb_Single_delayTrees_CME...
                                                (subtree,diffProcess_struct,...
                                                 delayProcess_struct,divisionMatrix_delayProcess,MODE,WL)
   
assert(length(subtree)==1);
DEBUG =0;

if MODE == 1 % USE THE FPT
    [mine_dist_birth_fpt, mine_dist_end_fpt, allCells_fpt] =probability_of_delay_tree_backwards(subtree,delayProcess_struct,divisionMatrix_delayProcess,MODE,WL);

    p_delaytrees =   calcIntegrationFactor_fastReplacement(mine_dist_end_fpt{1},subtree,delayProcess_struct,WL);
    
    if DEBUG
        p_delaytrees_slowerButRight =calcIntegrationFactor_slowButDefinitlyCorrect(mine_dist_end_fpt{1},subtree,delayProcess_struct,WL);
        assert(all(abs(p_delaytrees-p_delaytrees_slowerButRight)<eps))
    end

elseif MODE==2 %we dont know the exact onset timepoint, so in the leaves we integrate over all onsets
    [mine_dist_birth_integrated, mine_dist_end_fpt_integrated, allCells] =probability_of_delay_tree_backwards(subtree,       delayProcess_struct,divisionMatrix_delayProcess,MODE,WL);
    %figure;plot_distribution_onTree(allCells,mine_dist_birth_integrated,mine_dist_end_fpt_integrated,[0 30])
    error('dont use, speedup as for mode =1 not implemented')
    %now consider all possible times of differentiation in the root
    maxT = max(subtree.diffPoints);
    p_delaytrees5 = zeros(1,length(subtree.diffPoints));
    for i = 1:length(p_delaytrees5)

        tmp_time = maxT-subtree.diffPoints(i); %how much is the time between the end_dist and the diff timepoint
        % careful, if its only a branch/single cell we have to work with
        % the fpt, and under MODE 2 the integrated version
        if length(mine_dist_end_fpt_integrated)==1
            
            hasOnset = subtree.(['hasOnset_' WL])(subtree.allCells==1);
            if hasOnset
                assert(sum(mine_dist_end_fpt_integrated{1})==1,'the end of the leaves must be a delta dist') 

                ix_t = find(delayProcess_struct.timescale==tmp_time);
                stateIX = delayProcess_struct.absorbingState_ix;

                if ix_t ==1 % instant onset
                    p_delaytrees5(i) = delayProcess_struct.transitionDensities(1,stateIX,ix_t);
                else

                    %now we dont exactly know when the onset occured we
                    %just now it happend at or before the last timepoint
        %                     tmp_fpt = delayProcess_struct.transitionDensities(1,stateIX,2:ix_t);
        %                     tmp_fpt2 = delayProcess_struct.transitionDensities(1,stateIX,1:ix_t-1);
        %                     p_delaytrees5(i)= sum(tmp_fpt-tmp_fpt2,3); % the first passage time is the "time deriv", therefore the difference. afterwards we sum/integrate over all timepoint!
        %                     error('this is redundant sum and diff = identiy. ')
                    p_delaytrees5(i)=delayProcess_struct.transitionDensities(1,stateIX,ix_t);
                    
                end
            else
                % the tree is a single cell and it idint have any marker
                % onset
                %whats the prob to have no onset after tmp_time, statring
                %from state 1
                ix_t = find(delayProcess_struct.timescale==tmp_time);
                p_delaytrees5(i)=  1-delayProcess_struct.transitionDensities(1,delayProcess_struct.absorbingState_ix,ix_t);
                %its just 1- the prob to reach the absorbing state (which would see as marker)
            end

        else %just propagate the distribution
            transMat= delayProcess_struct.transitionDensities(:,:,delayProcess_struct.timescale==tmp_time);
            tmp = transMat*mine_dist_end_fpt_integrated{1}';
            p_delaytrees5(i) = tmp(1);
        end
    end
    p_delaytrees=p_delaytrees5;
else
   error('mode not 0 or 1') 
end
end



