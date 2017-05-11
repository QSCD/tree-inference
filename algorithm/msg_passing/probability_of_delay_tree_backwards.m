%return the distribution of the cells at birth and death, estimated by
%assuming that the leaves are in the absorbing state


% MODE Sspecifies how the prob is evaluated: 
% MODE = 0 -> at leaves, it just propagates backwards the delta peak
% MODE = 1 -> it uses the FPT at leaves
% MODE = 2 -> we dont exactly know where it turned on, because the
%             timepoints before the anotated onset are just not quantified
%             therefore we have to integrate over all timepoint in the
%             leaves from between the insptected timepoint in the beginning
%             to the insptectd tps in the end


% =====================
% CHANGELOG
% 
% 05/03/2013
% * changed the call of "probability_of_delay_tree_backwards" to not only
%   pass  the whole subtree struct instead of only subtree.originalSubtree
%   so that one can use the precalc stuff inside
% 29/04/2013
% * added the only one daughter if clause

function [distribution_at_birth distribution_at_end allCells ] =probability_of_delay_tree_backwards(tree_struct,processStruct,divisionMat,MODE,WL)

    if MODE == 0 
        error('mode = 0 is only for debugging my stuff')
    end
    VERBOSE = 0;

    %lets start a the leaves as here we know in which state we are: the
    %absorbing state of the delay process
    % then go upwards:
    % within a cell we calc the distribution at the birth of the cell given the
    % end state (which is the abosrbing state for the leaves)
    % then we join sister, as we require that both sisters started from the
    % same state (that of the mother)

    statespace = processStruct.statespace;
    absState_ix = processStruct.absorbingState_ix;
    borderState = processStruct.borderStates;assert(length(borderState)==1,'currently we allow for only on state at the border towards the absorbing state!')
% =====================
% CHANGELOG 05/03/2013
%     allCells = unique(tree.cellNr);
%     allCells_generation = floor(log2(allCells));
%     leaves = tUtil_getLeaves(tree);
%     [~ ,leave_ix] = ismember(leaves,allCells);
%........................
    allCells = tree_struct.allCells;
    hasOnset_perCell = tree_struct.(['hasOnset_' WL]);
    lifetime_perCell = tree_struct.lifetime;

    allCells_generation = floor(log2(allCells));
%======================

    unique_gens = sort(unique(allCells_generation),'descend');
    
    %here we memorize the distribution right after its birth
    distribution_at_birth = cell(1,length(allCells));
    distribution_at_end = cell(1,length(allCells));
 
    
    
    % alle generationen von der höchsten aus durchgehen
    for i =1:length(unique_gens)
       currentGen =  unique_gens(i);
       
       %whats the cells in the current Gen
       cells_in_current_gen = allCells(allCells_generation==currentGen);
       
       %for each cell, join the distribution of its children
       for j = 1:length(cells_in_current_gen)
           
           currentCell = cells_in_current_gen(j);
           currentCell_ix = allCells==currentCell; %index for saving the distribution
           %what about its daughters
           lD = currentCell.*2;
           rD = currentCell.*2+1;
           hasleft = any(allCells==lD  );
           hasright = any(allCells==rD  );
           
           %how long did it live
           current_liftetime = lifetime_perCell(allCells==currentCell);

           transMat= processStruct.transitionDensities(:,:,current_liftetime==processStruct.timescale); % the corresponding transition matrix
           
           %join daughters:
           if ~hasleft & ~hasright 
               %doesnt have daughters, its a leave and we know it had a
               %deltaPeak at its end if it has a marker onset
               % or it has a uniform distribution if we dont have an onset
               
               hasOnset = hasOnset_perCell(allCells==currentCell);
               if MODE == 0
                   error('dont use this, it does not take into account the FPT')
                   if VERBOSE,warning('this is not true, we know more: the system has just entered the absolring state and wasnt there before.');end
                   
                   if hasOnset
                       %now, move this delta peak backward in time, and determine where the cell could have started
                       tmp_dist = transMat(:,borderState==statespace); %its just the colum of the matrix
                       distribution_at_birth{currentCell_ix} = tmp_dist';
                   else
                       %no onset so we have no clue where we are in
                       %statespace, except that were not in the targetstate
                       %-> uniform dist P(y)=c, but targetstate = 0
                       % it is P(x) = Sum_y P(x->)P(y)
                       % which is jsut matrix multiplcation with a col vec
                       unifD = [ones(length(statespace),1) 0]./length(statespace);
                       tmp_dist = transMat(:,borderState==statespace)*unifD;
                   end
                   
               elseif MODE ==1  
                   %we're using the FPT for the observed onset
                   if hasOnset
                       ix_t = find(processStruct.timescale==current_liftetime);

                       % stateIX = borderState==statespace; warning('why border and not absorbing')
                       %===================
                       % actually, for the FPT we want the derivate/difference
                       % in time for the ABSOrBING STATE,hence
                       stateIX = absState_ix;

                       if ix_t ==1 % instant onset
                           fptDist = transMat(:,stateIX);
                           warning('instant onset, this might be buggeD!!')
                       else
                           transMat_before= processStruct.transitionDensities(:,:,ix_t-1);
                           absorbing_t = transMat(:,stateIX);
                           absorbing_tMinus1 = transMat_before(:,stateIX);

                           fptDist = absorbing_t-absorbing_tMinus1;
                       end
                       distribution_at_birth{currentCell_ix} = fptDist';
                   else
                       %if we didnt see an onset, we calc the probability
                       %that the onset didnt occur until now
                       % THATS just 1-prob(x->target) (the prob to not end up in the absorbing state)
                       distribution_at_birth{currentCell_ix} = 1-transMat(:,absState_ix)';
%                        error('has not been tested!!')
                   end
                   
               elseif MODE ==2
                   
                   if hasOnset
                       ix_t = find(processStruct.timescale==current_liftetime);
                       stateIX = absState_ix;

                       if ix_t ==1 % instant onset
                           fptDist = transMat(:,stateIX);
                           warning('instant onset, this might be buggeD!!')
                       else

                           %now we dont exactly know when the onset occured we
                           %just now it happend at or before the last timepoint
            %                            tmp_fpt = processStruct.transitionDensities(:,stateIX,2:ix_t);
            %                            tmp_fpt2 = processStruct.transitionDensities(:,stateIX,1:ix_t-1);
            %                            fptDist_old= sum(tmp_fpt-tmp_fpt2,3); % the first passage time is the "time deriv", therefore the difference. afterwards we sum/integrate over all timepoint!
            %                            error('this is redundant sum and diff = identiy. just be careful what happens at the absorbinf state')
            %                            
                           %integrating the fpt is just the same as the
                           %trnasition dens itself
                           fptDist = processStruct.transitionDensities(:,stateIX,ix_t);
                           fptDist(stateIX)=0; %if we are in the absorbing state already, there shouldnt be a leave anyway
                           
                            %this tells u how likely it is to have an onset at
                            %or before "ix_t" as a function of the starting
                            %state
                       end
                       distribution_at_birth{currentCell_ix} = fptDist';
                   else
                       %if we didnt see an onset, we calc the probability
                       %that the onset didnt occur until now
                       % THATS just 1-prob(x->target) (the prob to not end up in the absorbing state)
                       %no need to integrate here!!
                       %actually not true here, as we dontobserve what
                       %happens in between, a cell could go over the
                       %threshold but be lower again at the end
                       distribution_at_birth{currentCell_ix} = 1-transMat(:,absState_ix)';
                    
                   end
               else
 
                  error('unknown mode. must be 0,1,2') 
               end
               
               if hasOnset
                   % the distribution at the end was the delta
                   tmpDelta= zeros(1,size(transMat,1));
                   tmpDelta(borderState==statespace)=1;
                   distribution_at_end{currentCell_ix}  = tmpDelta;
               else
                   distribution_at_end{currentCell_ix}  =  [ones(1,length(statespace)) 0]./length(statespace);
               end
               
           %has only one
           elseif (hasleft & ~hasright)| (~hasleft & hasright)
%                error('my stuff not implemented')

%% added this part on 29.April 2013
                if hasleft
                    daughterNr = lD;
                    daughter_ix = daughterNr == allCells;
                else
                    daughterNr = rD;
                    daughter_ix = daughterNr == allCells;
                end
                lDDist_mine = distribution_at_birth{daughter_ix};
                motherDist_end = lDDist_mine;
               
                mine_Dist=transMat*motherDist_end';
                
               distribution_at_birth{currentCell_ix}=mine_Dist';
               distribution_at_end{currentCell_ix}=motherDist_end;  
%% =========================================               
               
%                 warning(['strange tree, cell ' num2str(currentCell) ' has only one daughter'] );
%                 if hasleft
%                     daughterNr = lD;
%                     daughter_ix = daughterNr == allCells;
%                 else
%                     daughterNr = rD;
%                     daughter_ix = daughterNr == allCells;
%                 end
%                 
%                 % we just push the daughter state to the mother
% if VERBOSE                
%                 warning('no divsion assumed at all!!')
% end                
%                 motherDist_end = distribution_at_birth{daughter_ix};
%                 distribution_at_end{currentCell_ix} = motherDist_end;
%                 
%                 % now we "propagate towards birth of the mother, see the
%                 % case below for details
%                 birth_dist = propagate_to_birth(motherDist_end,transMat,statespace);
%                 distribution_at_birth{currentCell_ix} =birth_dist;

           %has two daughter, join their distrubtions
           elseif hasleft & hasright % for non leaves both MODES are the same
               
               % WARNING: AT THE MOMENT WE ASSUME NO DIVISION AT ALL
               if VERBOSE
                   warning('no divsion assumed at all!!')
               end
               %join both daughters giving the joint, conditionend on the
               %current cell
               lD_ix = lD == allCells;
               rD_ix = rD == allCells;
               
               %mine, for non leavesits just them same, but no
               %normalisation
               lDDist_mine = distribution_at_birth{lD_ix};
               rDDist_mine = distribution_at_birth{rD_ix};
               
% here incorporate the cell division process

               
               motherDist_end = lDDist_mine.*rDDist_mine;
%                figure;plot([lDDist_mine;rDDist_mine;motherDist_end]')
               
               % do the propagation/integration!!
               % = sum_x p(x|x_mother) * currentCell_endDist(x),
               %this is just amtrix mulitply
               %however, more explicitly:
               %     mine_Dist=transMat*diag(end_dist); mine_Dist = sum(mine_Dist,2) %multiply each column i with end_dist(i) then sum over the cols to integrate it out
               %which is the same as:
               mine_Dist=transMat*motherDist_end';
               
               distribution_at_birth{currentCell_ix}=mine_Dist';
               distribution_at_end{currentCell_ix}=motherDist_end;
           else
               error('blub')
           end
           
           %probability conserved within a cell? 
           %however htis doesnt make sense if we look at the leave cells
           %and the FPT, as the fpt as a function of state does not sum to
           %one
           % NOTE: the FPT as a funtion of time does sum to one!!
%            if hasleft | hasright
%                prob_diff = sum(distribution_at_birth{currentCell_ix})-sum(distribution_at_end{currentCell_ix});
%                assert(abs(prob_diff)<1e-5,['some prob got lost in cell ' num2str(currentCell) '. parameteres:' num2str(processStruct.parameters) ' ' num2str(statespace(end))])
%            end

%%%%%%%%%%%%%%%%%55
%%%%%% I GUESS THIS DOES NOT MAKE SENSE AT ALL
%%%%%% its not like propagating a distribution backwards, its just
%%%%%% evaluating the factorisation of a joint distribution, in particular
%%%%%% on step in this algorihm just calculates p_start(x) = sum_y p_end(y)
%%%%%% * transition(y|x) , which doesnt really say mass conservation


       end
    end
end

