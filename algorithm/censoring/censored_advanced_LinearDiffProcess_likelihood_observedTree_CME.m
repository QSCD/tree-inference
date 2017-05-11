% exponential differentiation but molecular model for marker
%
% accounts for left censoring

% [p_obsTree, D_struct_all, prob_DelayTrees_Time]
% TODO: D_Struct is not correct, censoring is missing in there:
function [p_obsTree] = censored_advanced_LinearDiffProcess_likelihood_observedTree_CME(transformedTreesOverall,k_diff,k_delay,calcDivisionMatrix_func,PRECOMPILED,MODE,WL, logProbFLAG)

if ~exist('logProbFLAG','var')
    logProbFLAG = false;
end

if ~exist('WL','var')
    WL  = 'wavelength_2';
end

if (any([k_diff k_delay]<0) | any(isinf([k_diff k_delay])) | any(isnan([k_diff k_delay])))
    warning('parameters funny')
    p_obsTree = 0;
    return
end

assert(length(k_diff)==3,'WE NOW NEED 3 PARAMETERS, offset slope and mass before movie start')


%***********
% Part 0
% across all trees, gather at which points we need to evaluate (solve the CME for) 
% the two process
% then calc the CME solution once for both processes and use the obtained
% matrices in part 1...
timeRange_firstProcess = transformedTreesOverall.timeRange_firstProcess;
timeRange_secondProcess = transformedTreesOverall.timeRange_secondProcess;
transformedTrees =transformedTreesOverall.transformedTrees;

%now precalculate the CME solution at these timepoints
diff_offset = k_diff(1);
diff_slope = k_diff(2);
diff_beforeMovie =k_diff(3); assert(diff_beforeMovie<=1&diff_beforeMovie>=0) %the prob to diff before the movie starts


threshold_second = k_delay(3);
statespace_second = 0:threshold_second;

tic
secondProcess_struct = precalc_TransitionDensities_CME(statespace_second,timeRange_secondProcess,k_delay(1:2),PRECOMPILED);
divisionMatrix_second = calcDivisionMatrix_func(statespace_second);
time_to_solveODE=toc;

% so, now do the evaluation of the tree probabilites for each transformed
% tree seperatly
p_obsTree = zeros(1,length(transformedTrees));
log_p_obsTree = zeros(1,length(transformedTrees));
prob_DelayTrees_Time = cell(1,length(transformedTrees));
D_struct_all= cell(1,length(transformedTrees));
tic
for i =1:length(transformedTrees)
   
    currentTransformed = transformedTrees{i};
   
    %% ************
    % PART ONE
    % the diff prob is a linear function of time p_diff(t) = offset + slope*t
    % , but is not inherited!!
    %
    % for some tree, its just the product over all timepoints:
    % prod [1-p_diff(t)]^#timepoints
    %************
    
    % how to calc the prob that a given cell does not differentiate:
    % * a cell lives from t0 to t1
    % 
    % * the diff_rate is a function of time lambda(t)
    % 
    % * the "number of times differentiation occurs in [t0 t1]]" follows a poisson dist
    %       with parameter mu = \int_t0^t1 lambda(t) dt
    % 
    % * we are interested in "#differentiation in [t0 t1] = 0, which is
    %   P(k=0| mu) = exp(-mu) mu^k /k! = exp(-mu) = exp(- \int_t0^t1 lambda(t) dt)
    % 
    % * as lambda(t) = offset + slope*t its integral is
    %   offset*(t1-t0) + slope/2 * (t1^2 - t0^2)
    %
    % * so for each cell we just need to calculate 
    %   exp[-  (offset*(t1-t0) + slope/2 * (t1^2 - t0^2))]
    
    
    % =============
    % CHANGE 22/10
    % we have a starting time of the tree, that is: the tree already
    % existed previously to the root cell and we must account for the
    % unbserved upper part to not differentiate!
    % -> no diff in [0 starting time]
    if isfield(currentTransformed, 'startingTime')
        movieStart = currentTransformed.commonUpTreeCorrection_start_stop(1); % notice that this is not necces.  movie start (-> massive overlap if the trees com from the same huge tree), but the time until it merges with some other subtree 
        %note that the .commonUpTreeCorrection_start_stop(1) should be pretty close to
        %startingTime!!
        startingTime = currentTransformed.startingTime;
        additionalFactor = exp(-(diff_offset.*(startingTime-movieStart)+ 0.5.*diff_slope .* (startingTime.^2-movieStart.^2)));
    else
        additionalFactor=1;
        startingTime = 0;
    end   
    
    [PPrime,PPrime_cells,PPrime_rawHazardperCell] = calcP_Prime(currentTransformed.pprime_Structure,k_diff,currentTransformed);
    
    
    
    %% **************
    % PART Two: Prob of the delayTrees, already integrate out the exact time of
    % decision within the rootcell of the delaytree
    % result: matrix of subtree probabilites, depending on the state S of the
    % divided mother cell 
    % rows : different initial states
    % cols : the different subtrees
    %**************
    
    %===============
    % CHANGE 08/01/2014
    % actually we can do a lot of calculations outside of the loop over all subtrees up there:
    % we can propagate from the leaves towards the root only once
    [delayTree_probs_faster, delayTree_probs_vsTime_faster] = calc_delayTree_probs_FASTER(currentTransformed,secondProcess_struct,divisionMatrix_second,MODE,WL,diff_offset,diff_slope,startingTime);
    delayTree_probs = delayTree_probs_faster;
    prob_DelayTrees_Time{i} = delayTree_probs_vsTime_faster;
    %=====================

    
    %% ************
    % PART THREE: put diff tree and delay tree together
    % put it together using a more clever scheme which avoids the large sum over undiff_trees
    S_cells = currentTransformed.subtreeRoots;

    [p_faster, D_array,~, logP_faster,logD_array]= combine_S_and_PDiff_Prime(currentTransformed.pprime_Structure,delayTree_probs,S_cells,PPrime,PPrime_cells,PPrime_rawHazardperCell);
    D_struct.D_array = D_array;
    D_struct.logD_array = logD_array;
    D_struct.cells = PPrime_cells;
    D_struct.PPrime = PPrime;
    D_struct_all{i} = D_struct;

    
    %% the cesoring kicks in here,  there's an additional scenario, the whole tree is jsut delay and even censored at the start
    fullDelayTreeIx = currentTransformed.subtreeRoots==1; 
    
    fullDelayTree = currentTransformed.subtrees{fullDelayTreeIx};
    [distribution_at_start, distribution_at_end, allCells ] =probability_of_delay_tree_backwards(fullDelayTree,secondProcess_struct,divisionMatrix_second,MODE,WL);
    
    % we dont know the state at the beginning of the movie, just put
    % uniform prior on it, equivalent to just averaging over the states
    p_censored = mean(distribution_at_start{1}) ; 
    log_p_censored = log(mean(distribution_at_start{1})) ;
    
    % include the censoring possibility, downweight the other for not being censored:
    p_faster = p_faster * (1-diff_beforeMovie) + p_censored .* diff_beforeMovie;
    
    phi1 = logP_faster + log(1-diff_beforeMovie);
    phi2 = log_p_censored + log(diff_beforeMovie);
    logP_faster = logsumexp([phi1 phi2    ]);
    
    %% ************
    % done for this tree
    p_obsTree(i) = p_faster;
    
    log_p_obsTree(i) = logP_faster;

end
time_2_iterateOverTrees =toc;

% some workaround for p =0
% just replace them all by the smallest nonzero values
smallest = 1e-323;

p_obsTree(p_obsTree==0)=smallest;
p_obsTree(p_obsTree<0)=smallest;

if logProbFLAG %return  log probs if requested
    p_obsTree = log_p_obsTree;
    fprintf(1,'ODE Time: %.5f   TreeLoopTime: %.5f    Parameters [%.5f %.5f %.5f %.5f %.5f %.5f] log10like: %.5f \n',time_to_solveODE,time_2_iterateOverTrees, log(k_diff(1)) ,log(k_diff(2)) ,log(k_diff(3)) ,log(k_delay(1)),log(k_delay(2)),log(threshold_second),sum(p_obsTree))
else
    fprintf(1,'ODE Time: %.5f   TreeLoopTime: %.5f    Parameters [%.5f %.5f %.5f %.5f %.5f %.5f] log10like: %.5f \n',time_to_solveODE,time_2_iterateOverTrees, log(k_diff(1)) ,log(k_diff(2)) ,log(k_diff(3)) ,log(k_delay(1)),log(k_delay(2)),log(threshold_second),sum(log10(p_obsTree)))

end
end

%================
% calculate the probabily L(d|theta,eta), Eq.9 of the paper for each subtree d in the observed tree
%   this basically does the inference on the graph, and integrates to root of the subtree with the
%   differentiation process
%
% FASTER VERSION OF  calc_delayTree_probs(): does the inference only once, not for each d
function [delayTree_probs, delayTree_probs_vsTime] = calc_delayTree_probs_FASTER(currentTransformed,secondProcess_struct,divisionMatrix_second,MODE,WL,diff_offset,diff_slope,startingTime)

    subtrees = currentTransformed.subtrees;
    [~, distribution_at_end, allCells ] =probability_of_delay_tree_backwards(subtrees{1},secondProcess_struct,divisionMatrix_second,MODE,WL);
    
    %iterate over all cells/possible subtrees: 
    % 1. calc their probabilities, considering different length of the root cell
    % 2. fill up the remainder of the root cell with the undiff process (called prefactor prev)
    delayTree_probs = zeros(1,length(subtrees));
    delayTree_probs_vsTime = cell(1,length(subtrees)); % to later tell at which timepoitn the diff happend, we memorize for each delaytree its prob vs timepoint of diff
    for j = 1:length(subtrees)
        assert(allCells(j) == currentTransformed.subtreeRoots(j),'something is wrong in the ordering of the subtrees. this will srew up the results down here, which rely on the ordering of subtreeroots and allCells/distribution_at_end to be the same')

        currentST = currentTransformed.subtrees{j};
        %so far we only got half of the inference on the subtree, now we need to include the
        %integration factor at the root of the graph
        p_delaytrees =   calcIntegrationFactor_fastReplacement(distribution_at_end{j},currentST,secondProcess_struct,WL);
        prefactors =  calcPrefactors(diff_offset,diff_slope,currentST,startingTime);
        delayTree_probs(j) = sum(p_delaytrees.*prefactors);
        delayTree_probs_vsTime{j} =  p_delaytrees.*prefactors;
    end
end    

function overall =  logsumexp(Xi)
% calculcates log[ exp(X1) +...+ exp(Xn)]
    maxX = max(Xi);
    if maxX == -Inf % stupid exception: if all Xi are -Inf we subtract Inf-Inf -> NaN where instead it should be log(0)
        overall = -Inf;
    else
        overall = maxX + log( sum(exp(Xi-maxX)) );
    end
end

%=========================================    

% calcualtes the prefactor as a function of the state
% however in a memoryless model we dont have to worry about the state, so
% it just retuens one number per tree
%
% last argument is the starting time of the whole tree (not only the subtree!)
% needed if we want to calc the hazrard on uncentered time
function prefactors =  calcPrefactors(diff_offset,diff_slope,subtree,startingTime_ofWholeTree)

%=====================0
% CHANGE 22/10/2013
% account for the root of  the whole tree not starting at 0 but some other
% time!! therefore  " + startingTime_ofWholeTree"

    rootBirth = subtree.first_absoluteTime + startingTime_ofWholeTree;

    % the varoius point where diff could happen.
    % note that subtree.diffPoints is relative to the birth of the root
    % cell but diffPoints is in absolute time!!
    diffPoints = subtree.diffPoints+rootBirth;
    
    %stupid exception if the differentiation tree has only one timepoint
    % the prob to immediatly diff is 0!
    if length(diffPoints)==1
        prefactors =0;
        
    else % the normal case
        
        dT = diffPoints(2)-diffPoints(1);
        
        % we look for the prob that the diff happend in a interval [t*-dt, t*]
        % this prob is just the difference of the CDFs: CDF(t*)-CDF(t*-dt)
        % the CDF for an inhomogenious process is just 
        % 1-Exp[-mu], with mu being the integrated rate (Exp(-mu) is the prob of nothing happening): 
        % mu = \int_t0^t1 lambda(t)dt
        % as lambda(t) = a+b*t we can calc the integral

        %CDF t0 -> t*
        mu1 = diff_offset.*(diffPoints-rootBirth) + 0.5.*diff_slope.* (diffPoints.^2-rootBirth.^2);
        CDFtStar = 1-exp(-mu1);

        %CDF t0 -> t*-dt
        mu2 = diff_offset.*(diffPoints-dT-rootBirth)+ 0.5.*diff_slope.*((diffPoints-dT).^2 -rootBirth.^2);
        CDFtStar_Minus_dT = 1-exp(-mu2);

        prefactors = CDFtStar-CDFtStar_Minus_dT;
        % stupid excetion for instant diff at the first timepoint (above stuff gets negative)
        prefactors(1)=0;
    end
end
