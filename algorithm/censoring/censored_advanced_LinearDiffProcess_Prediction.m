% given the parameters, for each tree, what are the most likely  diff cells and timepoints
% this is somehow the maximum likelihood hidden tree!
% in fact it creates a new tree that has a field 'trueState', indicating
% undiff (0), diff(1) or marker (2)
%
% --------------------------------------------------------------------------------------
% with the advanced/faster/recursive way of calculating the likelihood
% we do not evaluate the likelihood of each hidden tree seperately (avoiding combinatorial explosion)
%
% this makes it more difficult to predict the most likely hidden tree.
% in fact, we can cleverly enumerate the possibilities recursively and 
% use a early stopping criterion (where we know that all recursion from now on will yield 
% only more improbable hidden trees than the current best hidden tree)
% typically saving a lot of recursion steps
% -----------------------------------------------------------------------------------------
%
% Parameters
%   - trees : the usual array of tree structures
%   - theta : a 1x5 vector containing the estimated parameters
%   - FATEMARKER : fieldname in trees, indicating which marker onset we're interested in
%   - timeframe : only needed for comparison with the old method
%
% Output
%   - mostLikelyHiddenTrees: array of tree structure identical to input trees, but added a field
%     trueState indicating undiff (0), diff(1) or marker (2)
%   - mostlikelyCells: array containing for each tree the cells the most likely differentiated (the cells where trueState switched 0->1)
%   - mostlikelyCells: same with timepoints (relative timepoints in the correspoding cell: TP=1 means first timepoint of the cell)
%   - levelOfConfidence: (DOESNT WORK YET) a value to indicate how unique this maximum likelihood
%     estimate is. its the fold difference between MLEprediction and next best
%     prediction 
%     ranges from Inf (we're very sure about the prediction) to 1 (we have no clue, the MLE and the second are equal)
%     this level of confidence is also added as a attribute to the tree
function [mostLikelyHiddenTrees, mostlikelyCells, mostLikelyTPs,  levelOfConfidence, marginals_per_cell]= censored_advanced_LinearDiffProcess_Prediction(trees,theta,FATEMARKER,timeframe,transformed)


error('not yet implemented to deal with the additional possbility/hidden tree of censoring!')

assert(size(theta,1)==1);
assert(size(theta,2) ==6, 'censored likelihood needs another parameter')

mostLikelyHiddenTrees = cell(1,length(trees));
mostlikelyCells = cell(1,length(trees));
mostLikelyTPs  = cell(1,length(trees));
levelOfConfidence = zeros(1,length(trees));

marginals_per_cell = cell(1,length(trees));
% we calculate for each tree on the fly this transformed variabvles which allow to apply the
% likelihood

%do the calculations
if ~exist('transformed','var') | isempty(transformed)
    transformed=  advanced_multiTree_preprocessing(trees,'wavelength_2',timeframe);
end

g = @(theta)censored_advanced_LinearDiffProcess_likelihood_observedTree_CME(transformed,[theta(1:3)],[theta(4:5) round(theta(6))],@calcDivisionMatrix_noDIVISION,1,1);
[~, D_structs, prob_DelayTrees_Time] = g(theta);

textprogressbar('Predicting trees ')
for i = 1:length(trees)
    textprogressbar(i*100/length(trees))
    cT = trees{i};
    theCells = sort(unique(cT.cellNr));
    assert(all(theCells==sort(D_structs{i}.cells)))
    %% only use the n-best
    nBest = 50; %some more to get the marginals right
        
    %new version with log space
    [log_tmp_decCells_top5, log_p_best_top5] = log_enumerate_hidden_trees_nBest(cT, D_structs{i}.logD_array, D_structs{i}.cells,D_structs{i}.PPrime,nBest,FATEMARKER);    
    [log_pbestHT, ix_best] = max(log_p_best_top5);
    decCells_best = log_tmp_decCells_top5{ix_best};    
    
    %figure out most likely TPs
    [dummy,subtree_indexing] = ismember(decCells_best,transformed.transformedTrees{i}.subtreeRoots);assert(all(dummy))
    possible_decidingTPs = prob_DelayTrees_Time{i}(subtree_indexing);
    [~,MLE_TPs] = cellfun(@(t)max(t),possible_decidingTPs);    
    
    %build the tree
    coloredTree=colorPredTrees_withTime(cT, decCells_best,MLE_TPs,FATEMARKER); 
    % figure;tUtil_drawMulitpleTrees(coloredTree,'trueState','timeframe','absoluteTime','treeLabel','tttFile','EXTEND_LEAVES_FLAG',true)

    mostLikelyHiddenTrees{i}  = coloredTree;
    mostlikelyCells{i} = decCells_best;
    mostLikelyTPs{i} = MLE_TPs;
    
    if length(log_p_best_top5) == 1
        levelOfConfidence(i) = Inf;
    else
        levelOfConfidence(i) = exp(log_p_best_top5(1)-log_p_best_top5(2));
    end
    
    %get a kind of marginal across cells
    nBest = 100; %some more to get the marginals right
    [log_tmp_decCells_top100, log_p_best_top100] = log_enumerate_hidden_trees_nBest(cT, D_structs{i}.logD_array, D_structs{i}.cells,D_structs{i}.PPrime,nBest,FATEMARKER);    

    %normalize it
    log_p_best_top100 = log_p_best_top100- logsumexp(log_p_best_top100);
    predCells = unique([log_tmp_decCells_top100{:}]); % all cells that poped up as differentiated in the topN hidden trees
    logmarginals = zeros(1,length(predCells));
    for j = 1:length(predCells)
        cC = predCells(j);
        [ix] = cellfun(@(c)any(c==cC),log_tmp_decCells_top100); % which hidden trees contain this one
        logprob = logsumexp(log_p_best_top100(ix)); % all these Hts contributed their probability to cell cC
        %wheres cC in the predCells array
        ix_CC =  predCells==cC;
        assert(logmarginals(ix_CC)==0)
        logmarginals(ix_CC) = logprob;
    end
    marginalTmp.cells = predCells;
    marginalTmp.marginal_prob = exp(logmarginals);
    marginals_per_cell{i} = marginalTmp;
end
textprogressbar(' Done')
end

function overall =  logsumexp(Xi)
    % calculcates log[ exp(X1) +...+ exp(Xn)]
    maxX = max(Xi);
    overall = maxX + log( sum(exp(Xi-maxX)) );
end