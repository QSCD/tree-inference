function unitTest___advanced_LinearDiffProcess_Prediction()

%load parameters and data for unitTest
[x_best, ~, ~] = qOpt_plotOptimizer('unitTestData/qEst_noUndiffRemoval.mat',20); 
theta = x_best(1,:);
dataFile = 'unitTestData/PuGata_transformed_ALLgen_noREMOVALofUNDIFF.mat';
data= load(dataFile);
trees = data.thinnedTrees_pruned;

g = @(theta)advanced_LinearDiffProcess_likelihood_observedTree_CME(data.transformed,[theta(1:2)],[theta(3:4) round(theta(5))],@calcDivisionMatrix_noDIVISION,1,1);
[~, D_structs, prob_DelayTrees_Time] = g(theta);

% for each tree compare preditions using the linear and logprob methods
for i = 1:length(data)
    cT = trees{i};
      
    nBest = 1000;

    %olf version with linear space
    [tmp_decCells_top5, p_best_top5] = enumerate_hidden_trees_nBest(cT,D_structs{i}.D_array,D_structs{i}.logD_array, D_structs{i}.cells,D_structs{i}.PPrime,nBest,'wavelength_2');    
    
    %new version with log space
    [log_tmp_decCells_top5, log_p_best_top5] = log_enumerate_hidden_trees_nBest(cT, D_structs{i}.logD_array, D_structs{i}.cells,D_structs{i}.PPrime,nBest,'wavelength_2');  
    
    assert(all(cellfun(@(a,b) all(a==b), tmp_decCells_top5, log_tmp_decCells_top5)), 'predictions of most likely tree differ')
    assert(all(abs(log(p_best_top5) - log_p_best_top5)<1e-14),'probs and logprobs different')
    
    % test D2marginal and its log version
    log_marginals = log_D2marginal(D_structs{i}.logD_array,D_structs{i}.cells);
    marginals = D2marginal(D_structs{i}.D_array,D_structs{i}.cells);
    assert(all(abs(exp(log_marginals)-marginals)<eps),'D2marginal and log_D2marignal are wrong')
end
disp('passed')
end