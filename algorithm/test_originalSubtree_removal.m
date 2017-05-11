%% im trying to remove the originalSubtree field form the transformed-structure
% as it takes up a lot of space
%
% currently the routines using originalSubtree have been modified to work without it and
% asserts have been inserted to ensure it works as before
% just run the below code and see if any of those asserts blows


nBest = 20;
[x_best, ~, ~] = qOpt_plotOptimizer('unitTestData/qEst_noUndiffRemoval.mat',nBest); 

theta = x_best(1,:);
dataFile = 'unitTestData/PuGata_transformed_ALLgen_noREMOVALofUNDIFF.mat';
data= load(dataFile);

tra = advanced_multiTree_preprocessing(data.thinnedTrees_pruned ,'wavelength_2','absoluteTime');
g_ = @(theta)advanced_LinearDiffProcess_likelihood_observedTree_CME(tra,[theta(1:2)],[theta(3:4) round(theta(5))],@calcDivisionMatrix_noDIVISION,1,1,'wavelength_2', true);
loglike_aug = g_(theta);