% main_testing.m
kdiff = [0.0001    0.4420].*1e-10;
theta =  [kdiff 0.0005    0.0000001   56.0000];
f(theta)


%% compare the old vs the new version
% based on some of the testSet trees
testSet = load('/home/michi/PhD/projects/MarkerDelay/code_git/MATLAB/linearDiff/toggleSwitch/lalehSwitch/toggleSwitchLaleh_testSet5000.mat');
trees = testSet.trees;
trees = trees(1:100); % take the first 100 trees
tran_old = testSet.transformed;
tran_old.transformedTrees = tran_old.transformedTrees(1:100);

% ======== calc the likelihood the old way!===========
f = @(theta)LinearDiffProcess_likelihood_observedTree_CME(tran_old,[theta(1:2)],[theta(3:4) round(theta(5))],@calcDivisionMatrix_noDIVISION,1,1);
p_old=f(theta);

% == calc the likelihood the new way
tran_new = advanced_multiTree_preprocessing(trees,'wavelength_2','absoluteTime');
g = @(theta)advanced_LinearDiffProcess_likelihood_observedTree_CME(tran_new,[theta(1:2)],[theta(3:4) round(theta(5))],@calcDivisionMatrix_noDIVISION,1,1);
p_new=g(theta);


assert(all(abs(p_old-p_new)<eps),'old and new method give different results!')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USE unitTest___lineaDiff_classic_advanced_likelihood()
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% additionally, test prediction
theta(end) = 10
[mostLikelyHiddenTrees, mostlikelyCells, mostLikelyTPs, levelOfConfidence]= advanced_LinearDiffProcess_Prediction(trees,theta,'wavelength_2','absoluteTime');
