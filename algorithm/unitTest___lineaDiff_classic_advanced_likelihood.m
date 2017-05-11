function unitTest___lineaDiff_classic_advanced_likelihood()

    unittest_data = load('/home/michi/PhD/projects/MarkerDelay/code_git/MATLAB/linearDiff/algorithm/advanced/unittest_data.mat');
    
    %%perform transformeation and see if it still matches
    
    
    %% perform likelihood evaluation and see if it still matches

    % == calc the likelihood the classical way
    f = @(theta)LinearDiffProcess_likelihood_observedTree_CME(unittest_data.tran_old,[theta(1:2)],[theta(3:4) round(theta(5))],@calcDivisionMatrix_noDIVISION,1,1);
    [p_old, prob_hiddenTrees_old]=f(unittest_data.theta);

    % == calc the likelihood the new way
    g = @(theta)advanced_LinearDiffProcess_likelihood_observedTree_CME(unittest_data.tran_new,[theta(1:2)],[theta(3:4) round(theta(5))],@calcDivisionMatrix_noDIVISION,1,1);
    p_new=g(unittest_data.theta);    
    
    % first make sure they still match one another
    assert(all(p_old==p_new),'classic and addvanced method give diff results')
    
    % also make sureit gives the same results as previously computed when i wrote this test
    assert(all(p_old==unittest_data.p_old),'classic method gives different results than previously')
    
    assert(all(p_new==unittest_data.p_new),'advanced method gives different results than previously')
 
    disp('---------------')
    disp('Passed the test for likelihood')
    disp('---------------')
    
    % =====================================
    % =====================================
    %for the given parameters its boring everthing decides early on
    unittest_data.theta(5) = 10;
    % =====================================
    % =====================================
    
    
    % =====================================
    %test it for a tree with one cell that is undiff
    theTree = unittest_data.trees{1};
    theTree = tUtil_treeIndexing(theTree,theTree.cellNr==1);% | theTree.cellNr==2 | theTree.cellNr==3
    
    theTree_transformed = multiTree_preprocessing({theTree},'wavelength_2','absoluteTime');
    f = @(theta)LinearDiffProcess_likelihood_observedTree_CME(theTree_transformed,[theta(1:2)],[theta(3:4) round(theta(5))],@calcDivisionMatrix_noDIVISION,1,1);
    [p_old, prob_hiddenTrees_old]=f(unittest_data.theta);

    
    % == calc the likelihood the new way
    theTree_transformed_new=  advanced_multiTree_preprocessing({theTree},'wavelength_2','absoluteTime')  ;
    g = @(theta)advanced_LinearDiffProcess_likelihood_observedTree_CME(theTree_transformed_new,[theta(1:2)],[theta(3:4) round(theta(5))],@calcDivisionMatrix_noDIVISION,1,1);
    p_new=g(unittest_data.theta);   
    
    
    [mostLikelyHTrees, mostlikelyCells, mostLikelyTPs, levelOfConfidence]= ...
        advanced_LinearDiffProcess_Prediction({theTree},unittest_data.theta,'wavelength_2','absoluteTime');
    
    % =====================================
    % test for all trees!
    [mostLikelyHTrees, mostlikelyCells, mostLikelyTPs, levelOfConfidence]= ...
        advanced_LinearDiffProcess_Prediction(unittest_data.trees,unittest_data.theta,'wavelength_2','absoluteTime');
    
end

