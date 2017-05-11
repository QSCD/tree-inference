% loading some simulated trees
load censoredTrees.mat

% doing some precalcualtions on the trees
% wavelength_2 is the marker which we observe
transformed = advanced_multiTree_preprocessing(censoredTrees,'wavelength_2','absoluteTime');


% parameters used to create the data
pCensor = 0.5;
k_diff = [exp([-11.3269  -24.9626]) pCensor];
k_delay = exp([-8.4624  -13.0085    2.8404]);
theta =  [k_diff k_delay];

% evaluate the likelihood on these trees
g = @(theta)censored_advanced_LinearDiffProcess_likelihood_observedTree_CME(transformed,[theta(1:3)],[theta(4:5) round(theta(6))],@calcDivisionMatrix_noDIVISION,1,1);
p_censor=g(theta);


