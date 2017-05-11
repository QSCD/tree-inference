function precompile_CME_ODEfiles(maximalState,folder)

% for each number of state up two maximalState, create an ode file
% the is the RHS of the CME with those states

stateRange = 1:maximalState;

for i = 1:length(stateRange)

    currentMaxState = stateRange(i);
    modelName = ['oneStageCME_CVODE_StateRange' num2str(0) '_' num2str(currentMaxState) ];

%     if ~exist([modelName '.mexa64'],'file')
        dummy_synth = 1;
        dummy_decay = 1;
        statespace = 0:currentMaxState;
        
        filename = [folder '/oneStageCME_CVODE_precalc_TransitionDensities_CME' num2str(0) '_' num2str(currentMaxState) '.txt'];
        createCVODE_ModelFile_CME_oneStage(dummy_synth,dummy_decay,max(statespace),modelName,filename);
        
        model = SBmodel(filename);
        SBPDmakeMEXmodel(model,[folder '/' modelName]);
        delete(filename)
end