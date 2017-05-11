% creates and solves the ode for the CME of a birth death process on the
% statespace
% PRECOMPILED : flag that tells wether to use precompiled mexFiles (must be in the path)
% or compile on the fly
% note that if we want to access a precompiled file that doesnt exist, that
% whole thing will crash
function processStruct= precalc_TransitionDensities_CME(statespace,timescale,theta,PRECOMPILED)

assert(timescale(1)==0,'first element of the timescale must be 0!!')
assert(~(any(theta<0) | any(isinf(theta)) | any(isnan(theta))),'parameters funny')
synth =theta(1);
decay =theta(2);

%create the A matrix used for Highams method as well as for simple matrix exp
extended_statespace = [statespace statespace(end)+1];
diagonalEntries = diag(-(synth+extended_statespace.*decay));
upperDiag = diag(extended_statespace(2:end).*decay,1);
lowerDiag = diag(synth(:, ones(length(extended_statespace)-1,1)),-1);
A = diagonalEntries+upperDiag+lowerDiag;
A(end,end) = 0;
A(end-1,end) = 0;


% MODE = 'ode';
MODE = 'expm';
% what mode is used to calculate the transition density
switch MODE
    case 'ode'
    %my old way, using the ode45 integrator
        [transitionDensities, CVODE_OPTIONS] = integrate_via_odeSolver(statespace,timescale,synth,decay,PRECOMPILED);
        processStruct.CVODE_OPTIONS = CVODE_OPTIONS; % additionally, remember the ode parameters
        
    case 'higham'
        % testing the method by Higham;
        %
        % first tests showed that its 4fold slower actually. i presume that the matlab code
        % of expmv is slow, because it includes a loop thats heavily used!!!
        % would be nice to utilize the C implementation
        %
        %     addpath('/home/mjr/PhD/projects/general/code/MATLAB/util/higham_matrixExp')   
        integrate_via_higham(A,statespace,timescale)
        
    case 'expm'
        % other idea: we need the whole matrix anyway:
        % above we calculate expm(At)* delta(x,x0) for all x0 in statespace
        % this could be done in a single operation by expm(At) * eye(): each column of eye corresponds to a
        % different delta function
        % hence it would be sufficient to just calc expm(At), and the resulting resulting matrix is just all
        % the solutions stacked to each other by coluumn
        %
        % furthermore we just need expm(A * dT) because all other timesteps can be obtained by matrix
        % multiplication: expm(A.* 2dT) = expm(A*dT)^2
        
        transitionDensities = integrate_via_expm(A,statespace,timescale);
        %check accuracy
        %         for i = 1:length(statespace)+1
        %             largestDiff = max(max(abs(squeeze(X2(i,:,:))-squeeze(transitionDensities(i,:,:)))));
        %             disp(largestDiff)
        %         end
    otherwise
        error('unkown mode. use ode, higham or expm as MODE')
end


processStruct.transitionDensities = transitionDensities;
processStruct.timescale = timescale;
processStruct.statespace = statespace;
processStruct.absorbingState_ix = size(transitionDensities,1); %which row/col correspondsto the absorbinf state
processStruct.borderStates = statespace(end); % the states that have transitions  into the absorbing state;

processStruct.parameters = [synth decay];
end

% =========================================================================================================
function [transitionDensities, CVODE_OPTIONS]= integrate_via_odeSolver(statespace,timescale,synth,decay,PRECOMPILED)
    %*******
    % CREATE THE ODE system using CVODE
    %*******
    filename = 'oneStageCME_CVODE_precalc_TransitionDensities_CME.txt';
    modelName = ['oneStageCME_CVODE_StateRange' num2str(statespace(1)) '_' num2str(statespace(end)) ];

    if PRECOMPILED==0
        %if its not compiled already, do it
        if ~exist([modelName '.mexa64'],'file')
            createCVODE_ModelFile_CME_oneStage(synth,decay,max(statespace),modelName,filename);
            model = SBmodel(filename);
            SBPDmakeMEXmodel(model)
        else
           warning('you choose not to precompile, but a precompiled file that matches the specifications is in the path') 
        end
    else
        % if PRECOMPILED==1 we assume that theres a mex-file in the path with is
        % named accorindg to the variable modelName
        % this file is then called in the loop below
        assert(exist([modelName '.mexa64'],'file')==3,'didnt find a precompile file')
    end

    % 1dim : initial
    % 2dim : target
    % 3dim : time
    transitionDensities=zeros(length(statespace)+1,length(statespace)+1,length(timescale));

    %integration error tolerance
    % CVODE_OPTIONS.abstol = 1e-6;
    % CVODE_OPTIONS.reltol = 1e-6;
    CVODE_OPTIONS.abstol = 1e-12;
    CVODE_OPTIONS.reltol = 1e-12;

    for i = 1:length(statespace)+1
        % starting in one single state
        ini = zeros(1,length(statespace)+1); %the absorbing state is also included
        ini(i)=1;

        %simdataMEX = oneStageCME_CVODE_StateRange0_50(timescale,ini,[synth decay]); %warning this function name must be equal to the content of modelName-varaible
        simdataMEX = feval(modelName,timescale,ini,[synth decay],CVODE_OPTIONS);  
        tmp_States = simdataMEX.statevalues';

        %sometimes the ode solver slips into negative
        tmp_States(tmp_States<0)=0;
        transitionDensities(i,:,:) = tmp_States;
    end
end

% =========================================================================================================
function integrate_via_higham(A,statespace,timescale)
    %create some intputs to the Higham method
    t0 = min(timescale);
    tmax = max(timescale);
    q = length(timescale)-1; %otherwise timescale is different

    %calc the solutions
    for i = 1:length(statespace)+1    

        % starting in one single state
        ini = zeros(1,length(statespace)+1); %the absorbing state is also included
        ini(i)=1;

        shiftFLAG = true;
        forcesimpleFLAG= true; %somehow works best when this is true. otherwise completely fails
    %                  expmv_tspan(A,b,  t0,tmax,q, prec   ,M, shift,    force_simple,   bal,prnt)
        [X,tvals,mv] = expmv_tspan(A,ini,t0,tmax,q,'double',[],shiftFLAG,forcesimpleFLAG);

        %compare accuracy to ODE solution
        largestDiff = max(max(abs(X-squeeze(transitionDensities(i,:,:)))));
        disp(largestDiff)
    end
end
% =========================================================================================================

function transitionDensities = integrate_via_expm(A,statespace,timescale)
    %simply uses the expm to get the integrator for dt and then use matrix multiplication
    dt = diff(timescale); assert(all(dt==dt(1)));
    dt = dt(1);
    expAt = expm(A.*dt);

    X2 = zeros(length(statespace)+1,length(statespace)+1,length(timescale)); % save all the transition densities over time
    X2(:,:,1) = eye(length(statespace)+1);
    %iterate over all times and multiple the prev state with the propgator (expAt) for dT
    for i = 2:length(timescale)
        X2(:,:,i) = X2(:,:,i-1)*expAt;
    end
    X2 = permute(X2,[2 1 3]); % for some reason each X2(:,:,i) has to be transposed to match the 'transitionDensities'
    
    transitionDensities = X2;

end

% =========================================================================================================

function sbtoolbox_nocompile_vERYSLOW()
% dont compile; at the moment this uses the matlab ode solvers
    
    createCVODE_ModelFile_CME_oneStage(synth,decay,max(statespace),modelName,filename);
    model = SBmodel(filename);
    % set the parameters, otherwise it uses the default ones
    model = SBparameters(model,{'k1','k2'},theta);
    
    
    % 1dim : initial
    % 2dim : target
    % 3dim : time
    transitionDensities=zeros(length(statespace)+1,length(statespace)+1,length(timescale));

    for i = 1:length(statespace)+1

        % starting in one single state
        ini = zeros(1,length(statespace)+1); %the absorbing state is also included
        ini(i)=1;

        simdata = SBsimulate(model,'ode45',timescale,ini);
        transitionDensities(i,:,:) = simdata.statevalues';
    end    

end
% =========================================================================================================

function solve_withoutCompile(statespace,timescale,theta)

% ---------
% CVODES initialization
% ---------

%%% time parameters
tstep=timescale(1);
tend=timescale(end);

rtol = 1.0e-5;
atol = 1.0e-3;

options = CVodeSetOptions('RelTol',rtol,...
                          'AbsTol',atol,...
                          'LinearSolver','GMRES',...
                          'UserData',data);

                      
CVodeInit(@rhs, 'BDF', 'Newton', 0, U0, options);

% ----------------
% Problem solution
% ----------------

[status,t,v,vS] = CVode((tstep:tstep:tend),'Normal');

% -----------
% Free memory
% -----------

CVodeFree;    

end