function qJob_advanced_LinearDiff_Delay(jobID)
    VERBOSE = 0;

%     maxNumCompThreads(1);
    initialisePath();
    
    if VERBOSE == 1
        path();
        disp('current WD')
        pwd
    end
    cd(['job',num2str(jobID) '/'])
    
    inputFile = ['inputs.mat'];
    load(inputFile);
    disp('finisihed reading the datafile')

    % load the thing with the transformed trees, but outside the loop, its
    % the same for all iteratuons
    load(inputs{1}.dataFile)
    
    for i =1:length(inputs)
    
        cI = inputs{i};
        
        lb = cI.loglowerBound;
        ub = cI.logupperBound;
        MODE = cI.MODE;
        warning (['using mode = ' num2str(MODE)])
        logprobFLAG = true;
        f = @(logParas,varargin)-sum(advanced_LinearDiffProcess_likelihood_observedTree_CME(transformed,exp(logParas(1:2)),[exp(logParas(3:4)) round(exp(logParas(5)))],@calcDivisionMatrix_noDIVISION,1,MODE, 'wavelength_2',logprobFLAG));
    
        xZeros = cI.xZeros;
        x0 = xZeros;
        
        assert(all(all(xZeros>=repmat(lb,size(xZeros,1),1)) & all(xZeros<=repmat(ub,size(xZeros,1),1))))
        assert(all(lb<=ub),'make sure that the upper bound is bigger than the lower bound')
        assert(size(xZeros,2) == length(lb))
        
        % lets check if the function is real at the point x0
        % if not break (otherwise the optimizer gets stuck here forever)
        obj = f(x0);
        if isinf(obj) | isnan(obj) | imag(obj)~=0
           disp('function not evaluable at x0, terminating')
           continue
        end
        
        x_results = zeros(size(xZeros)); assert(size(xZeros,1)==1,'only use one initial cond')
        fvals = zeros(1,size(xZeros,1));

        if isfield(cI,'optimOptions')
            searchOptions = cI.optimOptions;
        else
           warning('SUPPLY OPTIMSET, using default here') 
           searchOptions = optimset('PlotFcns',{@optimplotfval},'Display','final');
        end
        tic
        
        %WITH threshold optimisation: some cont/int variables
        [x,fval,exFLAG,~] = fminsearchbnd(f, x0, lb, ub, searchOptions);
        %             [x,fval,~, ~] = fmincon(f,x0,[],[],[],[],lb,ub,[], searchOptions); %@optimplotx
        toc
        
        x_results = x;
        fvals = fval;
        exitflag = exFLAG;
        save(['output' num2str(i) '.mat'],'x_results','fvals','xZeros','exitflag')
    end
    
    % stupid queue thing: one has to terminate matlab otherwise the
    % queueJobs gets stuck here
    if isdeployed()
        exit()
    end
end