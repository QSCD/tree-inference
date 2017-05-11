function createCVODE_ModelFile_CME_oneStage(synth,decay,threshold,modelName,filename)

file_1 = fopen(filename,'w');

fprintf(file_1,'********** MODEL NAME\n%s\n',modelName);
fprintf(file_1,'********** MODEL NOTES\n none');

%ODES
fprintf(file_1,'********** MODEL STATES\n');

% if the threshold is 0 its a bit funny as it goe sfrom 0 to sink
% immediatly
if threshold~=0

    fprintf(file_1,'d/dt(P0) = -k1 * P0 + k2 * P1\n'); % ok since it would be -(k1 + 0* k2) P0 + k2 * 1 * P1
     %LOOSING VIA SYNTH, GAIN VIA DEG FROM NEIGHTBOUR 
    for i = 1:threshold-1
        left = i-1;
        right = i+1;
        % ************               LOSS TERM                 GAIN TERMS
        fprintf(file_1,'d/dt(P%d) = -(k1+k2 * %d) * P%d       + k1*P%d + k2 * %d * P%d  \n',i,i,i   ,left,right,right);   
    end
    % we gain from th-1 via synth, butwe dont gain from th+1 (absoribing) via
    % deg
    fprintf(file_1,'d/dt(P%d) = -(k1 + k2 * %d) * P%d    +  k1 * P%d\n',threshold,threshold,threshold,threshold-1);
    fprintf(file_1,'d/dt(PA) = k1 * P%d\n',threshold);

    %sepcial case threshold =0
else
    fprintf(file_1,'d/dt(P0) = -k1 * P0\n'); % ok since it would be -(k1 + 0* k2) P0 + k2 * 1 * P1
    fprintf(file_1,'d/dt(PA) = k1 * P%d\n',threshold);
    
end
%PARAMETERS
fprintf(file_1,'\n********** MODEL PARAMETERS\n');
fprintf(file_1,'k1 = %f\nk2 = %f\n',synth,decay);
fprintf(file_1,'********** MODEL VARIABLES\n');
fprintf(file_1,'********** MODEL REACTIONS\n\n\n********** MODEL FUNCTIONS\n\n\n********** MODEL EVENTS\n\n\n********** MODEL MATLAB FUNCTIONS\n\n\n');
fclose(file_1);

end