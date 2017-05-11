function plotProcessDistributions(trees,timescale,color)
%if we simulated some diff/delay trees (we know where the diff point was), 
%extract the hazards and distributions from the trees


USE_SUBPLOTS =1;

% figure
[lambda_diff lambda_marker] = tUtil_extract_markerDelay_from_ToyData(trees,timescale);
[lambda_overall, ~, ~] = lambdaFromTrees_perTime(trees,'wavelength_2','absoluteTime',timescale);
% [newTimescale convolvedOnset] = convolve_to_observedMarker_contTime(P_from_Lambda(lambda_diff),P_from_Lambda(lambda_marker),timescale);

if ~exist('color','var')
    
    
    plot(timescale./3600,[P_from_Lambda(lambda_diff); P_from_Lambda(lambda_marker); P_from_Lambda(lambda_overall)]);
    % hold on;plot(newTimescale./3600,convolvedOnset,'r--')
else
    if USE_SUBPLOTS
        subplot(3,1,1),hold on        
        plot(timescale./3600,[P_from_Lambda(lambda_diff)],[color]);
        subplot(3,1,2),hold on
        plot(timescale./3600,[P_from_Lambda(lambda_marker)],[color]);
        subplot(3,1,3),hold on
        plot(timescale./3600,[P_from_Lambda(lambda_overall)],[color]);        
    else
        hold on
        plot(timescale./3600,[P_from_Lambda(lambda_diff)],[color '--']);
        plot(timescale./3600,[P_from_Lambda(lambda_marker)],[color 'x-']);
        plot(timescale./3600,[P_from_Lambda(lambda_overall)],[color '-']);
    end
end

legend({'diff','delay','overall'})
vline([12 24 36],'k--')
xlabel('Time (h)')
ylabel('Prob')
xlim([0 max(timescale./3600)])
title('FPT Distributions of the 2 processes in toy data')