% plot the resulting distributions of the delay estimation and compares to
% the true ones

% estimate_array: a cell, each element containing the best estimates of
%   some procedure (eg trees vs breanches)
% 

function plotResults_ofEstimation(estimate_array,colorArray,trueParameters,dt)
assert(size(trueParameters,1)==1)
%howe to plot the true curves
trueLineStyle = 'k--';
trueLineWidth = 3;

t_tmp = 0:dt:(140.*3600);
    %% compare the diff distributions
figure
subplot(3,1,1)
hold on

%iterate over all procedures to estimate
for k = 1:length(estimate_array)

    current_x_best = estimate_array{k};
    cColor = colorArray{k};
    
    %iterate over all estimates
    for i =1:size(current_x_best,1)
       plot(t_tmp./3600,P_from_lambda_linearLambda(current_x_best(i,1),current_x_best(i,2),t_tmp),cColor)
    end
end
%plot the true shape
plot(t_tmp./3600,P_from_lambda_linearLambda(trueParameters(1),trueParameters(2),t_tmp),trueLineStyle,'LineWidth',trueLineWidth)
% xlabel('Time (h)');
ylabel('Prob.')
title('Differentiation')
xlim([min(t_tmp) max(t_tmp)]./3600)
ylim([0 Inf])

%% compare the marker delay distribution
subplot(3,1,2)
hold on
%iterate over all procedures to estimate
for k = 1:length(estimate_array)

    current_x_best = estimate_array{k};
    cColor = colorArray{k};
    %iterate over all estimates
    for i =1:size(current_x_best,1)
       plot(t_tmp./3600,getFPT_distribution_oneStage(current_x_best(i,3:5),t_tmp) ,cColor)
    end
end
%plot the true shape
plot(t_tmp./3600,getFPT_distribution_oneStage([trueParameters(3:5)],t_tmp),trueLineStyle,'LineWidth',trueLineWidth)
% xlabel('Time (h)');
ylabel('Prob.')
title('Delay')
xlim([min(t_tmp) max(t_tmp)]./3600)
ylim([0 Inf])
%% compare the convolution of both
subplot(3,1,3)
hold on
%iterate over all procedures to estimate
for k = 1:length(estimate_array)

    current_x_best = estimate_array{k};
    cColor = colorArray{k};
    %iterate over all estimates
    for i =1:size(current_x_best,1)
        pd = P_from_lambda_linearLambda(current_x_best(i,1),current_x_best(i,2),t_tmp);
        pm = getFPT_distribution_oneStage(current_x_best(i,3:5),t_tmp);
        
        conV = conv(pd,pm);
        newT = 0:dt:(max(t_tmp).*2);
%         [newT conV] = convolve_to_observedMarker_contTime(pd,pm,t_tmp);
        plot(newT./3600,conV ,cColor)
    end
end
%plot the true shape
pd = P_from_lambda_linearLambda(trueParameters(1),trueParameters(2),t_tmp);
pm = getFPT_distribution_oneStage(trueParameters(3:5),t_tmp);
% [newT conV] = convolve_to_observedMarker_contTime(pd,pm,t_tmp);
conV = conv(pd,pm);

plot(newT./3600,conV ,trueLineStyle,'LineWidth',trueLineWidth)
xlabel('Time (h)');
ylabel('Prob.')
title('Convolution')
xlim([min(t_tmp) max(t_tmp)]./3600)
ylim([0 Inf])


end