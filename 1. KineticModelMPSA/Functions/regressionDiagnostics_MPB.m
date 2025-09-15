function regressionDiagnostics_MPB(model, file_modelName, FilePath)
% Plotting 5 diagnostic plots:
%     1) Residuals vs fitted
%     2) Fitted vs observed
%     3) Residuals histogram
%     4) Leverage plot
%     5) Cook's distance

%% Residuals vs fitted values
figure(1);
plotResiduals(model, 'fitted', ...
    'Marker','o','Color','black', 'MarkerSize', 5);
grid on;
title('Residuals vs Fitted values');
xlabel('Fitted Values');
ylabel('Residuals');
% ylim([-0.2, 0.8]);

%% Fitted vs Observed values
figure(2);
plot(model.Fitted, model.Variables.MPB, ...
    'ko', 'MarkerSize',5);
grid on;
hold on;
plot(model.Fitted, model.Fitted, 'r--', 'LineWidth',1); % line of perfect fit
title('Fitted vs Observed');
xlabel('Fitted Values');
ylabel('Observed Values');
    obsMin = min(model.Variables.MPB); fitMin = min(model.Fitted);
    obsMax = max(model.Variables.MPB); fitMax = max(model.Fitted);
%     obsMin = min(model.Variables.MPB_transform); fitMin = min(model.Fitted);
%     obsMax = max(model.Variables.MPB_transform); fitMax = max(model.Fitted);
ylim([min(obsMin, fitMin)*0.9, max(obsMax, fitMax)*1.05]);
xlim([min(obsMin, fitMin)*0.9, max(obsMax, fitMax)*1.05]);
% ylim([min(model.Variables.MPS)*0.85, max(model.Variables.MPS)*1.15])
% xlim([min(model.Variables.MPS)*0.85, max(model.Variables.MPS)*1.15])
hold off;

%% Residuals histogram
figure(3);
plotResiduals(model, 'histogram');
grid on;
title('Residuals Histogram');

%% Leverage plot
figure(4);
plotDiagnostics(model, 'leverage');
grid on;
title('Leverage plot');

%% Cook's Distance
figure(5);
plotDiagnostics(model, "cookd");
grid on;
title('Cook''s Distance Plot');

%% QQ plot
figure(6)
qqplot(model.Residuals.Raw)
grid on
title('Q-Q Plot of Regression Residuals')
xlabel('Theoretical Quantiles')
ylabel('Sample Quantiles')

%% Save plot

% prompt = 'Save figure? Y/N: ';
% SaveFile = input(prompt,'s');

% if SaveFile == 'Y'

    figure1 = figure(1);
    fileName_complete = sprintf(append(file_modelName, '_ResidualsVsFitted_%s'), ...
        datestr(now, 'yy-mm-dd HH-MM-SS'));
    figfile1 = fullfile(FilePath, fileName_complete);
    saveas(figure1, figfile1, 'pdf');
    
    figure2 = figure(2);
    fileName_complete = sprintf(append(file_modelName, '_FittedVsObersved_%s'), ...
        datestr(now, 'yy-mm-dd HH-MM-SS'));
    figfile2 = fullfile(FilePath, fileName_complete);
    saveas(figure2, figfile2, 'pdf');

    figure3 = figure(3);
    fileName_complete = sprintf(append(file_modelName, '_ResidualsHistogram_%s'), ...
        datestr(now, 'yy-mm-dd HH-MM-SS'));
    figfile3 = fullfile(FilePath, fileName_complete);
    saveas(figure3, figfile3, 'pdf');

    figure4 = figure(4);
    fileName_complete = sprintf(append(file_modelName, '_LeveragePlot_%s'), datestr(now, 'yy-mm-dd HH-MM-SS'));
    figfile4 = fullfile(FilePath, fileName_complete);
    saveas(figure4, figfile4, 'pdf');

    figure5 = figure(5);
    fileName_complete = sprintf(append(file_modelName, '_CooksDistance_%s'), datestr(now, 'yy-mm-dd HH-MM-SS'));
    figfile5 = fullfile(FilePath, fileName_complete);
    saveas(figure5, figfile5, 'pdf');

    figure6 = figure(6);
    fileName_complete = sprintf(append(file_modelName, '_qqPlot_%s'), datestr(now, 'yy-mm-dd HH-MM-SS'));
    figfile6 = fullfile(FilePath, fileName_complete);
    saveas(figure6, figfile6, 'pdf');

% else
% end

end