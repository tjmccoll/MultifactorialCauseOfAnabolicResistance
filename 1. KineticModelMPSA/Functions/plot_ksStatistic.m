function plot_ksStatistic(mpsa_parameter, passTest, pValue_sig, figNum, plotTitle, xLabels, figTitle, filePath)
    % Evaulation of the K-S statistic and plotting of the K-S statistic for
    % each parameter with p-values included for the statistically
    % significant parameters
    %   mpsa_parameter = a list of the randomly generated parameters used
    %   in the MPSA (i.e., a list of the k-values or concentrations for the
    %   n generated parameter sets for the MPSA trials)
    %   passTest = a list of the accepted or rejected parameter sets from
    %   the MPSA (accepted = 1, rejected = 0)
    %   pValue_sig = pre-defined level of statistical significance (e.g.,
    %   p<0.05)

    %% Evaulating the K-S statistic

    % Identifying the accepted or rejected parameter sets
    passTest_success = mpsa_parameter(:, passTest==1);
    passTest_fail = mpsa_parameter(:, passTest==0);

    % Evaulating the K-S statistic
    for k = 1:length(mpsa_parameter(:,1)) %number of parameters
        [h(k), pValue(k,1), ksStat(k)] = kstest2(passTest_success(k,:), ...
            passTest_fail(k,:), 'Alpha', pValue_sig);
    end
    
    %% Plotting the K-S statistic

    % K-S statistic with the highest pValue that is still statistically
    % significant
    highestSigValue = min(ksStat(pValue<pValue_sig & min(pValue)));

    figure(figNum)
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 30 15]); %x_width=10cm y_width=15cm
%     set(gca, 'XTickLabel', 'fontsize', 6)

    hy = bar(ksStat);
        ax=gca;
        ax.FontSize = 6;    
%     xlabel('Parameter', 'FontSize', 10);
    ylabel('K-S Statistic', 'FontSize', 10);
    title(append('Model acceptance: ', plotTitle), 'FontSize', 14)
    ylim([0 0.5]) % update
%     xticks([1,5:5:length(hy(1).XData)]) % update
    xticks(1:1:length(hy(1).XData))
    xticklabels(xLabels)

    % adding p-value text above bar if statistically significant
    for i=1:size(hy,2)
        X = hy(i).XData;
        Y = hy(i).YData;
        for j = 1:size(X,2)
            if pValue(j,i)<pValue_sig && pValue(j,i)>=0.0001
                text(X(1,j)-(2-i)*0.3*hy(i).BarWidth,Y(1,j)+.01, ...
                    num2str(pValue(j,i), 'p = %.4f'), ...
                    'Rotation',90);
            elseif pValue(j,i)<0.0001                
                text(X(1,j)-(2-i)*0.3*hy(i).BarWidth,Y(1,j)+.01, ...
                    'p < 0.0001', ...
                    'Rotation',90);
            end
        end
    end

    % adds a red, dashed horizontal line to delineate statistical 
    % significance
    hold on
    plot(xlim,[.98*highestSigValue .98*highestSigValue], '--r', 'LineWidth', 2)
    text(0, 1.15*highestSigValue, sprintf('p = %.2f', pValue_sig), 'Color','r', 'FontSize',8)

    %% Save plot
    
    fileName = append(figTitle, sprintf(' - %s', datestr(now, 'yy-mm-dd HH-MM-SS')));
    
    prompt = 'Save figure? Y/N: ';
    SaveFile = input(prompt,'s');

    if SaveFile == 'Y'
        figure1 = figure(figNum);
        figfile1 = fullfile(filePath, fileName);
        saveas(figure1, figfile1, 'epsc' );
    else
    end
end