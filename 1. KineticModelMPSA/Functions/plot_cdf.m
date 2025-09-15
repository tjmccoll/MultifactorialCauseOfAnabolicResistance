function [] = plot_cdf(tiledLayout, mpsaValues, mpsaPass, title_indiv, title_figure, fileName, filePath)
    % Function to plot the cumulative distribution function of the inputted
    % data
    %   tiledLayout = 2 values defining the number of tiles in the x and y
    %   axes
    %   mpsaValues = array of parameter values used in the MPSA
    %   mpsaPass = parameter sets that passed the MPSA criteria
    %   title_indiv = titles of individual CDFs
    %   title_figure = title of full figure
    %   fileName, FilePath = used to save figure

%% Identifying successful and unsuccessful parameter sets
passTest_success = mpsaValues(:, mpsaPass==1);
passTest_fail = mpsaValues(:, mpsaPass==0);

%% Creating CDF plot
clf
t = tiledlayout(tiledLayout(1), tiledLayout(2), "TileSpacing","tight");
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 40 30]);

for k = 1:size(mpsaValues,1)
    
    nexttile
    cdfplot(passTest_success(k,:))
    hold on
    cdfplot(passTest_fail(k,:))

    if min(mpsaValues(k,:)) == max(mpsaValues(k,:))
        xlim([-0.5 0.5])
    else
        xlim([min(mpsaValues(k,:)) ...
        max(mpsaValues(k,:))])
    end
%         if k <= 30    
%             xlim([min(concentration_mpsa_limit(k,:)) ...
%                 max(concentration_mpsa_limit(k,:))])
%         elseif k > 30
%             xlim([-0.5 0.5])
%         end
%     xlim([min(mpsaValues(k,:)) ...
%         max(mpsaValues(k,:))])
    xlabel('')
    ylabel('')
    title(title_indiv(k))

    title(t, title_figure)
    xlabel(t, 'Parameter Range')
    ylabel(t, 'Cumulative Frequency')

end
hold off

%% Saving figure
prompt = 'Save figure? Y/N: ';
SaveFile = input(prompt,'s');

if SaveFile == 'Y'
    figfile1 = fullfile(filePath, fileName);
    saveas(t, figfile1, 'epsc' );
else
end

end
