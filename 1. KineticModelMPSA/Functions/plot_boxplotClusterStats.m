function [] = plot_boxplotClusterStats(clusterParameters, fullParameterSet, boxplotTitles, fileName, filePath)


%%
nClusters = size(clusterParameters,3);

t = tiledlayout("flow", "TileSpacing","tight");
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 40 30]);

for j = 1:size(clusterParameters,1) 
    for i = 1:nClusters
        boxplot_clusterInput(:,i) = clusterParameters(j,:,i); % cols = cluster, rows = data
    end

    nexttile
    boxchart(boxplot_clusterInput, 'MarkerStyle','.', 'JitterOutliers','on');
    title(boxplotTitles(j))
    if min(fullParameterSet(j,:)) == max(fullParameterSet(j,:))
        ylim([-0.5 0.5])
    else
        ylim([min(fullParameterSet(j,:)) max(fullParameterSet(j,:))])
    end 
%     ylim([min(fullParameterSet(j,:)) max(fullParameterSet(j,:))])
    set(gca,'XTickLabel', strcat('c',string(1:nClusters)))
    hold on
    plot(mean(rmmissing(boxplot_clusterInput)),'.')
    hold on
end

%% Saving figure
prompt = 'Save figure? Y/N: ';
SaveFile = input(prompt,'s');

if SaveFile == 'Y'
    figfile1 = fullfile(filePath, fileName);
    saveas(t, figfile1, 'epsc' );
else
end

end
