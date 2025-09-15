function [cgo] = plot_clustergram(parameterSet, rowLabels, colLabels, title, dendrogramValue)
    % parameterSet = values being clustered
    % rowLabels = naming of the rows
    % colLabels = naming of the columns
    % title = title of clustergram
    % dendrogramValue = value to distinguish clusters

%%
cgo = clustergram(parameterSet, ...
    'Linkage', 'average', ...
    'Standardize', 'None', ...
    'Colormap', redbluecmap, ...
    'RowLabels', rowLabels, ...
    'ColumnLabels', colLabels, ...
    'Cluster', 'row');
set(0,'ShowHiddenHandles','on')

addXLabel(cgo, 'Parameter set (accepted)')
addYLabel(cgo, 'Parameter number')
addTitle(cgo, title)

% cgo.RowPDist = 'correlation';
cgo.ColumnPDist = 'correlation'; % 

cgo.Dendrogram=dendrogramValue; % 6 clusters of parameter sets

% % reducing label font size
% cgfig = findall(0,'type','figure', 'tag', 'Clustergram'); % Figure handle
% cgax = findobj(cgfig, 'type','axes','tag','HeatMapAxes'); % main axis handle
% % Change fontsize
% % cgax.Fontsize = 8;        
% cgax.FontSize = 2;

end