function [corMatrix] = plot_corMatrix(data, labels, title, fontSize, FilePath, FileName)
    % Function that calculates the correlation matrix of an inputted data
    % array and plots the correlation matrix as a heatmap
    %   data = data array that is used to create the correlation matrix
    %   labels = x and y axis labels for the heatmap (corresponds to the
    %   entered data)
    %   title = heatmap plot title
    %   fontSize = font size for axis labels
    %   FilePath, FileName = file path and name for saving the plot

%% Correlation matrix
corMatrix = corrcoef(data);

% creating empty upper, right triangle
write_ones = ones(size(corMatrix));
lowerTri = tril(write_ones);
corMatrix(~lowerTri) = nan;
removeOne = (corMatrix(:,:)==1);
corMatrix(removeOne)=nan;

%% Plotting correlation matrix heatmap

% heatmap limits
% maxCor_kValues = max(max(abs(corMatrix(:,:))));
maxCor_kValues = 0.75;

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 40 30]); %x_width=10cm y_width=15cm

h = heatmap(corMatrix, 'ColorLimits', [-maxCor_kValues*1.01 maxCor_kValues*1.01], ...
    'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ");
h.Colormap=redblueColormap(256);

h.XDisplayLabels = labels;
h.YDisplayLabels = labels;
h.Title = title;
h.FontSize = fontSize;

%% Save plot

prompt = 'Save figure? Y/N: ';
SaveFile = input(prompt,'s');

if SaveFile == 'Y'
    figfile1 = fullfile(FilePath, FileName);
    saveas(h, figfile1, 'epsc' );
else
end

end