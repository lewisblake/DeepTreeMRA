function [] = create_plots(data, predictionLocations, predictionMean, predictionVariance, verbose, displayPlots, savePlots, plotsFilePath)
%% CREATE_PLOTS creates the plots for the 'prediction' calculationType when plotting is true
%
%   Input: data, predictionLocations, predictionMean, predictionVariance,
% verbose, displayPlots, savePlots, plotsFilePath
%
%   Output: [] (empty vector)
%

if verbose
    disp('Beginning plotting')
end
% get the size of the data array
[~, nCols] = size(data);
% values is stored in different column depending if linear trend removed
if nCols == 3
    values = data(:,3);
elseif nCols == 4
    values = data(:,4);
end

% Figure 1
figure
if ~displayPlots
    set(gcf,'Visible', 'off');
end
scatter(data(:,1), data(:,2), 5, values,'square', 'filled');
colormap(parula)
colorbar
[cmin,cmax] = caxis;
caxis([cmin, cmax])
title('Observations')
xlabel('x coordinates')
ylabel('y coordinates')
if savePlots
    saveas(gcf, fullfile(plotsFilePath, 'observations'), 'png');
end

% Figure 2
figure
if ~displayPlots
    set(gcf,'Visible', 'off');
end
scatter(predictionLocations(:,1), predictionLocations(:,2), 5, predictionMean, 'square', 'filled');
colormap(parula)
colorbar
caxis([cmin, cmax])
title('Prediction mean')
xlabel('x coordinates')
ylabel('y coordinates')
if savePlots
    saveas(gcf, fullfile(plotsFilePath, 'prediction_mean'), 'png');
end

% Figure 3
figure
if ~displayPlots
    set(gcf,'Visible', 'off');
end
scatter(predictionLocations(:,1), predictionLocations(:,2), 5, predictionVariance, 'square', 'filled');
colormap(flip(autumn))
colorbar
title('Prediction variance')
xlabel('x coordinates')
ylabel('y coordinates')
if savePlots
    saveas(gcf, fullfile(plotsFilePath, 'prediction_variance'), 'png');
end

if verbose % Display progress indicators
    disp('Plotting completed')
end
end

