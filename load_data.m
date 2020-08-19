function [ data, regressionModel, domainBoundaries, predictionVector, theta, varEps ] = load_data(dataSource, nXGrid, nYGrid, offsetPercentage)
%% LOAD_DATA Loads data from various sources
%   Data files are loaded as a function of a variable called dataType
%
%   Input: dataSource
%
%   Output: data, regressionModel, domainBoundaries, predictionVector, theta, varEps
%%
switch dataSource
    
    case 'satellite'
        %% User Input
        % Change as needed
        load('./Data/satelliteData.mat')
        
        % Values of parameters of covariance function
        theta = [5.57,0.12]; varEps = 0.01;
        
    case 'simulated'
        %% User Input
        % Change as needed
        load('./Data/simulatedData.mat')
        
        % Values of parameters of covariance function
        theta = [8.13,0.72]; varEps = 0.1;
        
    case 'amsrSST'
        load('./Data/amsrDay.mat')
        theta = [2.117,1]; 
        varEps = 0.001;  
        
        %lon = lon(1:2:end);
        %lat = lat(1:2:end);
        %obs = obs(1:2:end);
        
    case 'modisSST'
        load('./Data/modisDay.mat')
        theta = [2.117,1];
        varEps = 0.001;
        
        lon = lon(1:3:end);
        lat = lat(1:3:end);
        obs = obs(1:3:end);
    otherwise
        error('Error. Specified datType is not a valid data set.');
end
disp('Loading data complete');

% Determine the boundaries of the domain spanded by the data.
xmin0 = min(lon);
xmax0 = max(lon);
ymin0 = min(lat);
ymax0 = max(lat);
domainBoundaries = [xmin0, xmax0, ymin0, ymax0];

% Make prediction grid
if nXGrid && nYGrid > 0 % If user defines a prediction grid
xPredictionVec = linspace(xmin0 + offsetPercentage, xmax0 - offsetPercentage, nXGrid);
yPredictionVec = linspace(ymin0 + offsetPercentage, ymax0 - offsetPercentage, nYGrid);
[xPredGridLocs, yPredGridLocs] = meshgrid(xPredictionVec, yPredictionVec);
else 
xPredGridLocs = [];
yPredGridLocs = [];
end


% Find observation locations.
logicalInd = ~isnan(obs);

% Declare predicition grid
predictionVector = [xPredGridLocs(:),yPredGridLocs(:)];

% Assign lon, lat and observations to data matrix.
data(:,1) = lon(logicalInd);
data(:,2) = lat(logicalInd);
data(:,4) = obs(logicalInd);

% Detrend data.
regressionModel = fitlm(data(:,1:2),data(:,4), 'linear');
residuals = table2array(regressionModel.Residuals(:, 1));
data(:,3) = residuals;
end
