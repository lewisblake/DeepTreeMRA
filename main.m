function [elapsedTime] = main()
%% Multi-resolution Approximation (MRA) main script
% Executes the MRA model. 
%
% Calculation options include prediction, parameter optimization and likelihood only evaluation
%
% calculationType specifies what is calculated. 
%   Option 'prediction' uses a default values for the parameters and just conducts 
%   the prediction. 
%   Option 'optimize' optimizes over the range, variance and measurement error 
%   and then predicts using the values obtained from the optimization.
%   Option 'likelihood' only calculates the likelihood.
%
% Documentation in these scripts will make references to Katzfuss, 2017.
%% User Input
% Run user_input.m script to get variables into workspace
user_input;

%% Validate User Input
validate_user_input(calculationType, NUM_LEVELS_M, NUM_PARTITIONS_J, NUM_KNOTS_r, offsetPercentage, NUM_WORKERS, nLevelsInSerial, nXGrid, nYGrid, displayPlots, savePlots, verbose, resultsFilePath, plotsFilePath);

%% Data Processing: Load data using load_data() function
[ data, regressionModel, domainBoundaries, predictionVector, theta, varEps ] = load_data(dataSource, nXGrid, nYGrid, offsetPercentage);

%% Build hierarchical grid structure using build_structure_in_parallel() function
[ knots, ~, nRegions, outputData, predictionLocations, indexMatrix ] = build_structure_in_parallel( NUM_LEVELS_M, ...
    NUM_PARTITIONS_J, NUM_KNOTS_r, domainBoundaries, offsetPercentage, NUM_WORKERS, nLevelsInSerial, data(:,1:3), predictionVector );

%% Switch clause for calculationType
switch calculationType
    case 'build_structure'
        %% Build Structure and report summary statistics about observations at finest resolution
        tic;
        compute_structure_statistics(outputData, NUM_WORKERS, resultsFilePath, plotsFilePath);
        elapsedTime = toc;
    case 'optimize'        
        %% Optimize
        isPredicting = false;
        fun = @(thetaOpt)MRA([thetaOpt(1) thetaOpt(2)], outputData, knots, ...
            NUM_LEVELS_M, NUM_PARTITIONS_J, nRegions, indexMatrix, isPredicting, nLevelsInSerial, NUM_WORKERS, verbose, thetaOpt(3));
        % Dummy values required by optimization routine
        A = []; b = []; Aeq = []; beq = [];
        % fmincon() optimizes over the bounds set
        tic; x = fmincon(fun, initalEstimate, A, b, Aeq, beq, lowerBound, upperBound);
        elapsedTime = toc;  % Unsuppress output to print to command window
        % Assign values from optimization to theta and varEps
        theta = [x(1) x(2)];
        varEps = x(3);
        % Save optimize results
        save([resultsFilePath, 'Optimization_Results'], 'theta', 'varEps');
    case 'prediction'
        %% Prediction
        isPredicting = true;
        tic;
        [ ~, predictions ] = MRA(theta, outputData, knots, ...
            NUM_LEVELS_M, NUM_PARTITIONS_J, nRegions, indexMatrix, isPredicting, nLevelsInSerial, NUM_WORKERS, verbose, varEps, predictionLocations);
        elapsedTime = toc;  % Unsurpress output to print to command window
        % Reformat data for plotting:     
        % Collect the distributed predictions, stack them on top of each
        % other, then convert this cell into a matrix. (LB: Stacking should work since regions are sequentially ordered however double check this.
        predictions = cell2mat(reshape(gather(predictions),[],1)); 
        predictionLocations = cell2mat(reshape(gather(predictionLocations),[],1));      
        predictionVariance = predictions(:,2);
        % Add the prediction from the regression
        predRegression = predict(regressionModel, predictionLocations);
        predictionMean = predictions(:,1) + predRegression;
        % Save prediction results
        save([resultsFilePath,'Prediction_Results_', dataSource], 'predictionLocations', 'predictionMean', 'predictionVariance');
        %% Plots
        if displayPlots || savePlots % If plotting
            create_plots(data, predictionLocations, predictionMean, predictionVariance, verbose, displayPlots, savePlots, plotsFilePath)
        end
    case 'likelihood'
        %% Likelihood
        isPredicting = false;
        tic;
        [ sumLogLikelihood ] = MRA(theta, outputData, knots, ...
            NUM_LEVELS_M, NUM_PARTITIONS_J, nRegions, indexMatrix, isPredicting, nLevelsInSerial, NUM_WORKERS, verbose, varEps);  % Unsuppress output to print to command window
        elapsedTime = toc; % Unsuppress output to print to command window
        if verbose % Display the sumLogLikelihood
            sumLogLikelihood
        end
        % Save likelihood results
        save([resultsFilePath, 'Likelihood_Results'], 'sumLogLikelihood');
    otherwise
        error('Undefined calculationType. Code is not executed.')
end
if verbose
    disp('MRA execution completed');
end
end
