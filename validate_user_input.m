function [] = validate_user_input(calculationType, NUM_LEVELS_M, NUM_PARTITIONS_J, NUM_KNOTS_r, offsetPercentage, NUM_WORKERS, NUM_LEVEL_ASSIGN_REGIONS_P, nXGrid, nYGrid, displayPlots, savePlots, verbose, resultsFilePath, plotsFilePath, fitRegressionModel, domainGeometry)
%% USER_INPUT_ERROR_HANDLING checks validity of user input
%
%   Input: NUM_LEVELS_M (double), NUM_PARTITIONS_J (double), NUM_KNOTS_r
%   (knots), offsetPercentage (double), nXGrid (double), nYGrid (double)
%
%   Output: [] (empty vector)
%
%
% Error handling on all input arguments except for dataSource which is
% checked within load_data.m
if ~strcmp(calculationType, {'prediction', 'likelihood', 'optimize', 'build_structure'})
    error('myfuns:validate_user_input: calulationType must either be prediction, likelihood, optimize, or build_structure')
elseif ~isa(NUM_LEVELS_M, 'double') || isinf(NUM_LEVELS_M) || (~isinf(NUM_LEVELS_M) && floor(NUM_LEVELS_M) ~= NUM_LEVELS_M)
    error('myfuns:validate_user_input: NUM_LEVELS_M must be a positive integer of class double.')
elseif ~isa(NUM_PARTITIONS_J, 'double') || (NUM_PARTITIONS_J ~= 2 && NUM_PARTITIONS_J ~= 4) || isinf(NUM_PARTITIONS_J) || (~isinf(NUM_PARTITIONS_J) && floor(NUM_PARTITIONS_J) ~= NUM_PARTITIONS_J)
    error('myfuns:validate_user_input: NUM_PARTITIONS_J must be either 2 or 4.')
elseif ~isa(NUM_KNOTS_r, 'double') || (~isinf(NUM_KNOTS_r) && floor(NUM_KNOTS_r) ~= NUM_KNOTS_r) || isinf(NUM_KNOTS_r)
    error('myfuns:validate_user_input: NUM_KNOTS_r must be a positive integer of class double.')
elseif ~isa(offsetPercentage, 'double') || (offsetPercentage < 0 || offsetPercentage >= 1) || isinf(offsetPercentage)
    error('myfuns:validate_user_input: offsetPercentage must be between 0 and 1.')
elseif ~isa(NUM_WORKERS, 'double') ||(~isinf(NUM_WORKERS) && floor(NUM_WORKERS) ~= NUM_WORKERS) || isinf(NUM_WORKERS) || NUM_WORKERS < 1 || abs(round(log(NUM_WORKERS)/log(NUM_PARTITIONS_J)) - log(NUM_WORKERS)/log(NUM_PARTITIONS_J)) % Last condition checks that NUM_WORKERS is power of J
    error('myfuns:validate_user_input: NUM_WORKERS must be a positive integer power of J of class double.')
elseif ~isa(NUM_LEVEL_ASSIGN_REGIONS_P, 'double') || (~isinf(NUM_LEVEL_ASSIGN_REGIONS_P) && floor(NUM_LEVEL_ASSIGN_REGIONS_P) ~= NUM_LEVEL_ASSIGN_REGIONS_P) || isinf(NUM_LEVEL_ASSIGN_REGIONS_P) || NUM_LEVEL_ASSIGN_REGIONS_P < 1 || NUM_LEVEL_ASSIGN_REGIONS_P > NUM_LEVELS_M
    error('myfuns:validate_user_input: NUM_LEVELS_ASSIGN_REGIONS_P must be a positive integer of class double less than or equal to NUM_LEVELS_M')
elseif ~isa(nXGrid, 'double') || (~isinf(nXGrid) && floor(nXGrid) ~= nXGrid) || isinf(nXGrid) || nXGrid < 1
    error('myfuns:validate_user_input: nXGrid must be a positive integer of class double.')
elseif ~isa(nYGrid, 'double') || (~isinf(nYGrid) && floor(nYGrid) ~= nYGrid) || isinf(nYGrid) || nYGrid < 1
    error('myfuns:validate_user_input: nYGrid must be a positive integer of class double.')
elseif ~islogical(displayPlots)
    error('myfuns:validate_user_input: displayPlots must be a boolean.')
elseif ~islogical(savePlots)
    error('myfuns:validate_user_input: savePlots must be a boolean.')
elseif ~islogical(verbose)
    error('myfuns:validate_user_input: verbose must be a boolean.')
elseif ~ischar(resultsFilePath)
    error('myfuns:validate_user_input: resultsFilePath must be a char. See default.')
elseif ~ischar(plotsFilePath)
    error('myfuns:validate_user_input: plotsFilePath must be a char. See default.')
elseif ~islogical(fitRegressionModel)
    error('myfuns:validate_user_input: fitRegressionModel must be a boolean.')
elseif ~strcmp(domainGeometry, {'plane', 'sphere'})
    error('myfuns:valuidate_user_input: domainGeometry must be either plane or sphere')
end

end
