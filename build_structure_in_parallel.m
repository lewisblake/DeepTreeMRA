function [ knots, partitions, nRegions, outputData, predictionLocations, indexMatrix ] = build_structure_in_parallel( NUM_LEVELS_M, ...
    NUM_PARTITIONS_J, NUM_KNOTS_r, domainBoundaries, offsetPercentage, NUM_WORKERS, NUM_LEVEL_ASSIGN_REGIONS_P, verbose, varargin )
% BUILD_STRUCTURE_IN_PARALLEL
%   This function builds the hierarchical domain partitioning defined by
%   the MRA in parallel across NUM_WORKERS
%% Check inputs and display progress check

% Check number of optional input arguments does not exceed two
numVarArgs = length(varargin); % LB: outputs 2
if numVarArgs > 2
    error('myfuns:build_structure:TooManyInputs', ...
        'requires at most 2 optional inputs');
end

if verbose
    % Display progress check
    disp('Begining to build hierarchical grid structure in parallel...');
end
% Optional argument that can be passed is data
optArgs(1:numVarArgs) = varargin;
[ data, predictionVector ] = optArgs{:};
%% Set finest knot level
if numVarArgs == 2
    finestKnotLevel = NUM_LEVELS_M-1;
    indexEndFinestKnotLevel = NUM_PARTITIONS_J^(finestKnotLevel)-1;
else
    finestKnotLevel = NUM_LEVELS_M;
    indexEndFinestKnotLevel = NUM_PARTITIONS_J^(finestKnotLevel)-1;
end

%% Calculate quantities of interest
mLevels = 0:NUM_LEVELS_M-1; % Create a vector of levels
nRegions = NUM_PARTITIONS_J.^mLevels; % Vector of regions (partitions) at each level
cummulativeRegions = cumsum(nRegions);


%% Calculate number of knots in each direction
if floor(sqrt(NUM_KNOTS_r)) == sqrt(NUM_KNOTS_r)  % Assign knots
    nKnotsX0 = sqrt(NUM_KNOTS_r); nKnotsX = sqrt(NUM_KNOTS_r); % Number of knots in x-direction
    nKnotsY0 = sqrt(NUM_KNOTS_r); nKnotsY = sqrt(NUM_KNOTS_r); % Number of knots in y-direction
else
    nKnotsX0 = ceil(sqrt(NUM_KNOTS_r)); nKnotsX = ceil(sqrt(NUM_KNOTS_r));
    nKnotsY0 = NUM_KNOTS_r/nKnotsX0; nKnotsY = NUM_KNOTS_r/nKnotsX;
end


% Make assumption that tiles are equally distributed across workers
nRegionsAtFinestLevelForEachWorker = (nRegions(NUM_LEVELS_M)/NUM_WORKERS);
% Calculate quantities needed for nTotalRegionsAssignedToEachWorker
maxLevelOnASingleRow = sum(nRegions <= NUM_WORKERS); % How many times indices from a level are assigned to a worker
counter = 1:(NUM_LEVELS_M - maxLevelOnASingleRow); % Count number of times indicies from a level are going to be repeated in each indexMatrix column
nTotalRegionsAssignedToEachWorker = maxLevelOnASingleRow + sum(NUM_PARTITIONS_J.^counter);


%% Create matrix to store continuous index for all regions
[indexMatrix] = create_indexMatrix( NUM_LEVELS_M, NUM_PARTITIONS_J, nRegions, NUM_WORKERS, NUM_LEVEL_ASSIGN_REGIONS_P);
% Find the index within the indexMatrix corresponding to the finest level at which the knots are not set to the data.
indexOfFinestKnotLevelWithinIndexMatrix = find(indexMatrix(:,end) == indexEndFinestKnotLevel);
%% Pre-allocate memory for codistributed arrays
spmd(NUM_WORKERS)
    codistributionScheme = codistributor1d(2); % Distribute across the second dimension
    % Create codistributed cell arrays
    knots = cell(nTotalRegionsAssignedToEachWorker,  NUM_WORKERS, codistributionScheme);
    outputData = cell(nRegionsAtFinestLevelForEachWorker, NUM_WORKERS, codistributionScheme);
    partitions = cell(indexOfFinestKnotLevelWithinIndexMatrix+1,  NUM_WORKERS, codistributionScheme);
    predictionLocations = cell(nRegionsAtFinestLevelForEachWorker, NUM_WORKERS, codistributionScheme);
end
%% Construct zeroth level
xMin0 = domainBoundaries(1); xMax0 = domainBoundaries(2);
yMin0 = domainBoundaries(3); yMax0 = domainBoundaries(4);
% Edge buffer added to xMax0 and yMax0 to
% include all observations on the boundary at the zeroth level
xMax0 = xMax0 + (offsetPercentage/2)*(xMax0 - xMin0);
yMax0 = yMax0 + (offsetPercentage/2)*(yMax0 - yMin0);
% Create the knots at the coarsest resolution
[ knotsX, knotsY ] = create_knots(xMin0, xMax0, nKnotsX0, yMin0, yMax0, nKnotsY0, offsetPercentage);
% Each branch at the coarsest resolution will have the same knots
knots(1,1:NUM_WORKERS) = {[knotsX(:), knotsY(:)]}; % Knots at coarsest resolution are part of the parental hierarchy for all branches
% Create the paritition at the coarsest resolution
[ xMin, xMax, yMin, yMax ] = create_partition(xMin0, xMax0, yMin0, yMax0, NUM_PARTITIONS_J);
% Each branch at the coarest resolution will be contained within same
% region
partitions(1,1:NUM_WORKERS) = {[ xMin, xMax, yMin, yMax ]};

vectorOfRegionsAtFirstParallelLevel = cummulativeRegions(NUM_LEVEL_ASSIGN_REGIONS_P) - nRegions(NUM_LEVEL_ASSIGN_REGIONS_P) +1 : cummulativeRegions(NUM_LEVEL_ASSIGN_REGIONS_P);
matrixOfRegionsAtFirstParallelLevel = reshape(vectorOfRegionsAtFirstParallelLevel, [], NUM_WORKERS);
%% Loop creating partitions and placing knots
spmd(NUM_WORKERS)
    %% First loop up until level M-1 creating partitions and placing knots
    for iRow = 2 : indexOfFinestKnotLevelWithinIndexMatrix
        % Find this region's index
        indexCurrent = indexMatrix(iRow, labindex);
        % Find this region's parent and assign it to an int
        [~, ~, indexParent] = find_parent(indexCurrent, nRegions, NUM_PARTITIONS_J);
        % Find this region's parent's location in indexMatrix
        thisIndexParentInIndexMatrix = sum(indexMatrix(:, labindex) <= indexParent);
        % Get partition coordinates of parent
        parentPartitionsLocalPart = getLocalPart(partitions(thisIndexParentInIndexMatrix, labindex));
        xMin = parentPartitionsLocalPart{:}(:,1); xMax = parentPartitionsLocalPart{:}(:,2); yMin = parentPartitionsLocalPart{:}(:,3); yMax = parentPartitionsLocalPart{:}(:,4);
        foundChildren = find_children(indexParent, nRegions, NUM_PARTITIONS_J);
        thesePartitionBoundaries = find(foundChildren <= indexCurrent, 1, 'last'); 
        thisXMin = xMin(thesePartitionBoundaries); thisXMax = xMax(thesePartitionBoundaries);
        thisYMin = yMin(thesePartitionBoundaries); thisYMax = yMax(thesePartitionBoundaries);
        % Create knots
        [knotsX,knotsY] = create_knots(thisXMin, thisXMax, nKnotsX, thisYMin, thisYMax, nKnotsY, offsetPercentage);
        % Create partitions
        [ xMinTemp, xMaxTemp, yMinTemp, yMaxTemp ] = create_partition(thisXMin, thisXMax, thisYMin, thisYMax, NUM_PARTITIONS_J);
        % Assign knots
        knots(iRow, labindex) = {[knotsX(:),knotsY(:)]};
        % Assign parititions
        partitions(iRow, labindex) = {[ xMinTemp, xMaxTemp, yMinTemp, yMaxTemp ]};        
    end
    
    %% Prep for assigning data at finest resolution
    thisWorkerParallelAssignment = matrixOfRegionsAtFirstParallelLevel(:, labindex);
    beginParallelAssignmentInIndexMatrix = sum(indexMatrix(:, labindex) <= thisWorkerParallelAssignment(1));
	endParallelAssignmentInIndexMatrix = sum(indexMatrix(:, labindex) <= thisWorkerParallelAssignment(end));
    thisWorkersAssignmentPartitions = partitions(beginParallelAssignmentInIndexMatrix: endParallelAssignmentInIndexMatrix, labindex);
    % now we need the extreme lon/lat for this workers assignment
    % following works because dimensions of matrices agree and therefore
    % cell2mat can be applied. if we were considering more general
    % paritioning methods, this would not work and we may need to consider
    % another approach like I did in the C++ buildStructure.cpp for
    % Math-540. But the below implementation is vectorized so keep for now.
    thisWorkersAssignmentPartitionsAsMatrix = cell2mat(getLocalPart(thisWorkersAssignmentPartitions));
    workerXMin = min(thisWorkersAssignmentPartitionsAsMatrix(:,1));
    workerXMax = max(thisWorkersAssignmentPartitionsAsMatrix(:,2));
    workerYMin = min(thisWorkersAssignmentPartitionsAsMatrix(:,3));
    workerYMax = max(thisWorkersAssignmentPartitionsAsMatrix(:,4));
    % selecting only the subset needed by this worker
    thisWorkersDataRows = data(:,1) >= workerXMin & data(:,1) <= workerXMax & data(:,2) >= workerYMin & data(:,2) <= workerYMax;
    thisWorkersData = data(thisWorkersDataRows,:);
    data = []; % force data out of memory since now allocated subsets to workers
    
    %% Loop through final resolution (i.e., level m = M) assigning data to the knots   
    % Special construct to find knots for finest resolution region
    if numVarArgs == 2 % If data is sent to build_structure_in_parallel
        for iRow = indexOfFinestKnotLevelWithinIndexMatrix + 1 : nTotalRegionsAssignedToEachWorker
            % Find this region's index
            indexCurrent = indexMatrix(iRow, labindex);
            % Find this region's parent and assign it to an int
            [~, ~, indexParent] = find_parent(indexCurrent, nRegions, NUM_PARTITIONS_J);
            % Find this region's parent's location in indexMatrix
            thisIndexParentInIndexMatrix = sum(indexMatrix(:, labindex) <= indexParent);
            % Get partition coordinates of parent
            parentPartitionsLocalPart = getLocalPart(partitions(thisIndexParentInIndexMatrix, labindex));
            xMin = parentPartitionsLocalPart{:}(:,1); xMax = parentPartitionsLocalPart{:}(:,2); yMin = parentPartitionsLocalPart{:}(:,3); yMax = parentPartitionsLocalPart{:}(:,4);
            foundChildren = find_children(indexParent, nRegions, NUM_PARTITIONS_J);
            thesePartitionBoundaries = find(foundChildren <= indexCurrent, 1, 'last');
            thisXMin = xMin(thesePartitionBoundaries); thisXMax = xMax(thesePartitionBoundaries);
            thisYMin = yMin(thesePartitionBoundaries); thisYMax = yMax(thesePartitionBoundaries);
            % find the rows of this worker's data which are contained within
            % the boundaries of this region
            thisRegionsDataRows = thisWorkersData(:,1) >= thisXMin & thisWorkersData(:,1) <= thisXMax & thisWorkersData(:,2) >= thisYMin & thisWorkersData(:,2) <= thisYMax;
            % assign the lon/lat of the observations in this region to be the
            % knots
            knots(iRow, labindex) = {thisWorkersData(thisRegionsDataRows, 1:2)};
            % assign the observations in this region to the assocaited cell of
            % outputData
            outputData(iRow - indexOfFinestKnotLevelWithinIndexMatrix, labindex) = {thisWorkersData(thisRegionsDataRows, 3)};
            thisWorkersData(thisRegionsDataRows,:) = []; % Eliminate the data that has already been assigned to a region, speeds up subsequent searching
        
            % Vinay's addition
            if ~isnan(predictionVector) % If predicting
                predictionIndex = find(predictionVector(:,1) >= thisXMin & predictionVector(:,1) < thisXMax & predictionVector(:,2) >= thisYMin & predictionVector(:,2) < thisYMax); % Find the predictionVector locations within this region
                predictionLocations(iRow - indexOfFinestKnotLevelWithinIndexMatrix, labindex) = {predictionVector(predictionIndex,:)}; % Assign predictionVector locations within this region to corresponding entry of predictionLocations codistributed cell
            else
                predictionLocations = NaN;
            end          
        end
    end
    
end
if verbose
    % Progress indicator
    disp('Building the hierarchical structure complete.')
end
