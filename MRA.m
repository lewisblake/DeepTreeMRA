function [ sumLogLikelihood, predictions ] = MRA( theta, outputData, knots, ...
    NUM_LEVELS_M, NUM_PARTITIONS_J, nRegions, indexMatrix, isPredicting, NUM_LEVELS_SERIAL_S, NUM_WORKERS, verbose, varargin)
%% MRA.m is the main Multi-resolution approximation function
%
% Input: theta, data, knots, MAX_LEVEL_M, NUM_PARTITIONS_J, nRegions, varargin
%
% Output: sumLogLikelihood, predictions
%%
% Check number of optional input arguments
numVarArgsIn = length(varargin);
if numVarArgsIn > 2
    error('myfuns:create_prior:TooManyInputs', ...
        'requires at most 2 optional inputs');
end
% set defaults for optional inputs
optionalArguments = {0 NaN};
% overwrite the ones specified in varargin.
optionalArguments(1 : numVarArgsIn) = varargin;
[ varEps, predictionLocations ] = optionalArguments{:};
% Calculate key quantities
%totalRegions = sum(nRegions);
cumulativeRegions = cumsum(nRegions);
nLevelToBeginInParallel = NUM_LEVELS_SERIAL_S + 1;
logLikelihoodSum = 0;
nRegionsAtFinestLevelForEachWorker = (nRegions(NUM_LEVELS_M)/NUM_WORKERS);
% 1/7 LB: Calculate quantities needed for nTotalRegionsAssignedToEachWorker
maxLevelOnASingleRow = sum(nRegions <= NUM_WORKERS); % How many times indices from a level are assigned to a worker
counter = 1:(NUM_LEVELS_M - maxLevelOnASingleRow);
nTotalRegionsAssignedToEachWorker = maxLevelOnASingleRow + sum(NUM_PARTITIONS_J.^counter);
nTotalRegionsInSerial = cumulativeRegions(NUM_LEVELS_SERIAL_S);
lastIndexOfSecondFinestLevel = nRegions(NUM_LEVELS_M)-1;
lastRowBeforeFinestLevel = find(indexMatrix(:,end)==lastIndexOfSecondFinestLevel);
lastIndexOfSerialLevel = nRegions(NUM_LEVELS_SERIAL_S+1)-1;
[lastRowInSerial, ~] = find(indexMatrix(:,:) == lastIndexOfSerialLevel,1);


%% Pre-allocate space for codistributed arrays
if verbose
    disp('In MRA.m: Pre-allocating space for objects ...')
end
spmd(NUM_WORKERS)
    codistributionScheme = codistributor1d(2); % Distribute across the second dimension
    % Create codistributed arrays. LB: CHECK SIZES OF THESE AGAIN ONCE COMPLETE!
    RpriorChol = cell(nTotalRegionsAssignedToEachWorker,  NUM_WORKERS, codistributionScheme);
    KcB = cell(nTotalRegionsAssignedToEachWorker,  NUM_WORKERS, codistributionScheme);
    AtildePrevious = cell(nTotalRegionsAssignedToEachWorker,  NUM_WORKERS, codistributionScheme);
    wtildePrevious = cell(nTotalRegionsAssignedToEachWorker,  NUM_WORKERS, codistributionScheme);
    %logLikelihood = nan(nTotalRegionsAssignedToEachWorker, nWorkersUsed, codistributionScheme);
    if isPredicting % If predicting, pre-allocate space for necessary quantities
        posteriorPredictionMean =  cell(nRegionsAtFinestLevelForEachWorker, NUM_WORKERS, codistributionScheme); % Only need to store values at finest resolution
        posteriorPredictionVariance = cell(nRegionsAtFinestLevelForEachWorker, NUM_WORKERS, codistributionScheme); % Only need to store values at finest resolution
        Btilde = cell(nRegionsAtFinestLevelForEachWorker, NUM_WORKERS, codistributionScheme); % Only need to store values at finest resolution
        predictions = cell(nRegionsAtFinestLevelForEachWorker, NUM_WORKERS, codistributionScheme); % Only need to store values at finest resolution
        RposteriorChol = cell(nTotalRegionsAssignedToEachWorker-nRegionsAtFinestLevelForEachWorker,  NUM_WORKERS, codistributionScheme); % LB: Doesn't need to include entries for finest resolution
        KcholA = cell(nTotalRegionsAssignedToEachWorker-nRegionsAtFinestLevelForEachWorker,  NUM_WORKERS, codistributionScheme); % LB: Doesn't need to include finest resolution
        Kcholw = cell(nTotalRegionsAssignedToEachWorker-nRegionsAtFinestLevelForEachWorker,  NUM_WORKERS, codistributionScheme); % LB: Doesn't need to include entries for finest resolution
    else % Workaround to pass correct object to create_prior
        predictionLocations = num2cell(nan(nTotalRegionsAssignedToEachWorker, NUM_WORKERS, codistributionScheme));
    end  
end
%% Pre-allocate space for cell arrays used in serial computations
% For creating the prior, except for knots, which are made by
% process_knots_for_serial()
%RpriorCholSerial = cell(nTotalRegionsInSerial,1);
%KcBSerial = cell(nTotalRegionsInSerial, 1);
% For calculating the posterior...

%% Create the prior distribution
% Assuming here that the last region to compute in serial is located in the last indexMatrix column
%knotsSubset = gather(knots(1:lastRowInSerial,:));
%indexMatrixSubset = indexMatrix(1:lastRowInSerial,:);

%% Process knots for serial computation - NOT sure if needed 8/24/20
%[ knotsSerial ] = process_knots_for_serial( knotsSubset, indexMatrixSubset, NUM_LEVELS_SERIAL_S, nRegions);

%%
if verbose
   disp('Creating the prior in parallel ...')
end
%% Main spmd block to create prior and parallel portion of posterio
spmd(NUM_WORKERS)
    %% First part of the MRA algorithm: Create prior from "coarsest to finest"
    mCounterIndex = 1; % LB: potential room for optimization.
    for iRow = 1 :  nTotalRegionsAssignedToEachWorker
        indexCurrent = indexMatrix(iRow, labindex);
        % Find the ancestry for this jRegion
        indexAncestry = find_ancestry(indexCurrent, nRegions, NUM_PARTITIONS_J);
        % Fill vector with location of indexAncestry in indexMatrix
        indexAncestryInIndexMatrixBool = ismember(indexMatrix(1:iRow, labindex), indexAncestry); % LB: this should work because the indices of the sub matrix do not change as we go down from top to bottom. Doesn't work when the submatrix goes from bottom to top in creating the posterior
        indexAncestryInIndexMatrix = find(indexAncestryInIndexMatrixBool ~= 0); % keep only nonzeros which correspond to the indices of indexAncestry in indexMatrix    
        
        % Get local part of the objects needed to create the prior
        knotsLocalPart = getLocalPart(knots([indexAncestryInIndexMatrix;iRow], labindex));        
        RpriorCholLocalPart = getLocalPart(RpriorChol(indexAncestryInIndexMatrix, labindex));       
        KcBLocalPart = getLocalPart(KcB(indexAncestryInIndexMatrix, labindex));

        if indexCurrent < nRegions(NUM_LEVELS_M) % If not dealing with an index on the finest resolution...
        % don't send create_prior any data or predictionLocations
        [thisRpriorChol, thisKcholBchol, ~, ~, ~] = create_prior(theta, ...
            NUM_LEVELS_M, knotsLocalPart, RpriorCholLocalPart,...
            KcBLocalPart, [], varEps, []);
        else
            dataLocalPart = getLocalPart(outputData(mCounterIndex, labindex));
            predictionLocationsLocalPart = getLocalPart(predictionLocations(mCounterIndex, labindex));
            % Create the prior
            [thisRpriorChol, thisKcholBchol, thisAtj, thiswtj, thisRetLikPred] = create_prior(theta, ...
                NUM_LEVELS_M, knotsLocalPart, RpriorCholLocalPart,...
                KcBLocalPart, dataLocalPart{:}, varEps, predictionLocationsLocalPart{:});

            AtildePrevious(iRow, labindex) = {thisAtj};
            wtildePrevious(iRow, labindex) = {thiswtj};
            if isPredicting  % If predicting
                posteriorPredictionMean(mCounterIndex, labindex) = {thisRetLikPred{1}};
                posteriorPredictionVariance(mCounterIndex, labindex) = {thisRetLikPred{2}};
                Btilde(mCounterIndex, labindex) = {thisRetLikPred{3}};
            else % If not predicting
                logLikelihoodSum = logLikelihoodSum + thisRetLikPred;
            end     
            mCounterIndex = mCounterIndex + 1;
        end
        RpriorChol(iRow, labindex) = {thisRpriorChol};
        KcB(iRow, labindex) = {thisKcholBchol};
    end
    if verbose && labindex == 1
        disp('Prior calculation complete.');
        disp('Creating the prior in parallel ...');
    end

    %% Second part of MRA: Caclulate the poster from "finest to coarsest".
    for iRow = lastRowBeforeFinestLevel:-1:maxLevelOnASingleRow
        index = indexMatrix(iRow, labindex);
        [indexChildren] = find_children(index, nRegions, NUM_PARTITIONS_J);
        % Fill vector with location of indexAncestry in indexMatrix
        indexChildrenInIndexMatrixBool = ismember(indexMatrix(:, labindex), indexChildren); % children must be further down in indexMtrix than iRow
        indexChildrenInIndexMatrix = find(indexChildrenInIndexMatrixBool ~= 0); % keep only nonzeros which correspond to the indices of indexAncestry in indexMatrix    
        
        
        RpriorCholj = getLocalPart(RpriorChol(iRow, labindex));
        wtildeChildren = getLocalPart(wtildePrevious(indexChildrenInIndexMatrix, labindex));
        AtildeChildren = getLocalPart(AtildePrevious(indexChildrenInIndexMatrix, labindex));
        
        % Calculate posterior_inference()
        [ wtildeCurrentj, AtildeCurrentj, logLikelihoodj, ...
            RposteriorCholj, Kcholwj, KcholAj ] = posterior_inference(RpriorCholj{:}, ...
            wtildeChildren, AtildeChildren);
        
        % LB: changed from wtileCurrent to wtildePrevious. Seems to be
        % working however need to check back later
        wtildePrevious(iRow,labindex) = {wtildeCurrentj};
        AtildePrevious(iRow, labindex) = {AtildeCurrentj};
        
        if isPredicting
            RposteriorChol(iRow, labindex) = {RposteriorCholj};
            Kcholw(iRow, labindex) = {Kcholwj};
            KcholA(iRow, labindex) = {KcholAj};
        else
            logLikelihoodSum = logLikelihoodSum + logLikelihoodj;
        end
    end
    
end
% Progress indicator
if verbose
   disp('Parallel section of the posterior complete.');
   % Serial section of the posterior
   disp('Calculating the serial section of the posterior ...');
end
%% "Serial Section"
% Now the question becomes what to do with the "serial" portion of the
% posterior. Want to minimize gathering from workers to clients. Also want
% to do "serial" portion as little as possible, so that we perofm in
% parallel until the number of regions at the level is equal to NUM_WORKERS

% Below is code to only gather the subsets of the distributed matrices
% needed on the client. nd2sub() is to get the sub-scripting index of a sub-matrix, and then sub2ind() is for taking these subscripts and relating them to linear indices in the entire matrix.
% Gather Rprior first since only need this level's entries
[~, indicesForRpriorChol] = unique(indexMatrix(1:maxLevelOnASingleRow, :),'first');
[rowSerialPortionForRpriorChol, colSerialPortionForRpriorChol] = ind2sub(size(indexMatrix(1:maxLevelOnASingleRow, :)), indicesForRpriorChol );
gatheredRprioChol = gather(RpriorChol(sub2ind(size(RpriorChol), rowSerialPortionForRpriorChol, colSerialPortionForRpriorChol)));

% Gather wtildePrevious, AtildePrevious
maxRowOfDistributedArraysToRetrieve = maxLevelOnASingleRow + NUM_PARTITIONS_J;
indexMatrixSubset = indexMatrix(1:maxRowOfDistributedArraysToRetrieve, :);
[~, indicesForwtildeAndAtilde] = unique(indexMatrixSubset,'first');
[rowSerialPortionForwtildeAndAtilde, colSerialPortionForwtildeAndAtilde] = ind2sub(size(indexMatrixSubset), indicesForwtildeAndAtilde );
% wtildePrevious and AtildePrevious same size, so use one to get linear
% indices for both
linearIndicesForwtildeAndAtilde = sub2ind(size(wtildePrevious), rowSerialPortionForwtildeAndAtilde, colSerialPortionForwtildeAndAtilde);
gatheredwtildePrevious = gather(wtildePrevious(linearIndicesForwtildeAndAtilde));
gatheredAtildePrevious = gather(AtildePrevious(linearIndicesForwtildeAndAtilde));
% Gather logLikelikehoodSum - composite object used to track likelihood in
% spmd block
% Intialize sumLogLikelihood to collect distributed logLikelihoodSum
gatheredSum = gather(logLikelihoodSum);
sumLogLikelihood = sum(cellfun(@sum, gatheredSum(:)));
gatheredSum = []; % force freeing memory. small so may remove.
% Loop through remaining regions in serial
for iLevel = maxLevelOnASingleRow:-1:1
    for jRegion = (cumulativeRegions(iLevel) - nRegions(iLevel) + 1) : cumulativeRegions(iLevel)
        [ indexChildren ] = find_children( jRegion, nRegions, NUM_PARTITIONS_J );
        % Calculate posterior quantities
        RpriorCholj = gatheredRprioChol{jRegion};
        wtildeChildren = gatheredwtildePrevious(indexChildren);
        AtildeChildren = gatheredAtildePrevious(indexChildren);
        
        % Calculate posterior_inference()
        [ wtildeCurrentj, AtildeCurrentj, logLikelihoodj, ...
            RposteriorCholj, Kcholwj, KcholAj ] = posterior_inference( RpriorCholj, ...
            wtildeChildren, AtildeChildren );
        
        gatheredwtildePrevious{jRegion} = wtildeCurrentj;
        gatheredAtildePrevious{jRegion} = AtildeCurrentj;
        
        [firstRowContainingThisRegion,firstColContainingThisRegion] = find(indexMatrixSubset(:,:)==jRegion, 1 );  
        nTimesEachIndexIsRepeatedThisRow = ceil(NUM_WORKERS/nRegions(iLevel)); % 1/7 LB: added ceil() function around calculation
    
        if isPredicting % If predicting
            RposteriorChol(firstRowContainingThisRegion, firstColContainingThisRegion:firstColContainingThisRegion + nTimesEachIndexIsRepeatedThisRow - 1) = {RposteriorCholj};
            Kcholw(firstRowContainingThisRegion, firstColContainingThisRegion:firstColContainingThisRegion + nTimesEachIndexIsRepeatedThisRow - 1) = {Kcholwj};
            KcholA(firstRowContainingThisRegion, firstColContainingThisRegion:firstColContainingThisRegion + nTimesEachIndexIsRepeatedThisRow - 1) = {KcholAj};
        else
            sumLogLikelihood = sumLogLikelihood + logLikelihoodj; % Update the sumLoglikelihood 
        end 
        
    end
end
% Force freeing up of memory
indexMatrixSubset = [];
if verbose
   disp('Serial section of the posterior complete.');
end

%% Spatial prediction
if isPredicting
    if verbose
       disp('Beginning the spatial prediction ...')
    end
    spmd(NUM_WORKERS)
        for iRow = lastRowBeforeFinestLevel+1:nTotalRegionsAssignedToEachWorker
            mCounterIndex = iRow - lastRowBeforeFinestLevel;
            if NUM_LEVELS_M > 0
                % Set up appropriate indicies
                index = indexMatrix(iRow, labindex);
                indexAncestry = find_ancestry(index, nRegions,NUM_PARTITIONS_J);
                
                % Fill vector with location of indexAncestry in indexMatrix
                indexAncestryInIndexMatrixBool = ismember(indexMatrix(1:iRow, labindex), indexAncestry); % LB: this should work because the indices of the sub matrix do not change as we go down from top to bottom. Doesn't work when the submatrix goes from bottom to top in creating the posterior
                indexAncestryInIndexMatrix = find(indexAncestryInIndexMatrixBool ~= 0); % keep only nonzeros which correspond to the indices of indexAncestry in indexMatrix    
        
                % Collect the appropriate inputs to make predictions at
                % this region
                thisPosteriorPredictionMean = getLocalPart(posteriorPredictionMean(mCounterIndex, labindex));
                thisPosteriorPredictionVariance = getLocalPart(posteriorPredictionVariance(mCounterIndex, labindex));
                thisBtilde = getLocalPart(Btilde(mCounterIndex, labindex));
                thisRposteriorChol = getLocalPart(RposteriorChol(indexAncestryInIndexMatrix, labindex));
                thisKcholA = getLocalPart(KcholA(indexAncestryInIndexMatrix, labindex));
                thisKcholw = getLocalPart(Kcholw(indexAncestryInIndexMatrix, labindex));
                % Use spatial_prediction() to make predictions. Place
                % output into a cell by wrapping function in {}.
                predictions(mCounterIndex,labindex) = {spatial_prediction(thisPosteriorPredictionMean{:}, ...
                    thisPosteriorPredictionVariance{:}, thisBtilde{:}, thisRposteriorChol, ...
                    thisKcholA, thisKcholw)};
            else
                predictions(mCounterIndex, labindex) = {[posteriorPredictionMean(mCounterIndex,labindex), posteriorVariance(mCounterIndex,labindex)]};
            end
        end        
    end
    if verbose
       disp('Spatial prediction complete.');
    end
end
   if verbose
      disp('MRA.m execution complete.');
   end
end
