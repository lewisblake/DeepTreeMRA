function [indexMatrix] = create_indexMatrix(NUM_LEVELS_M, NUM_PARTITIONS_J, nRegions, NUM_WORKERS, NUM_LEVELS_SERIAL_S)
%% Create the indexMatrix
%  Input: NUM_LEVELS_M, NUM_PARTITIONS_J, nRegions, nWorkersUsed, nTotalRegionsAssignedToEachWorker
%  nLevelsInSerial, 
% 
%  Output: indexMatrix
%

nLevelToBeginParallel = NUM_LEVELS_SERIAL_S +1;
% Calculate necessary quantities
cummulativeRegions = cumsum(nRegions);
vecOfRegionsAtFirstParallelLevel = nRegions(NUM_LEVELS_SERIAL_S+1):cummulativeRegions(NUM_LEVELS_SERIAL_S+1);
matrixOfRegionsAtFirstParallelLevel = reshape(vecOfRegionsAtFirstParallelLevel, [], NUM_WORKERS);
% Calculate quantities needed for nTotalRegionsAssignedToEachWorker
maxLevelOnASingleRow = sum(nRegions <= NUM_WORKERS); % How many times indices from a level are assigned to a worker
counter = 1:(NUM_LEVELS_M - maxLevelOnASingleRow);
nTotalRegionsAssignedToEachWorker = maxLevelOnASingleRow + sum(NUM_PARTITIONS_J.^counter);

% Allocate space for indexMatrix
indexMatrix = nan(nTotalRegionsAssignedToEachWorker, NUM_WORKERS);

% Loop through all workers
for iWorker = 1:NUM_WORKERS   
    if NUM_WORKERS == length(vecOfRegionsAtFirstParallelLevel)%length(vectorOfRegionsAtNLevelsInSerial) 
        disp('num workers is equal to num regions at level which we distribute objects')
        indexMatrix(nLevelToBeginParallel, iWorker) = vecOfRegionsAtFirstParallelLevel(iWorker); % Regions above will not take vertical space in indexMatrix
        indexMatrix(1:NUM_LEVELS_SERIAL_S, iWorker) = find_ancestry(indexMatrix(nLevelToBeginParallel,  iWorker), nRegions, NUM_PARTITIONS_J);
        indexMatrix(nLevelToBeginParallel:end, iWorker) = find_branch_children(indexMatrix(nLevelToBeginParallel, iWorker), nRegions, NUM_PARTITIONS_J, NUM_LEVELS_M);
        
    elseif NUM_WORKERS > length(vecOfRegionsAtFirstParallelLevel) % Report error
	error('Error in create_indexMatrix: NUM_WORKERS can not be larger than nRegions(nLevelsInSerial');    
    elseif NUM_WORKERS < length(vecOfRegionsAtFirstParallelLevel)
        % Assign the regions allocated to this worker in a vector
        thisWorkerParallelAssignment = matrixOfRegionsAtFirstParallelLevel(:, iWorker);
        % Find all parents of the first parallel level. Store in vector
        tempVecOfParentsOfFirstParallelLevel=[];
        for i=1:length(thisWorkerParallelAssignment)           
            thisRegion  = thisWorkerParallelAssignment(i);
            [~, ~, thisRegionParents]= find_parent(thisRegion, nRegions, NUM_PARTITIONS_J);
            tempVecOfParentsOfFirstParallelLevel = [tempVecOfParentsOfFirstParallelLevel thisRegionParents];
        end
        tempVecOfParentsOfFirstParallelLevel = unique(tempVecOfParentsOfFirstParallelLevel);
        % Take the transpose of vector of all parents
        tempVecOfParentsOfFirstParallelLevel = tempVecOfParentsOfFirstParallelLevel';
        % Loop through the parent regions of the first level at which to compute in parallel.
        % Find all decendants of the parallel regions parents -  these are
        % unique to each column. Moreover, find the ancestors of these regions and store them both in vectors with repeats       
        tempVecOfParentsDescendents=[];
        tempVecOfParentsAncestors=[];
        for j = 1: length(tempVecOfParentsOfFirstParallelLevel)
            % Parent descendents
            theseParentsDescendents = find_branch_children(tempVecOfParentsOfFirstParallelLevel(j), nRegions, NUM_PARTITIONS_J, NUM_LEVELS_M);
            theseParentsDescendents = theseParentsDescendents(2:end); % Cut off first entry as it contains the index which was given as input
            tempVecOfParentsDescendents = [tempVecOfParentsDescendents theseParentsDescendents]; 
            % Parent ancestors
            theseParentsAncestors = find_ancestry(tempVecOfParentsOfFirstParallelLevel(j), nRegions, NUM_PARTITIONS_J);
            tempVecOfParentsAncestors = [tempVecOfParentsAncestors theseParentsAncestors];
        end
        
        % Select only the unique entries from both vectors
        tempVecOfParentsDescendents = unique(tempVecOfParentsDescendents);
        tempVecOfParentsAncestors = unique(tempVecOfParentsAncestors);
        
        % Stack them all in a vector and store it into the indexMatrix
        indexMatrix(:, iWorker) = [tempVecOfParentsAncestors; tempVecOfParentsOfFirstParallelLevel; tempVecOfParentsDescendents];
      
    end    
end

end

