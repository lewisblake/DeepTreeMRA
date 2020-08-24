function [] = distribute_data_structures(knots, paritions, nRegions, ...
    outputData, predictionLocations, NUM_PARTITIONS_J)

nRegionsAtFinestLevel  = nRegions(NUM_LEVELS_M);

cumulativeRegions = cumsum(nRegions);
totalRegions = sum(nRegions);
knotsContainer = cell(NUM_LEVELS_M, totalRegions);
outputDataContainer = cell();

% Dummy Index to get into original cell structures
dummyIndex = (cumulativeRegions(NUM_LEVELS_M) - nRegions(NUM_LEVELS_M) + 1);
for iRegion = 1:nRegionsAtFinestLevel 
    indexAncestry = find_ancestry(dummyIndex, nRegions, NUM_PARTITIONS_J);
    knotsContainer(1:NUM_LEVELS_M, iRegion) = knots([indexAncestry; dummyIndex],1);
    outputDataContainer = outputData();
    
    
    
    dummyIndex = dummyIndex + 1;
end


end
