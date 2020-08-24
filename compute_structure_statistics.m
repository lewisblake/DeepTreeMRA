function [ ] = compute_structure_statistics(outputData, NUM_WORKERS, resultsFilePath, plotsFilePath)
%% COMPUTE_STRUCTURE_STATISTICS
% This function calculates and reports summary statistics regarding the
% distribution of observations in regions at the finest level. By default, summary
% statistics are piped into a text file in the Results folder.
%
%   Input: knots, partitions, nRegions, outputData, indexMatrix
%
%   Output: [] (empty vector)
%
%% 

% Collect summary statistics on each worker in parallel
spmd(NUM_WORKERS)
    thisWorkersData = getLocalPart(outputData(:, labindex));
    [numObs, ~] = cellfun(@size, thisWorkersData);
    maxNumObsThisWorker = max(numObs);
    minNumObsThisWorker = min(numObs);
    meanNumObsThisWorker = mean(numObs);
    stdNumObsThisWorker = std(numObs);
    medianNumObsThisWorker = median(numObs);
    numZerosThisWorker = length(numObs) - nnz(numObs); % LB 3/17: this needs to be tested
end

% Gather summary statistics from all workers
gatheredMaxObs = gather(maxNumObsThisWorker);
gatheredMinObs = gather(minNumObsThisWorker);
gatheredMeanObs = gather(meanNumObsThisWorker);
gatheredStdObs = gather(stdNumObsThisWorker);
gatheredMedianObs = gather(medianNumObsThisWorker);
gatheredNumObs = gather(numObs);
gatheredNumZeros = gather(numZerosThisWorker);

% Compute summary statistics across workers
maxNumObsAllWorkers = max([gatheredMaxObs{:}]);
minNumObsAllWorkers = min([gatheredMinObs{:}]);
meanNumObsAllWorkers = mean([gatheredMeanObs{:}]);
stdNumObsAllWorkers = std([gatheredStdObs{:}]);
medianNumObsAllWorkers = median([gatheredMedianObs{:}]);
gatheredNumObsMat = [gatheredNumObs{:}];
gatheredNumZeros = sum(cellfun(@sum, gatheredNumZeros(:)));

%Create and save histogram object
hist = histogram(gatheredNumObsMat);
set(gca,'YScale','log')
saveas(hist, [plotsFilePath ,'numObsHist.fig'], 'fig'); % Currently set up for Unix-like systems.

%Create output text file
summaryStatsFileID = fopen(fullfile(resultsFilePath, 'structureSummaryStats.txt'), 'w');

% Write to file
fprintf(summaryStatsFileID, 'Summary statistics for observations assigned to each region across all workers:\n');
fprintf(summaryStatsFileID, 'Max: ');
fprintf(summaryStatsFileID,'%d \n',maxNumObsAllWorkers);
fprintf(summaryStatsFileID, 'Min: ');
fprintf(summaryStatsFileID, '%d \n', minNumObsAllWorkers);
fprintf(summaryStatsFileID, 'Mean: ');
fprintf(summaryStatsFileID, '%6.2f \n', meanNumObsAllWorkers);
fprintf(summaryStatsFileID, 'Standard Deviation: ');
fprintf(summaryStatsFileID, '%6.2f \n', stdNumObsAllWorkers);
fprintf(summaryStatsFileID, 'Median: ');
fprintf(summaryStatsFileID, '%6.2f \n', medianNumObsAllWorkers);
fprintf(summaryStatsFileID, 'Number of regions with zero observations: ');
fprintf(summaryStatsFileID, '%6.2f \n', gatheredNumZeros);

% Close output textfile
fclose(summaryStatsFileID);

end

