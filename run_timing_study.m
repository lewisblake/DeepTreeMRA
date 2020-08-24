% Start cluster
cluster = parcluster('NCAR Cheyenne - R2018a')

% Create Communcaiting Jobs
j1 = createCommunicatingJob(cluster);
j2 = createCommunicatingJob(cluster);
j3 = createCommunicatingJob(cluster);
j4 = createCommunicatingJob(cluster);
j5 = createCommunicatingJob(cluster);

% Create tasks
createTask(j1, @main, 1, {});
createTask(j2, @main, 1, {});
createTask(j3, @main, 1, {});
createTask(j4, @main, 1, {});
createTask(j5, @main, 1, {});

% Submit jobs
submit(j1); submit(j2); submit(j3); submit(j4); submit(j5);

% Wait for the jobs to finish
wait(j1); wait(j2); wait(j3); wait(j4); wait(j5);

% Fetch outputs
out1 = fetchOutputs(j1);
out2 = fetchOutputs(j2);
out3 = fetchOutputs(j3);
out4 = fetchOutputs(j4);
out5 = fetchOutputs(j5);


% Calculate mean and standard deviatation
outVector = [out1{:} out2{:} out3{:} out4{:} out5{:}];
outBar = mean(outVector);
outStd = std(outVector);

% Display the mean and standard deviation
disp(['The sample mean is: ', num2str(outBar)]);
disp(['The standard deviation is: ', num2str(outStd)]);

% Delete the jobs
delete(j1); delete(j2); delete(j3); delete(j4); delete(j5);