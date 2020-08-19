function [ M_integer, nRegionsFinest, M, nRegions, totalRegions ] = find_num_levels_suggested( NUM_DATA_POINTS_n, NUM_KNOTS_r, NUM_PARTITIONS_J )
%% FIND_M_REQUIRED is a stand-alone function that estimates number of levels
%   Finds the number of levels for a given number of observations, number
%   of knots per region (NUM_KNOTS_r) and number of partitions per region (NUM_PARTITIONS_J).
%   The idea is to work backwards from the finest level by making the average
%   number of observations per region similar to the number of knots.
%
%	Input: NUM_DATA_POINTS_n, NUM_KNOTS_r, NUM_PARTITIONS_J
%
%	Output: M_integer, nRegionsFinest, M, nRegions, totalRegions

%%
% Find number of regions required at finest level
nRegionsFinest = NUM_DATA_POINTS_n/NUM_KNOTS_r;

% Solve for M which is m at the finest level using the J^M is the number of
% regions at the finest level

M = log(nRegionsFinest)/log(NUM_PARTITIONS_J);
M_integer = ceil(M);

% Comment: Matlab doesn't have a direct command to solve a logarithm for
% any base, so we use log(x)/log(b) in place of log_base_b(x)

mLevels = 0:M_integer; % Vector of levels
nRegions = NUM_PARTITIONS_J.^mLevels; % Regions (partitions) by level
totalRegions = sum(nRegions); % Total number of regions over all levels
end