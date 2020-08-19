function [ index ] = find_index( level, tileNum, nRegions )
%% FIND_INDEX finds continuous index i
%   This function finds the continuous index given the level and
%   tile as inputs
%
%	Input: level, tileNum, nRegions
%
%	Output: index, continous index 
%%
index = sum(nRegions(1 : level-1)) + tileNum;
end