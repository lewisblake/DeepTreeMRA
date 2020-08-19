function [indexsOfAllChildrenOnThisBranch] = find_branch_children(index, nRegions, ...
    NUM_PARTITIONS_J, NUM_LEVELS_M)
%% FIND_BRANCH_CHILDREN
% This function finds the indicies of all the children for a specific
% index. This function goes in the opposite direction from find_ancestry.
% That is, given an index, find all further indexes that have the given
% index within their ancestry.
%
%   INPUT: index, nRegions, NUM_PARTITIONS_J, NUM_LEVELS_M 
%
%   OUTPUT: indexsOfAllChildrenOnThisBranch (vector)
%
%%
indexsOfAllChildrenOnThisBranch = [index]; % All regions have this index as common ancestor, so we start here.
[thisLevel, ~] = find_level_tile(index, nRegions); 
iLevel = thisLevel:NUM_LEVELS_M;

for iThisLevel = 1:length(iLevel) % Loop iThisLevel through counter index of levels finer than thisLevel
    indexOfAllChildrenForThisLevel = [];
    % Store the indices for appropriate level in a vector
    indicesOfRegionsAtThisLevel = indexsOfAllChildrenOnThisBranch(NUM_PARTITIONS_J^(iThisLevel-1): NUM_PARTITIONS_J^(iThisLevel)-1);
    % Then, use jThisRegion to index into this vector and get the actual index
    % needed to calculate find_children. Append this to SOMETHING
    for jThisRegion = 1:length(indicesOfRegionsAtThisLevel) % Loop jThisRegion through all region indexs?
        thisRegion = indicesOfRegionsAtThisLevel(jThisRegion);        
        thisRegionsChildren = find_children(thisRegion, nRegions, NUM_PARTITIONS_J);   
        indexOfAllChildrenForThisLevel = [indexOfAllChildrenForThisLevel; thisRegionsChildren];
    end
    % Append the indicies for children of this level to the vector
    % containing the indicies of all children on this branch
    indexsOfAllChildrenOnThisBranch = [indexsOfAllChildrenOnThisBranch; indexOfAllChildrenForThisLevel];
end

indexsOfAllChildrenOnThisBranch = sort(indexsOfAllChildrenOnThisBranch); % Sort final vector.
