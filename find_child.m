function [ indexChildren ] = find_child( index, nRegions, NUM_PARTITIONS_J )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
        [ level, tile ] = find_level_tile( index, nRegions );
        levelChildren = level + 1;
        tileChildren=(NUM_PARTITIONS_J * (tile-1) + 1) : (NUM_PARTITIONS_J * tile);
        indexChildren=zeros(length(tileChildren),1, 'int64');
        for k = 1 : length(tileChildren)
            indexChildren(k) = find_index( levelChildren, tileChildren(k), nRegions );
        end

end

