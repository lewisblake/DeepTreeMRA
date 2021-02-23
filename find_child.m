function [ indexChildren ] = find_child( index, nRegions, NUM_PARTITIONS_J )
% FIND_CHILDREN
%   this function finds the children regions for a given region.
%   That is, for a given region, it finds the indices of paritions of the
%   given region at finer resolutions.
        [ level, tile ] = find_level_tile( index, nRegions );
        levelChildren = level + 1;
        tileChildren=(NUM_PARTITIONS_J * (tile-1) + 1) : (NUM_PARTITIONS_J * tile);
        indexChildren=zeros(length(tileChildren),1, 'int64');
        for k = 1 : length(tileChildren)
            indexChildren(k) = find_index( levelChildren, tileChildren(k), nRegions );
        end

end

