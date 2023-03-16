%% TS generator
% generate time-series for multiple ROIs 


% input
% TS_all = [# of all voxels * # of Time points] 
% pts    = [3 * # of voxels in ROI]

% output
% TS = [# of voxels in ROI * # of Time points]


function TS = TSGenerator_faster(TS_all, Coord_all, pts)
    
    
    TS = zeros(length(pts(1,:)),length(TS_all(1,:)));
    
    Coord_all_merge = 10000*Coord_all(1,:) + 100*Coord_all(2,:) + Coord_all(3,:) ;
    pts_merge = 10000*pts(1,:) + 100*pts(2,:) + pts(3,:);
    
    
    % Search voxels 
    for iVox = 1:length(pts(1,:))
        if ~isempty(find(Coord_all_merge == pts_merge(iVox)))
            IndVox = find(Coord_all_merge == pts_merge(iVox)); 
            TS(iVox, :) = TS_all(IndVox,:);
        end
	end
	TS = TS';
end

