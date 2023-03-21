%% TS generator
% generate time-series for multiple ROIs 


% input
% TS_all = [# of all voxels * # of Time points] 
% pts    = [3 * # of voxels in ROI]

% output
% TS = [# of voxels in ROI * # of Time points]


function TS = TSGenerator(TS_all, Coord_all, pts)

    TS = zeros(length(pts(1,:)),length(TS_all(1,:)));
    
    % Search voxels 
    TSindex = 0; 
    for iVox = 1:length(pts(1,:))
        for iVox_all = 1:length(Coord_all(1,:))
            tempIndex = Coord_all(:,iVox_all)-pts(:,iVox);
            if (tempIndex(1,1)==0) & (tempIndex(2,1)==0) & (tempIndex(3,1)==0)
                TSindex = TSindex + 1; 
                TS(TSindex, :) = TS_all(iVox_all,:);
            end
        end
	end
	TS = TS';
end

