function output = stmap_process(DataIn)

output = zeros(size(DataIn{1}{1})); 

for ix = 1:length(output(:,1))
    for iy = 1:length(output(1,:))
        tempX = [] ; 
        for iSub = 1:length(DataIn)
            if ~isempty(DataIn{iSub})
                for iTr = 1:length(DataIn{iSub})
                    tempX = [tempX DataIn{iSub}{iTr}(ix,iy)];
                end
            end
        end
        output(ix,iy) = nanmean(tempX); 
    end
end

        












% for iSub = 1:length(DataIn) 
%     if ~isempty(DataIn{iSub})
%         k_temp = zeros(size(DataIn{iSub}{1})); 
%         for iTr = 1:length(DataIn{iSub})
%             k_temp = k_temp + DataIn{iSub}{iTr}; 
%         end
%         k_temp = k_temp/length(DataIn{iSub}); 
%     end
%     
%     output_ind{iSub} = k_temp; 
%     output = output + output_ind{iSub}; 
% end
% 
% output = output/length(DataIn) ;

    
        