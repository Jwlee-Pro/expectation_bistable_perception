function plot_patch(x, data, data_std)


    if length(data(:,1)) == 1 
        % 1 X nT case 
         vertex = [[x ; ((data_std) + data)] [fliplr(x) ; fliplr((-(data_std) + data))]];
         face = 1:length(vertex); 

    else 
        % nT X 1 case 
         vertex = [[x' ; ((data_std') + data')] [fliplr(x') ; fliplr((-(data_std') + data'))]];
         face = 1:length(vertex); 
    end

%     patch('Faces',face,'Vertices',vertex','Facecolor',[0.75 0.75 0.75],'FaceAlpha',0.8,'line','none');
%     patch('Faces',face,'Vertices',vertex','Facecolor',[0.82 0.82 0.82],'line','none');
    patch('Faces',face,'Vertices',vertex','Facecolor',[0 0 0]+0.75);
    
end
