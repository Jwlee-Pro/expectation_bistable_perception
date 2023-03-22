function nifti_mapping(val_in, crits, ROI_Gray, GrayCoord, save_loc, save_name)
[nx,ny,nz,nt] = size(ROI_Gray.img);

tempI = nan(nx,ny,nz); 
for iVox = 1:length(val_in)
    if crits(iVox)==1
        tempI(GrayCoord(iVox,1), GrayCoord(iVox,2), GrayCoord(iVox,3)) = val_in(iVox); 
    end
end
ROI_Gray.img = tempI; 
save_nii(ROI_Gray,[save_loc '/' save_name '.nii']); 




