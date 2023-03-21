sigROI	=	cell(1,12);
% nii = load_nii('wrinplane.nii');
cd ../SVMresults/SD4perc
load wGM_SVM.mat
load sigcoords.mat


for TR		=	1:12
	sigROI{TR}	=	zeros(49,58,48);
	for ROI		=	1:size(ppcrrt_sorted,1)		
		sigROI{TR}(ppcrrt_sorted(ROI,13),ppcrrt_sorted(ROI,14),ppcrrt_sorted(ROI,15))	=	100*ppcrrt_sorted(ROI,TR);
	end
% 	sigROI{TR}	=	int16(sigROI{TR});
end


sigROI_TR01		=	nii;	sigROI_TR02		=	nii;	sigROI_TR03		=	nii;	sigROI_TR04		=	nii;
sigROI_TR05		=	nii;	sigROI_TR06		=	nii;	sigROI_TR07		=	nii;	sigROI_TR08		=	nii;
sigROI_TR09		=	nii;	sigROI_TR10		=	nii;	sigROI_TR11		=	nii;	sigROI_TR12		=	nii;

sigROI_TR01.img		=	sigROI{1};
sigROI_TR02.img		=	sigROI{2};
sigROI_TR03.img		=	sigROI{3};
sigROI_TR04.img		=	sigROI{4};
sigROI_TR05.img		=	sigROI{5};
sigROI_TR06.img		=	sigROI{6};
sigROI_TR07.img		=	sigROI{7};
sigROI_TR08.img		=	sigROI{8};
sigROI_TR09.img		=	sigROI{9};
sigROI_TR10.img		=	sigROI{10};
sigROI_TR11.img		=	sigROI{11};
sigROI_TR12.img		=	sigROI{12};

save_nii(sigROI_TR01,'sigROI_TR01')
save_nii(sigROI_TR02,'sigROI_TR02')
save_nii(sigROI_TR03,'sigROI_TR03')
save_nii(sigROI_TR04,'sigROI_TR04')
save_nii(sigROI_TR05,'sigROI_TR05')
save_nii(sigROI_TR06,'sigROI_TR06')
save_nii(sigROI_TR07,'sigROI_TR07')
save_nii(sigROI_TR08,'sigROI_TR08')
save_nii(sigROI_TR09,'sigROI_TR09')
save_nii(sigROI_TR10,'sigROI_TR10')
save_nii(sigROI_TR11,'sigROI_TR11')
save_nii(sigROI_TR12,'sigROI_TR12')