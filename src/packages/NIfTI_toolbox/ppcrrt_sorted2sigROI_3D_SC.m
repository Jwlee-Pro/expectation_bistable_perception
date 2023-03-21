clear all

%% preprocessing
clear all
proportion_of_correct	=	[];
ppcrrt					=	[];
load wGM_SVM_preslS1.mat
proportion_of_correct	=	[proportion_of_correct; proportion_of_correct_S];
ppcrrt					=	[ppcrrt; ppcrrt_S];	
load wGM_SVM_preslS2.mat
proportion_of_correct	=	[proportion_of_correct; proportion_of_correct_S];
ppcrrt					=	[ppcrrt; ppcrrt_S];	
load wGM_SVM_preslS3.mat
proportion_of_correct	=	[proportion_of_correct; proportion_of_correct_S];
ppcrrt					=	[ppcrrt; ppcrrt_S];	
save('wGM_SVM.mat','proportion_of_correct','ppcrrt')

%
clear all
proportion_of_correct	=	[];
ppcrrt					=	[];
load wGM_SVM_preslC1.mat
proportion_of_correct	=	[proportion_of_correct; proportion_of_correct_C];
ppcrrt					=	[ppcrrt; ppcrrt_C];	
load wGM_SVM_preslC2.mat
proportion_of_correct	=	[proportion_of_correct; proportion_of_correct_C];
ppcrrt					=	[ppcrrt; ppcrrt_C];	
load wGM_SVM_preslC3.mat
proportion_of_correct	=	[proportion_of_correct; proportion_of_correct_C];
ppcrrt					=	[ppcrrt; ppcrrt_C];	
save('wGM_SVM.mat','proportion_of_correct','ppcrrt')

%% predictability of preS and preC at significant ROIs in respect to current choicr

clear all
sigROI	=	cell(1,12);
nii = load_nii('wrinplane.nii');
cd ../SVMresults/SD4perc/preS
load wGM_SVM.mat
cd ../
load sigcoords.mat

ppcrrt_sorted(:,13:15)			=	1;
ppcrrt_sorted(find(ppcrrt_sorted<1 & ppcrrt_sorted>0))	=	1;
discardROI						=	find(mean(ppcrrt_sorted(:,1:12),2)==0);
ppcrrt_sorted(discardROI,:)		=	[];

ppcrrt_sorted_S					=	proportion_of_correct.*ppcrrt_sorted;
for TR		=	1:12
	sigROI{TR}	=	zeros(49,58,48);
	for ROI		=	1:size(ppcrrt_sorted_S,1)		
		sigROI{TR}(ppcrrt_sorted_S(ROI,13),ppcrrt_sorted_S(ROI,14),ppcrrt_sorted_S(ROI,15))	=	100*ppcrrt_sorted_S(ROI,TR);
	end
% 	sigROI{TR}	=	int16(sigROI{TR});
end 

cd ./preS
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

save_nii(sigROI_TR01,'sigROI_TR01_preS')
save_nii(sigROI_TR02,'sigROI_TR02_preS')
save_nii(sigROI_TR03,'sigROI_TR03_preS')
save_nii(sigROI_TR04,'sigROI_TR04_preS')
save_nii(sigROI_TR05,'sigROI_TR05_preS')
save_nii(sigROI_TR06,'sigROI_TR06_preS')
save_nii(sigROI_TR07,'sigROI_TR07_preS')
save_nii(sigROI_TR08,'sigROI_TR08_preS')
save_nii(sigROI_TR09,'sigROI_TR09_preS')
save_nii(sigROI_TR10,'sigROI_TR10_preS')
save_nii(sigROI_TR11,'sigROI_TR11_preS')
save_nii(sigROI_TR12,'sigROI_TR12_preS')

%%
clear all

sigROI	=	cell(1,12);
nii = load_nii('wrinplane.nii');
cd ../SVMresults/SD4perc/preC
load wGM_SVM.mat
cd ../
load sigcoords.mat

ppcrrt_sorted(:,13:15)			=	1;
ppcrrt_sorted(find(ppcrrt_sorted<1 & ppcrrt_sorted>0))	=	1;
discardROI						=	find(mean(ppcrrt_sorted(:,1:12),2)==0);
ppcrrt_sorted(discardROI,:)		=	[];

ppcrrt_sorted_C					=	proportion_of_correct.*ppcrrt_sorted;
for TR		=	1:12
	sigROI{TR}	=	zeros(49,58,48);
	for ROI		=	1:size(ppcrrt_sorted_C,1)		
		sigROI{TR}(ppcrrt_sorted_C(ROI,13),ppcrrt_sorted_C(ROI,14),ppcrrt_sorted_C(ROI,15))	=	100*ppcrrt_sorted_C(ROI,TR);
	end
% 	sigROI{TR}	=	int16(sigROI{TR});
end

cd ./preC
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

save_nii(sigROI_TR01,'sigROI_TR01_preC')
save_nii(sigROI_TR02,'sigROI_TR02_preC')
save_nii(sigROI_TR03,'sigROI_TR03_preC')
save_nii(sigROI_TR04,'sigROI_TR04_preC')
save_nii(sigROI_TR05,'sigROI_TR05_preC')
save_nii(sigROI_TR06,'sigROI_TR06_preC')
save_nii(sigROI_TR07,'sigROI_TR07_preC')
save_nii(sigROI_TR08,'sigROI_TR08_preC')
save_nii(sigROI_TR09,'sigROI_TR09_preC')
save_nii(sigROI_TR10,'sigROI_TR10_preC')
save_nii(sigROI_TR11,'sigROI_TR11_preC')
save_nii(sigROI_TR12,'sigROI_TR12_preC')


%% predictability of preS and preC at significant ROIs in respect to previous stimulus and choice

clear all
sigROI	=	cell(1,12);
nii = load_nii('wrinplane.nii');
cd ../SVMresults/SD4perc/preS
load wGM_SVM.mat

[sigcoords,ppcrrt_sorted] = ppcrrt2sigcoords(ppcrrt,proportion_of_correct);
ppcrrt_sorted(:,13:15)	=	1;
ppcrrt_sorted_S					=	proportion_of_correct.*ppcrrt_sorted;
for TR		=	1:12
	sigROI{TR}	=	zeros(49,58,48);
	for ROI		=	1:size(ppcrrt_sorted_S,1)		
		sigROI{TR}(ppcrrt_sorted_S(ROI,13),ppcrrt_sorted_S(ROI,14),ppcrrt_sorted_S(ROI,15))	=	100*ppcrrt_sorted_S(ROI,TR);
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

save_nii(sigROI_TR01,'sigROI_TR01_preS')
save_nii(sigROI_TR02,'sigROI_TR02_preS')
save_nii(sigROI_TR03,'sigROI_TR03_preS')
save_nii(sigROI_TR04,'sigROI_TR04_preS')
save_nii(sigROI_TR05,'sigROI_TR05_preS')
save_nii(sigROI_TR06,'sigROI_TR06_preS')
save_nii(sigROI_TR07,'sigROI_TR07_preS')
save_nii(sigROI_TR08,'sigROI_TR08_preS')
save_nii(sigROI_TR09,'sigROI_TR09_preS')
save_nii(sigROI_TR10,'sigROI_TR10_preS')
save_nii(sigROI_TR11,'sigROI_TR11_preS')
save_nii(sigROI_TR12,'sigROI_TR12_preS')

%
clear all

sigROI	=	cell(1,12);
cd ../../../Anatomy
nii = load_nii('wrinplane.nii');
cd ../SVMresults/SD4perc/preC
load wGM_SVM.mat

[sigcoords,ppcrrt_sorted] = ppcrrt2sigcoords(ppcrrt,proportion_of_correct);
ppcrrt_sorted(:,13:15)	=	1;
ppcrrt_sorted_C					=	proportion_of_correct.*ppcrrt_sorted;
for TR		=	1:12
	sigROI{TR}	=	zeros(49,58,48);
	for ROI		=	1:size(ppcrrt_sorted_C,1)		
		sigROI{TR}(ppcrrt_sorted_C(ROI,13),ppcrrt_sorted_C(ROI,14),ppcrrt_sorted_C(ROI,15))	=	100*ppcrrt_sorted_C(ROI,TR);
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

save_nii(sigROI_TR01,'sigROI_TR01_preC')
save_nii(sigROI_TR02,'sigROI_TR02_preC')
save_nii(sigROI_TR03,'sigROI_TR03_preC')
save_nii(sigROI_TR04,'sigROI_TR04_preC')
save_nii(sigROI_TR05,'sigROI_TR05_preC')
save_nii(sigROI_TR06,'sigROI_TR06_preC')
save_nii(sigROI_TR07,'sigROI_TR07_preC')
save_nii(sigROI_TR08,'sigROI_TR08_preC')
save_nii(sigROI_TR09,'sigROI_TR09_preC')
save_nii(sigROI_TR10,'sigROI_TR10_preC')
save_nii(sigROI_TR11,'sigROI_TR11_preC')
save_nii(sigROI_TR12,'sigROI_TR12_preC')