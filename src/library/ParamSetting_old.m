global session cmap

% Subject information 
session.ID{1} = 'EY110726';
session.ID{2} = 'GC110616';
session.ID{3} = 'HO110622';
session.ID{4} = 'HS110615';
session.ID{5} = 'IR110728';
session.ID{6} = 'SP110729';
session.ID{7} = 'SL110705'; 
session.ID{8} = 'RB110627'; 
session.ID{9} = 'EY110713'; 
session.ID{10} = 'HS110713'; 
session.ID{11} = 'IR110722'; 
session.ID{12} = 'RB110617'; 
session.ID{13} = 'SL110621'; 
session.ID{14} = 'SP110818'; 
% session.ID{15} = 'EJC_1';
% session.ID{16} = 'EJC_2';
% session.ID{17} = 'JHR_1';
% session.ID{18} = 'JHR_2';
% session.ID{19} = 'XJWL001_1';
% session.ID{20} = 'XKWC045_1';
% session.ID{21} = 'XKWC045_2';


% Colormap
cmap_bone = copper(8); 
cmap.allmodel = [cmap_bone(end-5:end,:); [255 111 207]/255 ; [154 206 81]/255]; 
cmap.choice=[102 204 255;204 204 204; 255 204 102]./255;
cmap_c = [0.4 0.8 1; 0.8 0.8 0.8; 1 0.8 0.4]; 
cmap_clight = [0.7 0.9 1; 0.8 0.8 0.8; 1 0.9 0.7];