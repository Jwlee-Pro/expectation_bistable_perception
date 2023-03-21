global session cmap session2 session_merge

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
session.ID{15} = 'EJC_1';
session.ID{16} = 'EJC_2';
session.ID{17} = 'JHR_1';
session.ID{18} = 'JHR_2';
session.ID{19} = 'XJWL001_1';
session.ID{20} = 'XKWC045_1';
session.ID{21} = 'XKWC045_2';



session2.ID{1} = 'XJWL001';
session2.ID{2} = 'XJWL002';
session2.ID{3} = 'XKWC045';
session2.ID{4} = 'XJWL005';
session2.ID{5} = 'XJWL006';
session2.ID{6} = 'XJWL007';
session2.ID{7} = 'XJWL008';
session2.ID{8} = 'XJWL009';
session2.ID{9} = 'XKWC068';

        

session_merge.ID{1}  = 'EY';                session_merge.index{1}  = {'EY110726'; 'EY110713'} ;    % fmri
session_merge.ID{2}  = 'GC';                session_merge.index{2}  = 'GC110616' ;                  % fmri
session_merge.ID{3}  = 'HO';                session_merge.index{3}  = 'HO110622' ;                  % fmri
session_merge.ID{4}  = 'HS';                session_merge.index{4}  = {'HS110615'; 'HS110713'} ;    % fmri
session_merge.ID{5}  = 'IR';                session_merge.index{5}  = {'IR110728'; 'IR110722'} ;    % fmri
session_merge.ID{6}  = 'SP';                session_merge.index{6}  = {'SP110729'; 'SP110818'} ;    % fmri
session_merge.ID{7}  = 'SL';                session_merge.index{7}  = {'SL110705'; 'SL110621'} ;    % fmri
session_merge.ID{8}  = 'RB';                session_merge.index{8}  = {'RB110627'; 'RB110617'} ;    % fmri
session_merge.ID{9}  = 'EJC';               session_merge.index{9}  = {'EJC_1'; 'EJC_2'} ;          % fmri
session_merge.ID{10} = 'JHR';               session_merge.index{10} = {'JHR_1'; 'JHR_2'} ;          % fmri
session_merge.ID{11} = 'XJWL001';           session_merge.index{11} = {'XJWL001_1'} ;               % fmri
session_merge.ID{12} = 'XKWC045';           session_merge.index{12} = {'XKWC045_1'; 'XKWC045_2'};   % fmri1, fmri2

session_merge.ID{13} = 'XJWL002';           session_merge.index{13} = 'XJWL002';                    % pupil
session_merge.ID{14} = 'XJWL005';           session_merge.index{14} = 'XJWL005';                    % pupil
session_merge.ID{15} = 'XJWL006';           session_merge.index{15} = 'XJWL006';                    % pupil
session_merge.ID{16} = 'XJWL007';           session_merge.index{16} = 'XJWL007';                    % pupil
session_merge.ID{17} = 'XJWL008';           session_merge.index{17} = 'XJWL008';                    % pupil
session_merge.ID{18} = 'XJWL009';           session_merge.index{18} = 'XJWL009';                    % pupil
session_merge.ID{19} = 'XKWC068';           session_merge.index{19} = 'XKWC068';                    % pupil
session_merge.ID{20} = 'XJWL001_pupil';     session_merge.index{20} = 'XJWL001';                 % pupil
session_merge.ID{21} = 'XKWC045_pupil';     session_merge.index{21} = 'XKWC045';                  % pupil


% Colormap
cmap_bone = copper(8); 
cmap.allmodel = [cmap_bone(end-5:end,:); [255 111 207]/255 ; [154 206 81]/255]; 
cmap.choice=[102 204 255;204 204 204; 255 204 102]./255;
cmap_c = [0.4 0.8 1; 0.8 0.8 0.8; 1 0.8 0.4]; 
cmap_clight = [0.7 0.9 1; 0.8 0.8 0.8; 1 0.9 0.7];