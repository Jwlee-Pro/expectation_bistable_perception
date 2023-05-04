%% Locus of expectation in medial OFC

% Author: Joonwon Lee 
% Date: 2023.5.3 (last revised)

clear all; close all; clc; 


%% [0] Module details
runModules = {};
runModules.compute_corr = 1; 
runModules.criterion_filter = 1; 


module_preprocess
module_voxelselection
module_svr_decoding
module_timecourse_extraction
module_behavior_tests
module_glm_summary



