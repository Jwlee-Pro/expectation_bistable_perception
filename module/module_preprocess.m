%% Preprocessing module (overview)
% 1. Standardize the BOLD signal 
% 2. Slow fMRI fluctuation correction

M_preproc.HPF.on = 1 ;

if M_preproc.HPF.on == 1
    M_preproc.HPF.freq_crit = nan; 
end


%% Standardize the BOLD signal 
% Run-wise mean centering 
% Normalize with std 


%% Slow fMRI fluctuation correction
% - This is not related to the neural response
% - Potential cause: MRI temporature, fatigue, ...

% Non-neural activity: outside of brain / csf
% Regress out, subtraction 



































