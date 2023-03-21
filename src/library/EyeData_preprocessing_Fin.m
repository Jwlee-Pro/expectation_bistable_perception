function corrected = EyeData_preprocessing_Fin(Eyelink_data, param, Sbj_info)
    

% Sampling rate
    sampleRate = Eyelink_data.RECORDINGS(1,1).sample_rate;
    
    
% Lock start time to 0
    baseNum = ceil(param.tStimOnOff{1,1}(1,1)*500); 
    
    
    
% average Stimulus On and Off duration 
    avg_OnTime = mean(param.tStimOnOff{1,1}(:,2)-param.tStimOnOff{1,1}(:,1));   % mean stimulus on time    
    temp_OffDuration = [];
    for iT = 1:length(param.tStimOnOff{1,1}(:,1))-1
        valOff = param.tStimOnOff{1,1}(iT+1,1)-param.tStimOnOff{1,1}(iT,2);
        temp_OffDuration = [temp_OffDuration;valOff];
    end
    avg_OffTime = mean(temp_OffDuration);

    total_time = (param.tStimOnOff{1,1}(end,2)+avg_OffTime);
    total_samples_eyetracker = total_time * sampleRate; 
%     fprintf('Maximum error in eyetracker data is %gms\n',abs(length(Eyelink_data.FSAMPLE.time)-total_samples_eyetracker)*1000/500); % ms

% Raw time data 
    corrected.time = ((Eyelink_data.FSAMPLE.time - Eyelink_data.FSAMPLE.time(1))*1)'; % time axis in ms
    
    
% Gaze position data and pupil area data 
    fixData.gx = Eyelink_data.FSAMPLE.gx(1,:)';
    fixData.gy = Eyelink_data.FSAMPLE.gy(1,:)';
    fixData.pa = Eyelink_data.FSAMPLE.pa(1,:)';

    gazeEye = horzcat(fixData.gx, fixData.gy);


% Apply butterworth filter spec. lowpass, 10hz
    Wn = 500/2; 
    [btb, bta] = butter(3, 10/Wn, 'low');

    
    corrected.Pa = single(filtfilt(btb, bta, double(fixData.pa(:,1)))); % smoothing the pupil size

    
% Regression: gaze position with pupil area + blink detection 
    rawPa = fixData.pa(:,[1 1]); % 
    tmpPa = [corrected.Pa corrected.Pa];
    tmpRes = [];
    for ii = 1:2
        % blink & error detection
        % blink criterions: (1)Absence of pupil area, (2)rapid change in pupil
        % area, (3)Wrong value from gaze data(too high). 
        Crit_pa = 200; % Pupil area should be larger than Crit_pa when fully exposed (at least!)
        Crit_diffpa = 700;
        Crit_g = 100000; 
        blinkPt = find((rawPa(:,ii) < Crit_pa) | vertcat((abs(diff(rawPa(:,ii))) > Crit_diffpa),0)|(abs(gazeEye(:,ii))) > Crit_g);

        % mark samples +- 200ms around the blinks
        flagConta = false(length(corrected.Pa),1);
        for ij = 1:length(blinkPt)
            if blinkPt(ij) <= 100
                flagConta(1:blinkPt(ij)) = true;
            elseif blinkPt(ij) >= (length(tmpPa)-100)
                flagConta(blinkPt(ij):end) = true;
            else
                flagConta(blinkPt(ij)+(-100:100)) = true;
            end
        end
        
        % merging nearby blinks 
        BlinkStart = find(diff(flagConta) == 1); 
        for iBl_temp = 2:length(find(diff(flagConta) == 1))
            difftemp = diff(flagConta); 
            BlinkEnd_last = find(difftemp(1:BlinkStart(iBl_temp)) == -1,1,'last'); 
            if BlinkStart(iBl_temp)-BlinkEnd_last < 100
                flagConta(BlinkEnd_last-10:BlinkStart(iBl_temp)+10) = 1;
            end
        end

        regPa = (tmpPa(:,ii) - nanmean(tmpPa(~flagConta,ii))) / nanstd(tmpPa(~flagConta,ii));
        matPa = [ones(length(tmpPa),1), regPa, regPa.^2];
        [b, bint, resi, rint, rstat] = regress(gazeEye(~flagConta,ii), matPa(~flagConta,:));
        b_values{ii} = b;
        % resi = tmpEye(:,ii) - (matPa(:,2)*b(2)) - (matPa(:,3)*b(3));
        resi = gazeEye(:,ii) - matPa*b; 
        resi(flagConta) = nan;
        gazeEye(flagConta,ii) = nan; % original (just for comparison)
        
        tmpRes = horzcat(tmpRes, resi);
    end
    corrected.Gx = tmpRes(:,1);
	corrected.Gy = tmpRes(:,2);
    corrected.BlinkIndex = flagConta; 
    
    for iT = 1:length(corrected.Gx)
        if isnan(corrected.Gx(iT))|isnan(corrected.Gy(iT))
            corrected.Pa(iT) = nan;
        end
    end


% Fill in the blank period (linear interpolation) 
    corrected.Pa_padded = Blankpadding(corrected.Pa);
    corrected.Gx_padded = Blankpadding(corrected.Gx);
	corrected.Gy_padded = Blankpadding(corrected.Gy);

    
% Calculate Blink rate 
    RateWindow = 500*5; 
%     % version 1: average at mid point ( -Window/2 ~ +Window/2 )
%     corrected.BlinkRate = zeros(length(corrected.time),1); 
%     for iF = 1:length(corrected.Pa)-(RateWindow)-1
%         blinkIndex_frag = flagConta(iF:iF+RateWindow); 
%         tempBlink_frag = diff(blinkIndex_frag);
%         if (blinkIndex_frag(1)==1) & (blinkIndex_frag(RateWindow+1)==1)
%             corrected.BlinkRate(iF+(RateWindow/2)+1) = (60/(RateWindow/500))*floor(sum(abs(tempBlink_frag))/2)+1; 
%         else
%             corrected.BlinkRate(iF+(RateWindow/2)+1) = (60/(RateWindow/500))*floor(sum(abs(tempBlink_frag))/2);
%         end
%     end

    % version 2: average considering past history ( -Window ~ 0 )
    corrected.BlinkRate = zeros(length(corrected.time),1); 
    for iF = 1:length(corrected.Pa)-(RateWindow)-1
        blinkIndex_frag = flagConta(iF:iF+RateWindow); 
        tempBlink_frag = diff(blinkIndex_frag);
        if (blinkIndex_frag(1)==1) & (blinkIndex_frag(RateWindow+1)==1)
            corrected.BlinkRate(iF+RateWindow+1) = (60/(RateWindow/500))*floor(sum(abs(tempBlink_frag))/2)+1; 
        else
            corrected.BlinkRate(iF+RateWindow+1) = (60/(RateWindow/500))*floor(sum(abs(tempBlink_frag))/2);
        end
    end
    
    tempBlinkconta = corrected.BlinkIndex; 
% Calculate microsaccade rate

    [resultMS, corrected] = extractMicrosaccades_mono (corrected);
    corrected.BlinkIndex = tempBlinkconta; 
    
    corrected.saccades = resultMS; 
%     RateWindow = 500*2; 
    corrected.MicrosaccadeRate = zeros(length(corrected.time),1); 
    for iF = 1:length(corrected.Pa)-(RateWindow)-1
        blinkIndex_frag = corrected.sacCandi3(iF:iF+RateWindow); 
        tempBlink_frag = diff(blinkIndex_frag);
        if (blinkIndex_frag(1)==1) & (blinkIndex_frag(RateWindow+1)==1)
            corrected.MicrosaccadeRate(iF+RateWindow+1) = (1/(RateWindow/500))*floor(sum(abs(tempBlink_frag))/2)+1; 
        else
            corrected.MicrosaccadeRate(iF+RateWindow+1) = (1/(RateWindow/500))*floor(sum(abs(tempBlink_frag))/2);
        end
    end
            
% % Figure plotting
%     figure(1000); 
%     subplot(5,2,Sbj_info(3)+(2*0)); 
%     plot(corrected.time, corrected.Pa_padded); 
%     xlim([min(corrected.time) max(corrected.time)]); 
%     if Sbj_info(3) == 1
%         title(['Subject #' num2str(Sbj_info(1)) ', iRun = ' num2str(Sbj_info(2)) ', Rivalry'],'fontsize',14,'fontweight','bold'); 
%     else
%         title(['Subject #' num2str(Sbj_info(1)) ', iRun = ' num2str(Sbj_info(2)) ', Replay'],'fontsize',14,'fontweight','bold'); 
%     end
%     
%     subplot(5,2,Sbj_info(3)+(2*1)); 
%     plot(corrected.time, corrected.Gx_padded);
%     xlim([min(corrected.time) max(corrected.time)]); 
% 
%     subplot(5,2,Sbj_info(3)+(2*2)); 
%     plot(corrected.time, corrected.BlinkRate);
%     xlim([min(corrected.time) max(corrected.time)]); 
%     
%     subplot(5,2,Sbj_info(3)+(2*3)); 
%     plot(corrected.time, corrected.MicrosaccadeRate);
%     xlim([min(corrected.time) max(corrected.time)]); 
%     
%     drawnow;
%     
% %     % Baseline time correction 
% %     corrected.time = corrected.time(baseNum:end)-corrected.time(baseNum); 
% %     corrected.Pa = corrected.Pa(baseNum:end); 
% %     corrected.Gx = corrected.Gx(baseNum:end); 
% %     corrected.Gy = corrected.Gy(baseNum:end); 
% %     corrected.BlinkIndex = corrected.BlinkIndex(baseNum:end); 
% %     corrected.BlinkRate = corrected.BlinkRate(baseNum:end); 
% %     corrected.Pa_padded = corrected.Pa_padded(baseNum:end); 
% %     corrected.Gx_padded = corrected.Gx_padded(baseNum:end); 
% %     corrected.Gy_padded = corrected.Gy_padded(baseNum:end); 
    
end



















