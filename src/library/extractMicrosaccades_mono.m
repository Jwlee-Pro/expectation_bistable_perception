function [saccades, eyeData] = extractMicrosaccades_mono (eyeData)

% assumes eyeData contains:
% high resolution gx, gy, pa

delta = 1/500; % 500 hz
veloThresh = 6; % of the median STD

lenData = length(eyeData.Gx);
axisTime = (.5:lenData) * delta; % secs

%% pupil-based  blink detection
% Di Stasi et al., 2013 ejn
% We identified and removed blink periods as portions of the raw data where
% pupil information was missing. We also removed portions of data where
% very fast decreases and increases in pupil area occrued
% (> 50 units/sample); such periods are probably semi-blinks during which
% the pupil is never fully occuluded (Troncoso et al., 2008; McCamy et al.,
% 2012). We added 200 ms before and after each blink or semi-blink to
% eliminate the initial and final parts during which the pupil was still
% partially occuluded (Troncoso et al., 2008; McCamy et al., 2012).


blinkPt = find(...
    (eyeData.Pa(:,1) < 100) | vertcat((abs(diff(eyeData.Pa(:,1))) > 50),0) );%| ... left eye
%     (eyeData.Pa(:,2) < 100) | vertcat((abs(diff(eyeData.Pa(:,2))) > 50),0)); % right eye

% add 200ms (100 samples) before/after each blink or semi-blink
eyeData.blinkConta = false(lenData,1);
for ii = 1:length(blinkPt)
    if blinkPt(ii) <= 100
        eyeData.blinkConta(1:blinkPt(ii)) = true;
    elseif blinkPt(ii) >= (length(eyeData.blinkConta)-100)
        eyeData.blinkConta(blinkPt(ii):end) = true;
    else
        eyeData.blinkConta(blinkPt(ii)+(-100:100)) = true;
    end
end

eyeData.blinkConta = eyeData.blinkConta | isnan(eyeData.Gx(:,1));
% eyeData.blinkConta = eyeData.blinkConta | isnan(eyeData.Gx(:,1)) | isnan(eyeData.Gx(:,2));


%% microsaccade detection
% Engbert & Kliegl, 2003; Engbert & Mergenthaler, 2006
% 1. the time series of eye positions was transformed to velocities by
% v_n = (x_n+2 + x_n+1 - x_n-1 - x_n-2) / 6*dt
% which represents a moving average of velocities over 5 data sample to
% suppress noise.
% 2. computation of velocity thresholds for the detection algoriths was
% based on the median of the velocity time series to protect the analysis
% from noise. A multiple of the standard deviation of the velocity
% distribution was used as the detection threshold.
% 3. microsaccades occur binocularly
% 4. combine microsaccades that are less than 20ms apart

eyeData.velx = single(zeros(size(eyeData.Gx)));
eyeData.vely = single(zeros(size(eyeData.Gy)));

% velocity & threshold computation
% v_n = (x_n+2 + x_n+1 - x_n-1 - x_n-2) / 6*dt
eyeData.velx = (circshift(eyeData.Gx,-2) + circshift(eyeData.Gx,-1) + ...
    - circshift(eyeData.Gx,1) - circshift(eyeData.Gx,2)) / (6*delta);
eyeData.velx([1 2 end-1 end],:) = 0;
eyeData.velx(eyeData.blinkConta,:) = 0;
eyeData.velxTh = sqrt( nanmedian(eyeData.velx.^2) - (nanmedian(eyeData.velx).^2) );

eyeData.vely = (circshift(eyeData.Gy,-2) + circshift(eyeData.Gy,-1) + ...
    - circshift(eyeData.Gy,1) - circshift(eyeData.Gy,2)) / (6*delta);
eyeData.vely([1 2 end-1 end],:) = 0;
eyeData.vely(eyeData.blinkConta,:) = 0;
eyeData.velyTh = sqrt( nanmedian(eyeData.vely.^2) - (nanmedian(eyeData.vely).^2) );

eyeData.vel2d = hypot(eyeData.velx, eyeData.vely);

% applying an elliptic threshold in 2d velocity space
eyeData.velInTh = zeros(size(eyeData.vel2d));
eyeData.ngx = zeros(size(eyeData.Gx));

eyeData.ngy = zeros(size(eyeData.Gy));
eyeData.velInTh(:,1) = hypot( ...
eyeData.velx(:,1)./eyeData.velxTh(1), ...
eyeData.vely(:,1)./eyeData.velyTh(1));
eyeData.ngx(:,1) = eyeData.Gx(:,1) - nanmedian( eyeData.Gx(:,1));
eyeData.ngy(:,1) = eyeData.Gy(:,1) - nanmedian( eyeData.Gy(:,1));


% binocularity check
eyeData.sacCandi = ...
    (eyeData.velInTh(:,1) > veloThresh);
    
%     and(eyeData.velInTh(:,1) > veloThresh, eyeData.velInTh(:,2) > veloThresh) & ... % they should occur in both
%     (eyeData.velx(:,1).*eyeData.velx(:,2) + ...
%     eyeData.vely(:,1).*eyeData.vely(:,2) > 0);  ... % moving to the same direction ... largely


% if saccade candidates occur within 30ms (10 samples), merge them
eyeData.sacCandi2 = false(size(eyeData.sacCandi));
candIdx = find(eyeData.sacCandi);
for iPt = 2:length(candIdx)
    if (candIdx(iPt) - candIdx(iPt-1)) <= 15
        eyeData.sacCandi2(candIdx(iPt-1):candIdx(iPt)) = true;
    else
        eyeData.sacCandi2(candIdx(iPt)) = true;
    end
end

% count the clusters & apply the duration filter: at least 3 samples
candIdx = find(eyeData.sacCandi2);
eyeData.sacCandi3 = false(size(eyeData.sacCandi));
for iPt = 1:length(candIdx)
    if iPt > 2
        flag1 = (candIdx(iPt)-candIdx(iPt-2) == 2);
    else flag1 = false;
    end
    if and(iPt > 1, iPt < length(candIdx))
        flag2 = (candIdx(iPt+1)-candIdx(iPt-1) == 2);
    else flag2 = false;
    end
    if iPt < length(candIdx)-2
        flag3 = (candIdx(iPt+2)-candIdx(iPt) == 2);
    else flag3 = false;
    end
    eyeData.sacCandi3(candIdx(iPt)) = (flag1 | flag2 | flag3);
end

% assign id to each saccades, then calculate the statistics
candIdx = find(eyeData.sacCandi3);
saccIdx = zeros(size(candIdx));
saccCnt = 1;
saccIdx(1) = 1;
for ii = 2:length(saccIdx)
    if diff(candIdx(ii+[-1 0])) < 2
        saccIdx(ii) = saccCnt;
    else
        saccCnt = saccCnt + 1;
        saccIdx(ii) = saccCnt;
    end
end
% currCnt is the number of the detected microsaccades
% 7+3 columns: start, duration, peak velocity, delta x, delta y, amplitude, phi
% 3 manual confirmations 

eyeData.saccInit = false(size(eyeData.sacCandi));
saccades = zeros(saccCnt, 10);

for iSac = 1:saccCnt
    saccRng = candIdx(saccIdx == iSac);
    eyeData.saccInit(min(saccRng)) = true; % tag like spike
    saccades(iSac, 1) = axisTime(min(saccRng)); % start time
    saccades(iSac, 2) = length(saccRng) * 2; % ms, duration
    saccades(iSac, 3) = mean(max(eyeData.vel2d(saccRng,:))); % peak velocity
    
    delXX = mean(eyeData.ngx(saccRng(end),:) - eyeData.ngx(saccRng(1),:)); 
    delYY = mean(eyeData.ngy(saccRng(end),:) - eyeData.ngy(saccRng(1),:)); 
    saccades(iSac, 4) = delXX; 
    saccades(iSac, 5) = delYY; 
    saccades(iSac, 6) = hypot(delXX, delYY); % amplitude
    saccades(iSac, 7) = atan2(delYY, delXX);
end

figStat = figure(448); clf;
set(gcf,'color','w');
set(figStat, 'Position', [500 500 700 340]);
subplot(1,2,1); 
set(gca, 'FontSize',14);
nn = hist(saccades(:,6), 0:.1:2); 
bar(0:.1:2, nn); hold on;
plot(median(saccades(:,6))*[1 1], [0 max(nn)+30], 'r-');
xlim([0 2]);
xlabel('Amplitude (deg)');
ylim([0 max(nn)+30]);
title(['Count: ' num2str(saccCnt) ', Freq: ' num2str((500*saccCnt)/lenData,3) ' hz']);
box off;

subplot(1,2,2);
rose(saccades(:, 7));
title('Saccade direction', 'FontSize',14);
box off;


% stat filter
%validSac = and(saccades(:,6) < 1, saccades(:,7) < 2);


% 
% %% microsaccades confirmation
% for iRep = 1:3
%     testIdx = randperm(saccCnt);
%     for ii = 1:saccCnt
%         curSacc = testIdx(ii);
%         saccRng = candIdx(saccIdx == curSacc);
% 
%         % display the figure with a key input
%         callstr = 'set(gcbf,''Userdata'',double(get(gcbf,''Currentcharacter''))) ; uiresume ' ;
%         fh = figure(...
%             'name','microsaccade manual check-up', ...
%             'keypressfcn',callstr, ...
%             'windowstyle','modal',...
%             'numbertitle','off', ...
%             'position',[200 200 1000 700],...
%             'userdata','timeout') ; clf;
%         
%         
%         % velocity threshold
%         subplot(3,3,[1 2]);
%         plot(saccades(curSacc, 1)*[1 1], [0 30], 'r-'); hold on;
%         plot(axisTime, eyeData.velInTh(:,1), 'b-');
% %         plot(axisTime, eyeData.velInTh(:,2), 'k-');
%         plot([0 max(axisTime)], veloThresh*[1 1], 'r--');
%         for iSac = 1:saccCnt
%             plot(saccades(iSac, 1)*[1 1], [0 30], 'r:');
%         end
%         xlim(saccades(curSacc, 1) + [-.7 .7]);
%         ylim([0 15]);
%         ylabel('Norm velocity', 'FontSize',14);
%         box off;
%         
%         subplot(3,3,3); % stats
%         text(0, 4, ['duration: ' num2str(saccades(curSacc, 2)) ' ms'], 'FontSize',16);
%         text(0, 3, ['peak velo: ' num2str(saccades(curSacc, 3),3) ' deg/s'], 'FontSize',16);
%         text(0, 2, ['amplitude: ' num2str(saccades(curSacc, 6),3) ' deg'], 'FontSize', 16);
%         text(0, 1, ['angle: ', num2str(180*saccades(curSacc, 7)/pi,3) ' deg'], 'FontSize', 16);
%         axis off;
%         xlim([-.5 5]);
%         ylim([0.5 4.5]);
%                 
%         % horizontal gaze
%         subplot(3,3,[4 5]);
%         plot(saccades(curSacc, 1)*[1 1], [-10 10], 'r-'); hold on;
%         plot(axisTime, eyeData.ngx(:,1), 'b-');
% %         plot(axisTime, eyeData.ngx(:,2), 'k-');
%         for iSac = 1:saccCnt
%             plot(saccades(iSac, 1)*[1 1], [-10 10], 'r:');
%         end
%         xlim(saccades(curSacc, 1) + [-.7 .7]);
%         ylim([min(min(eyeData.ngx(saccRng,:)))-.3  max(max(eyeData.ngx(saccRng,:)))+.3]);
%         ylabel('Horizontal (deg)', 'FontSize',14);
%         box off;
%         
%         subplot(3,3,6); % zoomed
%         plot(axisTime, eyeData.ngx(:,1), 'bo-'); hold on;
% %         plot(axisTime, eyeData.ngx(:,2), 'ko-');
%         plot(saccades(curSacc, 1)*[1 1], [-10 10], 'r-', 'LineWidth',2);
%         plot((saccades(curSacc, 1)+saccades(curSacc, 2)/1000)*[1 1], [-10 10], 'r-');
%         xlim(saccades(curSacc, 1) + [-.012 saccades(curSacc, 2)/1000+.012]);
%         ylim([min(min(eyeData.ngx(saccRng,:)))-.2  max(max(eyeData.ngx(saccRng,:)))+.2]);
%         ylabel('Horizontal (deg)', 'FontSize',14);
%         box off;
%                 
%         % vertical gaze
%         subplot(3,3,[7 8]);
%         plot(saccades(curSacc, 1)*[1 1], [-10 10], 'r-'); hold on;
%         plot(axisTime, eyeData.ngy(:,1), 'b-');
% %         plot(axisTime, eyeData.ngy(:,2), 'k-');
%         for iSac = 1:saccCnt
%             plot(saccades(iSac, 1)*[1 1], [-10 10], 'r:');
%         end
%         xlim(saccades(curSacc, 1) + [-.7 .7]);
%         ylim([min(min(eyeData.ngy(saccRng,:)))-.3  max(max(eyeData.ngy(saccRng,:)))+.3]);
%         ylabel('Vertical (deg)', 'FontSize',14);
%         box off;
%         
%         subplot(3,3,9); % zoomed
%         plot(axisTime, eyeData.ngy(:,1), 'bo-'); hold on;
% %         plot(axisTime, eyeData.ngy(:,2), 'ko-');
%         plot(saccades(curSacc, 1)*[1 1], [-10 10], 'r-', 'LineWidth',2);
%         plot((saccades(curSacc, 1)+saccades(curSacc, 2)/1000)*[1 1], [-10 10], 'r-');
%         xlim(saccades(curSacc, 1) + [-.012 saccades(curSacc, 2)/1000+.012]);
%         ylim([min(min(eyeData.ngy(saccRng,:)))-.2  max(max(eyeData.ngy(saccRng,:)))+.2]);
%         ylabel('Vertical (deg)', 'FontSize',14);
%         box off;
%         
%         while 1
%             uiwait;
%             resp = get(fh, 'Userdata');
%             if resp == 122 % 'Z' key in -> confirm
%                 saccades(curSacc,7+iRep) = 1;
%                 break;                
%             end
%             if resp == 120 % 'X' key in -> reject
%                 saccades(curSacc,7+iRep) = -1;
%                 break;
%             end
%         end
%         delete(fh); % move to next        
%     end
%     
% end
%     
%     
%     
    
    
end
