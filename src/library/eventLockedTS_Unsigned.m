% Event-related time-series cutting & stacking 
function eResult = eventLockedTS_Unsigned(TS, TS_time, event_time, event_type, sf, event_time_all, Windowlength)

% input: TS [1 X nT], TS_time [1 X nT], event_time [1 X nEvent], sf (sampling frequency)
% output: stack of event-locked TS [nEvent X nW]


% Define window length 
eResult.WindowLength_TS = ceil(Windowlength/sf); % in seconds

switch event_type
    case 'start' % transition start 
        eResult.WindowInd = (1:eResult.WindowLength_TS)-round(eResult.WindowLength_TS/4); 
    case 'mid'   % transition mid
        eResult.WindowInd = (1:eResult.WindowLength_TS)-round(eResult.WindowLength_TS/2);
    case 'end'   % transition end
        eResult.WindowInd = (1:eResult.WindowLength_TS)-round(3*eResult.WindowLength_TS/4);
end
eResult.Ind0point = find(eResult.WindowInd == 0); % event timing


eResult.stack_TS = [] ; 
eResult.valInd = [] ; 
eResult.stack_info = [] ; 
eResult.stack_TS_t = [] ;
for iTr = 1:length(event_time)
    
    % Convert absolute timing to TS index 
    tempEvent = find(abs(TS_time - event_time(iTr)) == min(abs(TS_time - event_time(iTr)))); 
    if length(tempEvent) >= 2 
        tempEvent = tempEvent(1); % take the smallest one
        event_time_TS(iTr) = tempEvent; 
    elseif length(tempEvent) == 0
        fprintf('no corresponding time points'); 
        event_time_TS(iTr) = nan; 
    else 
        event_time_TS(iTr) = tempEvent;
    end
    
    
    % Extract Time-series
    tempRange = eResult.WindowInd + event_time_TS(iTr); 

    if ~isnan(event_time_TS(iTr)) & (tempRange(1) > 0) & (tempRange(end) < length(TS))
        eResult.stack_TS = [eResult.stack_TS; TS(tempRange)];
        eResult.stack_info = [eResult.stack_info; event_time_all(iTr,:)]; 
        eResult.stack_TS_t = [eResult.stack_TS_t; TS_time(tempRange)]; 
        eResult.valInd = [eResult.valInd 1]; 
    else 
        eResult.valInd = [eResult.valInd 0]; 
    end
end
    
    
end




            
            
                
                