function [MPR, ReversalRate] = GenDominanceReversalRate(vectEvent_signed)

% analysis of dominance duration
ccnt=0;
MPR=[]; % Markov Renewel Process

for iT=1:length(vectEvent_signed)
	if(iT==1)
		vectEvent_cum(iT)=vectEvent_signed(iT);
% 		ccnt=ccnt+1;
% 		MPR(ccnt,:)=[iT,vectEvent_cum(iT)];
	elseif (1 < iT) & (iT < length(vectEvent_signed))
		if(vectEvent_signed(iT)==vectEvent_signed(iT-1))
			vectEvent_cum(iT)=vectEvent_cum(iT-1)+vectEvent_signed(iT);
		else
			vectEvent_cum(iT)=vectEvent_signed(iT);
			ccnt=ccnt+1;
			MPR(ccnt,:)=[iT-1,vectEvent_cum(iT-1)];
		end
	elseif (iT==length(vectEvent_signed))
		if(vectEvent_signed(iT)==vectEvent_signed(iT-1))
			vectEvent_cum(iT)=vectEvent_cum(iT-1)+vectEvent_signed(iT);
		else
			vectEvent_cum(iT)=vectEvent_signed(iT);
		end
		ccnt=ccnt+1;
		MPR(ccnt,:)=[iT-1,vectEvent_cum(iT)];
	end
end


% Computing Reversal Rate
numReversal=sum(abs(diff(vectEvent_signed)))/2;
durRun_sec=300;
% reversal per min
ReversalRate = numReversal/(durRun_sec/60);


% % power-law distribution
% [valx(1), valx(2), valx(3)]=plfit(abs(MPR(:,2))); 
% powerLawProperty = valx ; % [alpha, xmin, L]



