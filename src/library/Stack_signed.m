function [Tr_all, Tr_t_all, NtotalTr] = Stack_signed(eResult, Tr_all, Tr_t_all, NtotalTr, LRInd, valInd, behavior)
    % Stack every fragments (this part is for neural data only; sampling frequency is the same)
	if ~isempty(eResult.stack_TS)
        LRInd = LRInd(find(valInd==1)); 
        for iTr = 1:length(LRInd)
            if LRInd(iTr) == 1 % L > R
                Tr_all = [Tr_all; eResult.stack_TS(iTr,:)];
                Tr_t_all = [Tr_t_all; eResult.stack_TS_t(iTr,:)];
                NtotalTr = 1 + NtotalTr;

            elseif LRInd(iTr) == 0 % R > L
                Tr_all = [Tr_all; -eResult.stack_TS(iTr,:)];
                Tr_t_all = [Tr_t_all; eResult.stack_TS_t(iTr,:)];
                NtotalTr = 1 + NtotalTr;
            end
        end
    end
end