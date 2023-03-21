function [Tr_all, Tr_t_all, NtotalTr] = Stack_unsigned(eResult, Tr_all, Tr_t_all, NtotalTr, valInd, behavior)
    % Stack every fragments (this part is for neural data only; sampling frequency is the same)
	if ~isempty(eResult.stack_TS)
        if length(valInd) == length(eResult.stack_TS(:,1))
            val = find(valInd==1); 
        else
            val = 1:length(eResult.stack_TS(:,1));
        end
		Tr_all = [Tr_all; eResult.stack_TS(val,:)];
		NtotalTr = length(eResult.stack_TS(val,1)) + NtotalTr;
    end
    
    if (behavior == 1) & isempty(Tr_t_all)
		Tr_t_all = [Tr_t_all; eResult.stack_TS_t(val,:)];
    end
end