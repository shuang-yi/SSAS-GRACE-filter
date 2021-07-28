function LLZ1 = LLZ_append(LLZ1,LLZ2)
% {{LLZ}} {{append}}
% no check in the size or time
if isempty(LLZ1)
    if iscell(LLZ2)
        LLZ1 = LLZ2{1};
        for ii = 2:numel(LLZ2)
            LLZ1 = LLZ_append(LLZ1,LLZ2{ii});
        end
    else
        LLZ1 = LLZ2;
    end
else
    i1 = LLZ1.Nmonth+1;
    if iscell(LLZ2)
        for ii = 1:numel(LLZ2)
            LLZ1 = LLZ_append(LLZ1,LLZ2{ii});
        end
    else
        i2 = LLZ2.Nmonth+i1 - 1;
        LLZ1.rg(:,:,i1:i2) = LLZ2.rg;
        LLZ1.tt(i1:i2) = LLZ2.tt;
    end
end

LLZ1 = LLZ1.toStandard;
end