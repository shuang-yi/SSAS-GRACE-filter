function [SH_filter, LLZ_noise, PAR]=fun_ssa_spatial_filter(SH, varargin)
% {{SSAS; spatial; filter; GRACE}}

PAR0 = struct('Njump',0,'M',40,'K',20,'iplot',0,'corr',0.9,...
    'imode','global;;regional','ran','global1', 'only_pair',1, ...
    'ilat_list',[],'max_shift',3,'Nstep',6,...
    'ilist',[],'igap_sumatra',1,'ind_nan',[],'ipolar_weight',0.1,'istep_detail',0,...
    'ito','gr','grid','block;;pole');

% corr=0.9: If two modes have a correlation value >=0.9, then they are paired
% max_shift: Two modes will be shifted forward/backward by max_shift to test their correlation
% Nstep: the process will iteratively search for coupled modes and mark them
%        as NS-stripe noise, but the maximum step is Nstep
% igap_sumatra=1: signals in Sumatra are filled by a gap-filling method
%             =0: no gap filling method will be applied
% ipolar_weight=0.1: extreme signals in Greenland/Antartica are weakened to 10%
%              =1:  no extra weight is applied

PAR = var_initial(PAR0, varargin);

ran = [0,360,-90,90];

% assign zeros values for coefficients with the degree < Njump
SH_cut = SH.setval(sprintf('n<=%d',PAR.Njump+1),0,0); 

% convert SH to LLZ
LLZ = LLZ_forward_ns(SH_cut, 'to', PAR.ito,'ran',ran,'grid',PAR.grid);

% ----------ind_nan ------------- {
if PAR.igap_sumatra == 1
    a = load('data/LLZ_mask_three_separate.mat'); 
    LLZ_mask2 = a.LLZ(1);
    ind = LLZ_mask2.rg>0.6;
    
    PAR.ind_nan = ind;
end
% ----------ind_nan ------------- }

% original observations, full degrees, no filter
LLZ_ori = LLZ_forward_ns(SH, 'to', PAR.ito,'ran',ran,'grid',PAR.grid);

LLZ_noise = [];

for ii = 1:SH.Nmonth
    fprintf('Processing monthly product %d\n',ii);
    
    [~,LLZ_noise0]=fun_ssa_spatial_whole_3deg(LLZ.extract(ii),PAR);
    
    LLZ_signal0 = LLZ_ori.extract(ii) - LLZ_noise0;
     
    SH_filter(ii) = LLZ_inverse_ns(LLZ_signal0,SH.Nmax, PAR.ito);
        
    LLZ_noise = LLZ_append(LLZ_noise,LLZ_noise0);
    
end

end