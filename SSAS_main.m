clear;
%% envirionment setting
addpath functions

%% read in data
% SH0 = load_SH('csr06','trange',[2004,12,2005,2]);
% SH0.save('csr06_gsm_2004-12_2005-02.mat');
a = load('data/csr06_gsm_2004-12_2005-02.mat'); % SH
SH0 = cSH(a.SH); % convert to the data format used here.
SH0 = SH0(3); % only use the first month

%% main parameters
% do not search coefficients with the degree <= Njump for NS stripe noise
Njump = 20; 

% window width for SSA
M = 40;

% only the leading K modes are searched. Since the noise may also exist in 
%  higher modes, the whole process will iteratively remove the paired modes
%  and search for new paired ones, until 1) a maximum Nstep is reached, 
%  or 2) no more paired modes are found
K = 20;

% -- other parameters
% the signals in greenland and antarctic are reduced to 10%, 
%   set to 1 if no extra weight is needed
ipolar_weight = 0.1;

% 1: create a gap region in sumatra 
igap_sumatra = 1;

%% main program
% SH0: monthly input of SH
% SH_filter: filtered SH
% LLZ_noise: identified NS-stripe noise by SSAS

[SH_filter, LLZ_noise]=fun_ssa_spatial_filter( SH0, ...
    'Njump',Njump, 'M',M, 'K',K, ...
    'ipolar_weight',ipolar_weight, 'igap_sumatra',igap_sumatra);

%% check the results
% compare your figures with these in the check/ folder. 


icheck = 1;
if icheck == 1 % map view
    figure('position',[1,1,1000,600]);
    imon = 1;
    SH_t0 = SH0(imon);
    SH_t1 = SH_filter(imon);
    cran = [-10,10];

    subplot(2,2,1)
    
    [LLZ]=LLZ_forward_ns(SH_t0, 'to', 'gr','resolution',1);
    mypcolor(LLZ.lon,LLZ.lat,LLZ.rg);
    caxis(cran);
    title('(a) Original map w/o filtering, in gravity disturbance (\mu Gal)')
    
    subplot(2,2,2)
    
    [LLZ]=LLZ_forward_ns(SH_t1, 'to', 'gr','resolution',1);
    mypcolor(LLZ.lon,LLZ.lat,LLZ.rg);
    caxis(cran);
    title('(b) After SSAS')
    
    subplot(2,2,3);
    dSH = SH_t0 - SH_t1;
    [LLZ]=LLZ_forward_ns(dSH, 'to', 'gr','resolution',1);
    mypcolor(LLZ.lon,LLZ.lat,LLZ.rg);
    caxis(cran);
    title('(c) a - b')
    
elseif icheck == 2 % Kaula curve
    figure;
    imon = 1;
    [y1,~] = SH0(imon).kaula;
    [y2,nn] = SH_filter(imon).kaula;
    semilogy(nn(2:end),y1(2:end));
    hold on;
    semilogy(nn(2:end),y2(2:end));
    hold off;
    legend('Original','Filtered');
    grid on
    xlabel('Degree');
    ylabel('Dimensionless')
    
end