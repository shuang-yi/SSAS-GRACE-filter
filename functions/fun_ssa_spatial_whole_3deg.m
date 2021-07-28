function [LLZ_signal,LLZ_noise]=fun_ssa_spatial_whole_3deg(LLZ,PAR)
%
%

% PAR0 = struct('ex',0,'M',120,'K',100,'corr',0.9,'only_pair',1);
% PAR = var_initial(PAR0, varargin);

% --------------------end of head-----------

% LLZ0 = LLZ;
mlon = LLZ.lon; lon_one = mlon(1,:); lon_one = lon_one(:);
mlat = LLZ.lat; lat_one = mlat(:,1); lat_one = lat_one(:);
Nlon = numel(lon_one);
Nlat = numel(lat_one);

ilat_ran = 2:3:Nlat;
ilon_ran = 2:3:Nlon;
%     M = PAR.M;
K = PAR.K;
% rweight=fun_weight(mlat(ilat_ran,:)', 1);
rweight = sqrt(cosd(mlat(ilat_ran,:)'));
rweight(rweight<0.01) = 0.01;

% rweight = ones(size(rweight)); fprintf('uniform weight\n');

ind_nonzero = rweight(1,:)~=0;

if PAR.ipolar_weight == 1
    a = load('data/LLZ_mask_three_separate.mat');
    LLZ_mask2 = a.LLZ;
    rweight2 = LLZ_mask2(2).rg+LLZ_mask2(3).rg;
    rweight2 = 1-rweight2*(1-PAR.ipolar_weight);
    LLZ.rg = LLZ.rg.*rweight2;
end

for istep = 1:PAR.Nstep
    rg = LLZ.rg;
%     rg_signal = rg*NaN;
%     rg_noise = rg*NaN;
rg_noise = rg(ilat_ran,:)*nan;
    
    
    X = rg(ilat_ran,:)';
    
    X_mean = mean(X);
    X = X - X_mean;    

    X = X.*rweight;
    
    X = X(:,ind_nonzero);
    % ------- gap filling --------{
    if ~isempty(PAR.ind_nan)
        if istep == 1
            PAR.ind_nan = PAR.ind_nan';
            PAR.ind_nan = PAR.ind_nan(:,ilat_ran(ind_nonzero));
        end
        
        if istep == 3
            ind_save = PAR.ind_nan;
            PAR.ind_nan = [];
            X(ind_save) = 0;
        end
        
        [LAMBDA, RC, EOF, PC, ilist] = mssa_longitude_missing(X, PAR);
        
        % ------- gap filling --------}
    else
        if exist('ind_save','var') && istep >=3
            X(ind_save) = 0;
        end
        [LAMBDA, RC, EOF, PC, ilist] = mssa_longitude(X, PAR);
    end
    
    if isempty(ilist) && istep > 1
        break;
    end
%     [ilist,inot]=paired_PC(PC,PAR.corr, PAR.max_shift);

    iplot = 0;
    if iplot == 1
        myfigure;
        for ii = 1:K
            if min(abs(ii-ilist)) == 0
                color = 'r';
            else
                color = 'k';
            end
            plot(PC(:,ii)-ii*0.2,'color',color);
            hold on;
        end
        hold off;
    elseif iplot == 2
        myfigure;
        for ii = 1:K
            ilat = 30;
            if min(abs(ii-ilist)) == 0
                color = 'r';
            else
                color = 'k';
            end
            plot(RC(:,ii,ilat)-ii*1,'color',color);
            hold on;
        end
        plot( sum(RC(:,inot,ilat),2) ,'linew',2)
        hold off;
    end
    
    if ~isempty(PAR.ilist)
        ilist = PAR.ilist;
        mydisp('ilist is pre-set\n');
    end
    
    iloc = 1:sum(ind_nonzero);
    iloc = iloc(ind_nonzero);
%     ilat_ran2 = ilat_ran(ind_nonzero);
    for ij = iloc
        y = sum(RC(:,ilist,ij),2);
        y = y';
        rg_noise(ij,:) = y;
        
        %     rg_signal(ij,:) = y+X_mean(ij);
        %     rg_noise(ij,:) = rg(ij,:)-rg_signal(ij,:);
    end
    rg_noise = rg_noise./rweight';
    rg_noise(~isfinite(rg_noise)) = 0;
    LLZ_noise0 = cLLZ(LLZ.lon(ilat_ran,ilon_ran),LLZ.lat(ilat_ran,ilon_ran),rg_noise(:,ilon_ran),LLZ.tt);
    SH = LLZ_inverse_ns(LLZ_noise0,60,'none','method','ls');
    LLZ_noise = LLZ_forward_ns(SH, 'to', 'none','grid','block');
    rg_noise = LLZ_noise.rg;
    
    rg_signal = rg - rg_noise;
    
    LLZ_signal = cLLZ(mlon,mlat,rg_signal,LLZ.tt);
%     LLZ_noise = cLLZ(mlon,mlat,rg_noise,LLZ.tt);
    if istep == 1
        LLZ_noise_save = LLZ_noise;
    else
        LLZ_noise_save = LLZ_noise + LLZ_noise_save;
    end
    LLZ = LLZ_signal;
    
end

LLZ_noise = LLZ_noise_save;

end
