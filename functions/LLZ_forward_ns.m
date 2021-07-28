function LLZ = LLZ_forward_ns(SH, varargin)

%[LLZ]=LLZ_forward_ns(SH, 'to', 'ewh','resolution',1,'iGauss', 300,'iPxMy',0,'iplot',1);


PAR0 = struct( 'ran',[], 'grid','pole;;block', 'resolution',1, ...
    'iGauss', 0,'iPxMy',0,'iFan',0,'iddk',0,...
    'smooth', [],...
    'to','ewh',...
    'from', 0,...
    'iplot', 0, ...
    'icran',[], ...
    'PAR', [], ...
    'ex', 0);

PAR = var_initial(PAR0, varargin);

if PAR.ex == 1
    fun_ex1;
    return
elseif PAR.ex == 2
    fun_ex2;
    return
end

% SH = CS_smooth(SH,PAR);

field = SH.reformat('sc');

Nmonth = numel(SH);
for ii = 1:Nmonth
    field(:,:,ii) = unit_convert(field(:,:,ii),'none',PAR.to);
end

% tic
[V_pot, theRAD, lamRAD] = gshs_multi(field, 'quant', 'none', ...
    'grid', PAR.grid, 'gridsize', 180/PAR.resolution, 'height', 0, 'sub_wgs84', false);
% toc

% tic
% for ii = 1:Nmonth
%     [V_pot1(:,:,ii), theRAD, lamRAD] = gshs_(field(:,:,ii), 'quant', 'none', ...
%         'grid', PAR.grid, 'gridsize', 180/PAR.resolution, 'height', 0, 'sub_wgs84', false);
% end
% toc


% mydisp('maximum difference: %.3e\n',max(abs(V_pot1(:) - V_pot(:))));
%%
lat = 90-rad2deg(theRAD);

lon = rad2deg(lamRAD);
% remove roundoff error {
% lat = sig_digits(lat,'decim',10);
lat = round(lat,10);
lon = round(lon,10);
% lon = sig_digits(lon,'decim',10);
% remove roundoff error };
[mlon,mlat] =meshgrid(lon,lat);

tt = [SH.tt];
% LLZ = struct('lon',mlon, 'lat',mlat,'rg', V_pot, 'tt',tt);

if isfield(SH(1),'info') && isfield('name',SH(1).info)
    sname = sprintf('%s / %s',SH(1).info.name,PAR.to);
else
    sname = [];
end
if strcmp(PAR.to,'water')
    sunit = 'm';
else
    sunit = [];
end
LLZ = cLLZ(mlon,mlat,V_pot,tt);
LLZ = LLZ.setinfo('name',sname,'unit',sunit);

if ~isempty(PAR.ran) && numel(PAR.ran) == 4
    LLZ = LLZ.cut(PAR.ran);
end


end