function SH=LLZ_inverse_ns(LLZ, lmax, from, varargin)

PAR0 = struct('method','fnm'); %'grid','pole;;block',
PAR = var_initial(PAR0, varargin);

% LLZ = LLZ_extend_global(LLZ,'fill_value',0);

f = flipud(LLZ.rg);
% f = LLZ.rg;


ran6 = LLZ.ran;

if all(LLZ.ran4 == [0,360-ran6(2),-90,90])
    grd = 'pole';
elseif all(LLZ.ran4 == [ran6(2)/2,360-ran6(2)/2,-90+ran6(5)/2,90-ran6(5)/2])
    grd = 'block';
else
    fprintf('input is no supported (ran6=%.2f,%.2f,%.2f,%.2f,%.2f,%.2f)\n',ran6);
    error(1);
end

% grd = PAR.grid;

% lmax = 60;
for ii = 1:size(f,3)
    cs= gsha(f(:,:,ii), PAR.method, grd, lmax);

    sc(:,:,ii) = cs2sc(cs);
    sc(:,:,ii) = unit_convert(sc(:,:,ii),from,'none');
end


SH = cSH(sc,LLZ.tt,'sc');
SH = SH.setinfo('name',LLZ.info.name,'unit','none');
end