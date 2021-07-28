function [tt] = read_time(filename,stype)
% {{time}} {{read}}
%{
sfile : a short name for filename

--- stype
1. 'csr':   20[012]\d{4}: 200xxxx or 201xxxx or 202xxxx, x is a number (0~9)
GSM-2_2006001-2006031_0031_UTCSR_0096_0005.gfc : csr
GSM-2_2006001-2006031_0031_EIGEN_G---_005a.gfc : gfz
kfilter_DDK1_GSM-2_2006032-2006059_0028_UTCSR_0096_0005 : DDK
GX-OG-_2-GSM+CSR-GSM-2_2002095-2002120_0021_UTCSR_0060_0005 : grcof2

2. 'daily':     20[01]\d-\d\d-\d\d: 200x-xx-xx
ITSG-Grace2014_2013-08-28.gfc

3. 'ymd_ymd': 20[01]\d_\d\d_\d\d_20[01]\d_\d\d_\d\d
2002_08_08_2002_08_17

4. 'yyyymm': (19|20)\d\d[01]\d
200101, 197503

5. 'yyyy-mm': (19|20)\d\d-[01]\d
2001-01

5.2. 'yyyy_mm': (19|20)\d\d-[01]\d
2001_01

6. 'yyyymmdd': (19|20)\d\d[01]\d[0123]\d
200101, 197503

7. 'swarm_ICGEM': 20\d{6}T\d{6}_20\d{6}T\{6}
20131201T000000_20131231T235959_0101
--- stype
%}



if strcmp(stype,'csr')
    pattern = '20[012]\d{4}';
    loc_time = regexp(filename,pattern);
    if numel(loc_time) ~= 2
        mydisp('Warning: cann''t read "%s", NaN is used\n',filename);
        tt = NaN;
    else    
        [tt] = for_csr(filename,loc_time);
    end
elseif strcmp(stype,'daily')
    pattern = '20[012]\d-\d\d-\d\d';
    loc_time = regexp(filename,pattern);
    if numel(loc_time) ~= 1
        error('Error: cann''t read %s\n',filename);
    end    
    [tt] = for_daily(filename,loc_time);
elseif strcmp(stype,'ymd_ymd')
    pattern = '20[012]\d_\d\d_\d\d';
    loc_time = regexp(filename,pattern);
    if numel(loc_time) ~= 2
        error('Error: cann''t read %s\n',filename);
    end    
    for ii = 1:2
        [tt0(ii)] = for_daily(filename,loc_time(ii));        
    end
    tt = mean(tt0);
elseif strcmp(stype,'yyyymm')
    pattern = '(18|19|20)\d\d[01]\d';
    [s,e] = regexp(filename,pattern);
    for ii = 1:numel(s)
        tt0 = sscanf(filename(s(ii):e(ii)),'%4d%2d');
        tt(ii) = tt0(1) + (tt0(2)-0.5)/12;
    end
elseif strcmp(stype,'yyyy-mm')
    pattern = '(18|19|20)\d\d-[01]\d';
    [s,e] = regexp(filename,pattern);
    for ii = 1:numel(s)
        tt0 = sscanf(filename(s(ii):e(ii)),'%4d-%2d');
        tt(ii) = tt0(1) + (tt0(2)-0.5)/12;
    end
elseif strcmp(stype,'yyyy_mm')
    pattern = '(18|19|20)\d\d_[01]\d';
    [s,e] = regexp(filename,pattern);
    for ii = 1:numel(s)
        tt0 = sscanf(filename(s(ii):e(ii)),'%4d_%2d');
        tt(ii) = tt0(1) + (tt0(2)-0.5)/12;
    end
elseif strcmp(stype,'yyyymmdd')
    pattern = '(18|19|20)\d\d[01]\d[0123]\d';
    [s,e] = regexp(filename,pattern);
    for ii = 1:numel(s)
        tt0 = sscanf(filename(s(ii):e(ii)),'%4d%2d%2d'); tt0 = tt0(:)';
        tt(ii) = time_transfer(tt0,1);
    end
elseif strcmp(stype,'swarm_ICGEM')
    % SW_OPER_EGF_SHA_2__20131201T000000_20131231T235959_0101
    pattern = '20\d{6}T\d{6}';
    [s,e] = regexp(filename,pattern);
    if numel(s) ~= 2
        error('input not recognized, should have 2 matched\n')
    end
    for ii = 1:2
        tt0 = sscanf(filename(s(ii):e(ii)),'%4d%2d%2d'); tt0 = tt0(:)';
        tt1(ii) = time_transfer(tt0,1);
    end
    tt = mean(tt1);
else
    error('%s %s\n','Wrong stype: ',stype);
end

end

% csr
function [tt] = for_csr(filename,loc_time)
% disp(sfile);
% example : 2006001-2006031_0031
a=sscanf(filename(loc_time(1):loc_time(1)+6),'%4d%3d');
yr1 = a(1);
day1 = a(2);
a=sscanf(filename(loc_time(2):loc_time(2)+6),'%4d%3d');
yr2 = a(1);
day2 = a(2);
tt = yr1 + ( day2+day1 + 365*(yr2-yr1) )/365/2;
end

% daily
function [tt] = for_daily(filename,loc_time)
% disp(sfile);
% example : 2013-08-28
a=sscanf(filename(loc_time:loc_time+9),'%4d%1s%2d%1s%2d');

yr = a(1);
mon = a(3);
day = a(5);
tt = time_transfer([yr,mon,day],1);
end

%{

%  example :   GAC-2_2003001-2003031_0025_UTCSR_0000_0005
a = sscanf(files.name{ifile},'%*6s%4d%3d%*1s%4d%3d'); %[y1,d1,y2,d2]
b = time_transfer([a(1:2)]);
b2 = time_transfer([a(3:4)]);
tt(ifile) = mean([b{1},b2{1}]);

%}