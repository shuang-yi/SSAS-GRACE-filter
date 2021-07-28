function LLZ0 = mlonmlat(itype)
% {{LLZ}} {{create}} {{condition}}

% mlonmlat([lon1, lon2, lat1, lat2]);
% mlonmlat([lon1, dlon, lon2, lat1, dlat, lat2])
% mlonmlat('global1');
% mlonmlat([lon,lat]);

% 1. itype = [lon1,lon2,lat1,lat2], or [lon1, dlon,lon2,lat1, dlat, lat2]
if isnumeric(itype) && numel(itype) ~= 1
    if numel(itype) == 4 % [lon1,lon2,lat1,lat2]
        [mlon,mlat] = meshgrid(itype(1):itype(2),itype(3):itype(4));
    elseif numel(itype) == 6 % [lon1, dlon, lon2, lat1, dlat, lat2]
        [mlon,mlat] = meshgrid(itype(1):itype(2):itype(3),itype(4):itype(5):itype(6));
    elseif numel(itype) == 2 % [lon,lat]
        mlon = itype(1);
        mlat = itype(2);
    else
        error('size of input = %d\n', numel(itype));
    end
    LLZ0 = cLLZ(mlon,mlat,zeros(size(mlon)),NaN);
    return
end

% 2. itype = 'china'
if ~isnumeric(itype)
    if strcmp(itype,'global361')
        itype = 6;
    elseif strcmp(itype,'global1')
        itype = 7;
    elseif strcmp(itype,'global0.5')
        itype = 8;
    elseif strcmp(itype,'global0.25')
        itype = 10;
    elseif strcmp(itype,'global0.5c')
        itype = 9;
    elseif strcmp(itype,'global1c')
        itype = 4;
    elseif strcmp(itype,'global0.1')
        itype = 101;
    elseif strcmp(itype,'cn') || strcmp(itype,'china')
        itype = 11;
    elseif strcmp(itype,'cn-L')
        itype = 16;
    elseif strcmp(itype,'tibet')
        itype = 1006;
    elseif strcmp(itype,'tibet-L')
        itype = 1007;
    elseif strcmp(itype,'tohoku')
        itype = 1012;
    elseif strcmp(itype,'sumatra')
        itype = 1013;
    elseif strcmp(itype,'box1')
        itype = 30;
    elseif strcmp(itype,'longyang')
        itype = 1022;
    elseif strcmp(itype,'GSFC_mascon')
        itype = 201;
    else
        error('unrecgnized itype: %s\n',itype);
    end
end

% 3. itype = 1/2/3 ...
dlon = 1;
dlat = 1;
ijump = 0;
if itype == 0
    disp('1: global, -179.5:1:179.5 / -89.5:1:89.5');
    disp('2: global, -180:1:180 / -90:1:90');
    disp('3: global, -179:1:180 / -89:1:90');
    disp('4: global, 0.5:1:359.5 / -89.5:1:89.5');
    disp('5: global, 0.5:2:359.5 / -89.5:2:89.5');
    disp('6: global, 0:1:360 / -90:1:90');
    disp('7: global, 0:1:359 / -90:1:90');
    disp('8: global, 0:0.5:360 / -90:0.5:90');
    disp('9: global, 0.25:0.5:359.75 / -89.75:0.5:89.75');
    disp('10: global, 0:0.25:360 / -90:0.25:90');
    disp('101: global, 0:0.1:360 / -90:0.1:90');
    disp('11: china, 65:1:140 / 15:1:55');
    disp('12: china, 65:0.5:140 / 15:0.5:55');
    disp('13: china, 65:0.2:140 / 15:0.2:55');
    disp('14: china, 65:0.1:140 / 15:0.1:55');
    disp('15: china, 65.5:1:140.5 / 15.5:1:55.5');
    disp('16: china-L, 60:1:145 / 10:1:60');
    disp('21: greenland, -90:1:-10 / 55:1:89');
    disp('30: box1, 0:1:10 / 0:1:10');
    disp('31: box1, 0:1:20 / 0:1:20');
    disp('201: GSFC mascon');
    disp('1001: yaluzangbu, 80:1:105 / 20:1:40');
    disp('1002: lhasa, 90:1:92 / 28:1:30');
    disp('1003: india-negative, 79:1:81 / 24:1:26');
    disp('1004: songhuajiang, 115:1:135 / 40:1:55');
    disp('1005: tibet-lake, 78:1:92 / 28:1:38');
    disp('1006: tibet, 70:1:105 / 26:1:40');
    disp('1007: tibet-larger, 65:1:110 / 15:1:50');
    disp('1008: tibet-2, 65:1:105 / 20:1:50');
    disp('1010: asia, 20:1:150 / 0:1:60');
    disp('1011: point in West kunlun, 81.037 / 35.392');
    disp('1012: Tohoku earthquake 2011, 130:0.5:154 / 26:0.5:50');
    disp('1013: Sumatra earthquake 2004, 84:0.5:104 / -6:0.5:14');
    disp('1014: tibet east, 95:2:105 / 28:2:38');
    disp('1020~1022: sanxia, longtanxia, longyangxia');
    disp('1030 / 1031: tianshan / tianshan negative center');
    disp('1040 : wang remote sensing, himalya + tianshan 0.5');
    disp('1041 : larger version of 1040, for Yinzhi-WHU');
    LLZ0 = NaN;
    return
elseif itype == 1
    ran = [-179.5,179.5, -89.5,89.5];
    
elseif itype == 2
    ran = [-180,180, -90,90];
    
elseif itype == 3
    ran = [-179,180, -89,90];
    
elseif itype == 4
    ran = [0.5,359.5, -89.5,89.5];
    
elseif itype == 5
    ran = [-179, 179, -89, 89];
    dlon = 2;
    dlat = 2;
    
elseif itype == 6
    ran = [0, 360, -90, 90];
    
elseif itype == 7
    ran = [0, 359, -90, 90];
    
elseif itype == 8
    ran = [0, 360, -90, 90];
    dlon = 0.5;
    dlat = 0.5;
    
elseif itype == 9
    ran = [0.25, 359.75, -89.75, 89.75];
    dlon = 0.5;
    dlat = 0.5;
elseif itype == 10
    ran = [0, 360, -90, 90];
    dlon = 0.25;
    dlat = 0.25;
elseif itype == 101
    ran = [0, 360, -90, 90];
    dlon = 0.1;
    dlat = 0.1;
elseif itype == 11  %china
    ran = [65,140,15,55];
    
elseif itype == 12  %china
    ran = [65,140,15,55];
    dlon = 0.5;
    dlat = 0.5;
    
elseif itype == 13  %china
    ran = [65,140,15,55];
    dlon = 0.2;
    dlat = 0.2;
    
elseif itype == 14  %china
    ran = [65,140,15,55];
    dlon = 0.1;
    dlat = 0.1;
    
elseif itype == 15  %china
    ran = [65,140,15,55]+0.5;
elseif itype == 16  %china-L
    ran = [60,145,10,60];
elseif itype == 21
    ran = [-90, -10, 55, 89];
    
elseif itype == 30
    ran = [0, 10, 0, 10];
    
elseif itype == 31
    ran = [0, 21, 0, 21];
elseif itype == 201 % GSFC mascon
    mypath('GSFC_mascon');
    LLZ0 = fun_gsfc_mascon('lonlat');
    
    ijump = 1;
    
elseif itype == 1001  %yaluzangbu
    ran = [80, 105, 20, 40];
    
elseif itype == 1002  %lhasa
    ran = [90, 92, 28, 30];
    
elseif itype == 1003  %india-negative
    ran = [79, 81, 24, 26];
    
elseif itype == 1004  %songhuajiang
    ran = [115, 135, 40, 55];
    
elseif itype == 1005  %tibet-lake
    ran = [78, 92, 28, 38];
    
elseif itype == 1006  %tibet
    ran = [70, 105, 26, 40];
    
elseif itype == 1007  %tibet-larger
    ran = [65, 110, 15, 50];
    
elseif itype == 1008  %tibet-M
    ran = [65, 105, 20, 50];
    
elseif itype == 1010  %asia
    ran = [20, 150, 0, 60];
    
elseif itype == 1011 % point in West kunlun,
    ran = [81.037, 81.037, 35.392, 35.392];
    
elseif itype == 1012  %Tohoku earthquake 2011
    ran = [130, 154, 26, 50];
    dlon = 0.5;
    dlat = 0.5;
    
elseif itype == 1013  %Sumatra earthquake 2004
    ran = [84, 104, -6, 14];
    dlon = 0.5;
    dlat = 0.5;
elseif itype == 1014  %east tibet
    %     ran = [95, 106, 27, 38]; % S
    %     ran = [95, 104, 25, 40]; % (S+M)/2
    ran = [94, 103, 28, 39]; % M
    %     ran = [90, 110, 23, 45]; % L
    dlon = 1;
    dlat = 1;
elseif itype == 1020  %sanxia
    addpath /Users/Yish/workplace/mypaper/P19_human-reserviors-in-China/work/function/ %load_data_P19.m
    
    range = load_data_P19('ran_reserviors');
    [ran]=range(1,:);
elseif itype == 1021  %sanxia
    addpath /Users/Yish/workplace/mypaper/P19_human-reserviors-in-China/work/function/ %load_data_P19.m
    
    range = load_data_P19('ran_reserviors');
    [ran]=range(2,:);
elseif itype == 1022  %longyangxia
    
    addpath /Users/Yish/workplace/mypaper/P19_human-reserviors-in-China/work/function/ %load_data_P19.m
    
    range = load_data_P19('ran_reserviors',10);
    [ran]=range(3,:);
elseif itype == 1030 % tianshan
    ran = [72, 95, 39, 49];
elseif itype == 1031 % tianshan, negative center
    ran = [82, 87, 42, 45];
elseif itype == 1032 % tianshan
    ran = [72, 95, 39, 49];
    dlon = 0.5;
    dlat = 0.5;
elseif itype == 1033 % tianshan
    ran = [65, 100, 34, 54];
elseif itype == 1033 % tianshan
    ran = [65, 100, 34, 54];
elseif itype == 1040 % wang, remote sensing, himalaya
    ran = [65, 100, 22, 45];
    %     ran = [75, 100, 26, 34];
    dlon = 0.5;
    dlat = 0.5;
elseif itype == 1041 % wang, remote sensing, himalaya
    ran = [65, 110, 20, 50];
    dlon = 0.5;
    dlat = 0.5;
else
    error('unrecgonized itype in mlonmlat');
end

if ijump == 0
    ran2(1:3) = [ran(1),dlon,ran(2)];
    ran2(4:6) = [ran(3),dlat,ran(4)];
    LLZ0 = mlonmlat(ran2);
end
end