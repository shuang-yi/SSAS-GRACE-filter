function [tt,sfilename] = generate_tt(t1,o_t2,idatenum)
% {{time}} {{create}} {{str}}
% t1 : start time, [year,month]
% t2 : end time
% input is ([y1,m1], [y2,m2]) or ([y1,m1, y2,m2]);

if nargin == 1
    if numel(t1) == 4
        o_t2 = t1(3:4);
    else
        o_t2 = t1;
    end
end

if nargin < 3
    idatenum = 0;
end

if numel(t1) == 1 
    t1 = time_transfer(t1,3);
    o_t2 = time_transfer(o_t2,3);
end


iyear1 = t1(1);
iyear2 = o_t2(1);
imonth1 = t1(2);
imonth2 = o_t2(2);

% Ntt = (iyear2-iyear1)*12 - imonth1 + imonth2 +1;
NN = (iyear2-iyear1+1)*12;
rtmp(1:NN,1) = 0;

ik = 0;
for iyear = iyear1:iyear2
    for imonth = 1:12
        ik = ik+1;
        rtmp(ik) = (imonth-0.5)/12 + iyear;% - iyear1;
    end
end

tt = rtmp(imonth1:end-12+imonth2);

if idatenum == 1
    tt = time_transfer(tt(:),-1);
end

% ============================================================
ik = 0;
for iyear = iyear1:iyear2
    for imonth = 1:12
        ik = ik+1;
        stmp{ik} = sprintf('%6.6d',iyear*100 + imonth);
    end
end
% stmp{1:2}
sfilename = stmp(imonth1:end-12+imonth2);


end