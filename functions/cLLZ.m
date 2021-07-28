classdef cLLZ
    % 1. create
    %  LLZ = cLLZ(LLZ0)
    %  LLZ = cLLZ(lon,lat,rg,tt)
    %  LLZ = cLLZ([lon(:),lat(:),rg(:)], N1, N2, tt)
    %  LLZ = cLLZ(N_lonlat);
    %  LLZ = cLLZ(LLZ_old, rg_new, tt_new)
    % 2. info
    %  addinfo(msg); dispinfo(imon); displonlat(Nline)
    %  Nlon; Nlat; vlat; vlon; Nmonth; ran; ran2; trange
    % 3. check
    %  irr; compatible(LLZ1,LLZ2);
    % 4. calculation
    %  minus, plus, mean
    % 5. reshape
    %  cut(ran); extract(ilist_or_trange); LLZ2std
    
    % {{LLZ; class}}
    
    properties
        lon
        lat
        rg
        tt
        info % type(normal/single/vector/irr); unit; name; history
    end
    
    methods
        % create {
        function obj = cLLZ(varargin)
            if nargin == 4 % two conditions
                if numel(varargin{2})>1 || numel(varargin{1}) == 1  % lon,lat,rgs ,tt
                    % numel(varargin{1}) == 1: only one grid
                    [lon,lat,rgs,tt] = deal(varargin{1:4});
                    obj.rg = rgs;
                    obj.lon = lon;
                    obj.lat = lat;
                    obj.tt = tt(:);
                    obj.info = datainfo;
                    
                else % [lon(:),lat(:),rg(:)], N1, N2, tt
                    [mat,N1,N2] = deal(varargin{1:3});
                    lon = reshape(mat(:,1),N1,N2);
                    lat = reshape(mat(:,2),N1,N2);
                    rg = reshape(mat(:,3),N1,N2);
                    obj = cLLZ(lon,lat,rg,varargin{4});
                end
            elseif nargin == 1
                if (isstruct(varargin{1}) && isfield(varargin{1},'lon')) || isa(varargin{1},'cLLZ')
                    % LLZ ['lon','lat','rg','tt'], or cLLZ
                    LLZ = varargin{1};
                    rg = LLZ.rg;
                    lon = LLZ.lon;
                    lat = LLZ.lat;
                    try
                        tt = LLZ.tt(:);
                    catch
                        tt = NaN;
                    end
                    
                    obj  = cLLZ(lon,lat,rg,tt);
                    
                    try % inherit the info, but it may fail
                        obj = obj.setinfo(LLZ.info);
                    catch
                        %                         warning('info is not read in correctly');
                    end
                else % Nlonlat
                    try
                        OBS = mlonmlat(varargin{1});
                        obj = cLLZ(OBS.lon,OBS.lat,zeros(size(OBS.lon)),NaN);
                    catch
                        error('input is not supported\n');
                    end
                    
                end
            elseif nargin == 3 % cLLZ, rg, tt
                obj = varargin{1};
                obj.rg = varargin{2};
                tt = varargin{3};
                obj.tt = tt(:);
            else
                error('unsupported input');
            end
            obj = obj.toStandard;
            
        end
        % create }
        
        %% info {
        % multiple msgs should be stored in cell
        function obj = setinfo(obj,varargin)
            % LLZ.setinfo('name','xx','unit','xx','type','xx','history',{'xx','yy'});
            % LLZ.setinfo(struct('name','xx'));
            % LLZ.setinfo(LLZ2);
            %
            % multiple history should be stored in cell
            
            % {{cLLZ; set; info}}
            
            for ii = 1:numel(obj)
                rtmp = obj(ii).info;
                if ~isa(rtmp,'datainfo') % in case the info type is the old version
                    rtmp = datainfo(rtmp);
                end
                
                if isa(varargin{1},'cLLZ') % input is another cLLZ
                    a = varargin{1};
                    if numel(a) == 1  % asign one to multiple
                        obj(ii).info = rtmp.set(a.info);
                    elseif numel(a) == numel(obj) % asign multiple to multiple
                        obj(ii).info = rtmp.set(a(ii).info);
                    end
                else % input is a structure or duple parameters
                    obj(ii).info = rtmp.set(varargin{:});
                end
            end
            
        end
        
        function print(obj, imon)
            % {{cLLZ; dispinfo}}
            
            if nargin == 1
                imon = 1;
            end
            
            fprintf('===== info =====\n');
            
            % if obj.Nmonth>1
            %     fprintf('more the one input, only the %d is used\n',imon);
            %     obj0 = obj.extract(imon);
            % end
            
            obj.info.print;
            
            
            fprintf(' - ran: %g, %g, %g,  %g, %g, %g\n',obj.ran);
            
            if obj.Nmonth > 1
                fprintf(' - size of rg: [%d, %d, %d]\n',size(obj.rg));
            else
                fprintf(' - size of rg: [%d, %d, 1]\n',size(obj.rg));
            end
            
            fprintf(' - range of rg in month %d:\n',imon);
            rg = obj.rg(:,:,imon); rg = rg(:);
            ind = ~isnan(rg);
            
            fprintf('   (min) %g; (mean) %g; (median) %g; (max) %g; (number_of_nan) %g\n',...
                min(rg(ind)), mean(rg(ind)), median(rg(ind)), max(rg(ind)), sum(~ind))
            
            fprintf(' - trange:')
            if all(isnan(obj.tt))
                fprintf(' %d * NaN\n',numel(obj.tt));
            else
                tt1 = time_transfer(obj.tt,3);
                fprintf(' from %d-%d-%d to %d-%d-%d\n',tt1(1,:),tt1(end,:));
            end
            fprintf('\n');
        end
        
        function displonlat(obj, Ndisp)
            % Ndisp: number of lines to disp
            
            if nargin == 1
                Ndisp = 3;
            end
            
            ran = obj.ran;
            
            % lon {
            fprintf('lon: 1 -> %d; dlon = %.3f\n',obj.Nlon,ran(2));
            for ii = 1:Ndisp
                one_line(obj.lon, ii, Ndisp)
            end
            
            fprintf('      ...\n')
            one_line(obj.lon, obj.Nlat, Ndisp) % the last row: Nlat
            % }
            
            % lat {
            fprintf('\nlat: 1 -> %d; dlat = %.3f\n', obj.Nlon, ran(5));
            for ii = 1:Ndisp
                one_line(obj.lat, ii, Ndisp)
            end
            
            fprintf('      ...\n')
            one_line(obj.lat, obj.Nlat, Ndisp)
            % }
        end
        
        function rval = Nlon(obj)
            rval = size(obj.lon,2);
        end
        
        function rval = Nlat(obj)
            rval = size(obj.lon,1);
        end
        
        function rval = Nmonth(obj)
            rval = size(obj.rg,3);
        end
        
        function vran = ran(obj, iprint)
            if nargin == 1
                iprint = 0;
            end
            lon0 = min(obj.lon(:));
            lon1 = max(obj.lon(:));
            lat0 = min(obj.lat(:));
            lat1 = max(obj.lat(:));
            dlon = (lon1-lon0)/(obj.Nlon-1);
            dlat = (lat1-lat0)/(obj.Nlat-1);
            vran = [lon0,dlon,lon1,lat0,dlat,lat1];
            
            if iprint == 1
                fprintf('%s.ran: (%g, %g, %g, %g, %g, %g)\n',inputname(1), vran);
            end
        end
        
        function vran = ran4(obj, iprint)
            if nargin == 1
                iprint = 0;
            end
            lon0 = min(obj.lon(:));
            lon1 = max(obj.lon(:));
            lat0 = min(obj.lat(:));
            lat1 = max(obj.lat(:));
            vran = [lon0,lon1,lat0,lat1];
            
            if iprint == 1
                fprintf('%s.ran: (%g, %g, %g, %g)\n',inputname(1), vran);
            end
        end
        
        function tran=trange(obj)
            tt = obj.tt;
            if isnan(tt(1))
                tran = NaN;
                return
            end
            t0 = time_transfer(tt(1),3);
            t1 = time_transfer(tt(end),3);
            tran = [t0(1:2), t1(1:2)];
        end
        
        function rval = vlat(obj)
            rval = obj.lat(:,1);
        end
        
        function rval = vlon(obj)
            rval = obj.lon(1,:)';
        end
        
        %% check
        
        function iYN=compatible(LLZ1, LLZ2)
            
            if all(LLZ1.ran == LLZ2.ran)
                iYN = 1;
            else
                iYN = 0;
            end
        end
        
        %% calculation
        function obj3 = minus(obj1,obj2)
            if compatible(obj1,obj2) ~= 1
                if obj1.Nlon*obj1.Nlat > obj2.Nlon*obj2.Nlat
                    ichoice = 2;
                    obj1 = LLZ_resample(obj1,'stype','spatial','lonlat',obj2);
                else
                    ichoice = 1;
                    obj2 = LLZ_resample(obj2,'stype','spatial','lonlat',obj1);
                end
                schoice = {'1st','2nd'};
                fprintf('In cLLZ-minus, size of two LLZs: (%d,%d), (%d,%d), converted to the %s one\n',...
                    obj1.Nlon,obj1.Nlat,obj2.Nlon,obj2.Nlat, schoice{ichoice});
            end
            
            if obj2.Nmonth ~= 1
                if obj1.Nmonth == obj2.Nmonth
                    rtmp = obj1.rg - obj2.rg;
                else
                    error('The sizes of two inputs are %d and %d',obj1.Nmonth,obj2.Nmonth);
                end
            else
                for ii = 1:obj1.Nmonth
                    rtmp(:,:,ii) = obj1.rg(:,:,ii) - obj2.rg;
                end
            end
            obj3 = cLLZ(obj1.lon,obj1.lat,rtmp,obj1.tt);
            obj3 = obj3.setinfo(obj1);
            
        end
        
        function obj3 = plus(obj1,obj2)
            if compatible(obj1,obj2) ~= 1
                if obj1.Nlon*obj1.Nlat > obj2.Nlon*obj2.Nlat
                    ichoice = 2;
                    obj1 = LLZ_resample(obj1,'stype','spatial','lonlat',obj2);
                else
                    ichoice = 1;
                    obj2 = LLZ_resample(obj2,'stype','spatial','lonlat',obj1);
                end
                schoice = {'1st','2nd'};
                fprintf('In cLLZ-plus, size of two LLZs: (%d,%d), (%d,%d), converted to the %s one\n',...
                    obj1.Nlon,obj1.Nlat,obj2.Nlon,obj2.Nlat, schoice{ichoice});
            end
            
            if obj2.Nmonth ~= 1
                if obj1.Nmonth == obj2.Nmonth
                    rtmp = obj1.rg + obj2.rg;
                else
                    error('The sizes of two inputs are %d and %d',obj1.Nmonth,obj2.Nmonth);
                end
            else
                for ii = 1:obj1.Nmonth
                    rtmp(:,:,ii) = obj1.rg(:,:,ii) + obj2.rg;
                end
            end
            
            obj3 = cLLZ(obj1.lon,obj1.lat,rtmp,obj1.tt);
            obj3 = obj3.setinfo(obj1);
            
        end
        
        function LLZ=mean(obj,o_trange)
            % {{LLZ; mean}}
            % calculate the mean rg between o_trange
            
            if nargin == 2
                obj = obj.extract(o_trange);
            end
            
            rg = mean(obj.rg,3);
            
            LLZ = cLLZ(obj,rg,mean([obj.tt]));
        end
        
        %% reshape
        function [LLZ,ind] = cut(LLZ_G,ran)
            % ran: (lon(1),lon(2),lat(1),lat(2))
            % ran: (lon0, dlon, lon1, lat0, dlat, lat1)
            
            if numel(ran) == 6
                ran = ran([1,3,4,6]);
            end
            mlonG = LLZ_G.lon;
            mlatG = LLZ_G.lat;
            ind = mlonG>=ran(1) & mlonG<=ran(2) & mlatG>=ran(3) & mlatG<=ran(4);
            Nlon = max(sum(ind,2));
            Nlat = max(sum(ind,1));
            Nrg = size(LLZ_G.rg,3);
            
            rg(Nlat,Nlon,Nrg) = 0;
            for irg = 1:Nrg
                rgtemp = LLZ_G.rg(:,:,irg);
                rg(:,:,irg) = reshape(rgtemp(ind),Nlat,Nlon);
            end
            lon = reshape(mlonG(ind),Nlat,Nlon);
            lat = reshape(mlatG(ind),Nlat,Nlon);
            
            tt = LLZ_G.tt;
            
            LLZ= cLLZ(lon,lat,rg,tt);
            LLZ = LLZ.setinfo(LLZ_G);
            LLZ = LLZ.setinfo('history',sprintf('cut: %s',num2str(ran)));
        end
        
        function LLZ=extract(LLZ, ilist_or_trange)
            % {{LLZ; extract}}
            % trange : [y1,m1,y2,m2];
            % ilist: [1:3,6]
            
            if numel(ilist_or_trange) == 4 && min(ilist_or_trange([1,3])) > 1000 && max(ilist_or_trange([2,4])) <=12
                trange = ilist_or_trange;
                msg = sprintf('extract by trange: %s',num2str(trange));
            elseif numel(ilist_or_trange) == 2 && min(ilist_or_trange(1)) > 1000 && max(ilist_or_trange(2)) <=12
                trange = ilist_or_trange;
                trange = [trange,trange];
                msg = sprintf('extract by trange: %s',num2str(trange));
            else
                ilist = ilist_or_trange;
                if numel(ilist) > 10
                    msg = sprintf('extract by ilist: %s, ..., %s',num2str(ilist(1:5)),num2str(ilist(end-4:end)) );
                else
                    msg = sprintf('extract by ilist: %s',num2str(ilist));
                end
            end
            
            if exist('trange','var') == 1 % trange is set
                tt = LLZ.tt;
                t1 = trange(1) + (trange(2)-1)/12;
                t2 = trange(3) + trange(4)/12;
                loc1 = sum(tt<=t1)+1;
                loc2 = sum(tt<t2);
                fprintf('location range: %d, %d\n',loc1,loc2);
                ilist = loc1:loc2;
            end
            
            LLZ.rg = LLZ.rg(:,:,ilist);
            LLZ.tt = LLZ.tt(ilist);
            LLZ = LLZ.setinfo('history',msg);
        end
        
        function LLZ=toStandard(LLZ)
            % {{standard}} {{LLZ}}
            %
            
            % --check{
            if LLZ.Nlon == 1 && LLZ.Nlat == 1
                LLZ.info.type = 'single';
            elseif LLZ.Nlon == 1 || LLZ.Nlat == 1
                LLZ.info.type = 'vector';
            else
                dlon = LLZ.lon(2:end,2:end) - LLZ.lon(1:end-1,1:end-1);
                dlat = LLZ.lat(2:end,2:end) - LLZ.lat(1:end-1,1:end-1);
                if max( max(dlon(:))-min(dlon(:)), max(dlat(:))-min(dlat(:)) ) > 1D-10
                    LLZ.info.type = 'irr'; % irregular distributed
                else
                    LLZ.info.type = 'normal';
                end
            end
            % --check}
            
            if ~strcmp(LLZ.info.type,'normal')
                fprintf('LLZ type is %s, no standard check\n',LLZ.info.type);
                return
            end
            
            % 1. rotate
            if LLZ.lon(2,1) ~= LLZ.lon(1,1) % need to rotate
                rgs = zeros(LLZ.Nlon,LLZ.Nlat,LLZ.Nmonth);
                LLZ.lon = LLZ.lon';
                LLZ.lat = LLZ.lat';
                
                for ii = 1:LLZ.Nmonth
                    rgs(:,:,ii) = LLZ.rg(:,:,ii)';
                end
                LLZ.rg = rgs;
                fprintf('LLZ is rotated to standard\n');
            end
            
            % 2. flip longitude
            if LLZ.lon(1,2) < LLZ.lon(1,1) % lon in descending order
                LLZ.lon = fliplr(LLZ.lon);
                for ii = 1:LLZ.Nmonth
                    LLZ.rg(:,:,ii) = fliplr(LLZ.rg(:,:,ii));
                end
                %     fprintf('LLZ is flipped in longitude to standard\n');
            end
            
            % 3. flip latitude
            if LLZ.lat(2,1) < LLZ.lat(1,1) % lat in descending order
                LLZ.lat = flipud(LLZ.lat);
                for ii = 1:LLZ.Nmonth
                    LLZ.rg(:,:,ii) = flipud(LLZ.rg(:,:,ii));
                end
                %     fprintf('LLZ is flipped in latitude to standard\n');
            end
            % tt
            tt = LLZ.tt;
            if numel(tt) ~= LLZ.Nmonth
                if all(isnan(tt))
                    LLZ.tt = ones(LLZ.Nmonth,1) *nan;
                else
                    error('tt are not compatible in LLZ');
                end
            end
        end
        
        % output
        function save(obj,opfile)
            % {{save; object; struct}}
            
            LLZ = classToStruct(obj);
            
            save(opfile,'LLZ');
            
        end
        
        function save_nc(obj,sfile)
            %{
write cLLZ in the format of netcdf
% {{netcdf; LLZ; write}}
            %}
            
            if nargin == 0
                fun_ex1;
                return;
            end
            % ----- end of head -----
            Nmonth = obj.Nmonth;
            
            Nlat = obj.Nlat;
            Nlon = obj.Nlon;
            
            nccreate(sfile,'lon',...
                'Dimensions', {'Nlat',Nlat,'Nlon',Nlon},...
                'FillValue','disable');
            ncwrite(sfile,'lon',obj.lon);
            ncwriteatt(sfile,'lon','description','Longitude (degree)');
            
            nccreate(sfile,'lat',...
                'Dimensions', {'Nlat',Nlat,'Nlon',Nlon},...
                'FillValue','disable');
            ncwrite(sfile,'lat',obj.lat);
            ncwriteatt(sfile,'lat','description','Latitude (degree)');
            
            nccreate(sfile,'rg',...
                'Dimensions', {'Nlat',Nlat,'Nlon',Nlon,'Ntime',Nmonth},...
                'FillValue','disable');
            ncwrite(sfile,'rg',obj.rg);
            ncwriteatt(sfile,'rg','description','Gridded observations');
            
            nccreate(sfile,'time',...
                'Dimensions', {'Ntime',Nmonth},...
                'FillValue','disable');
            ncwrite(sfile,'time',obj.tt);
            ncwriteatt(sfile,'time','description','Time in decimal years');
            
            ncwriteatt(sfile,'/','creation_date',datestr(now));
            ncwriteatt(sfile,'/','Name',obj.info.name);
            ncwriteatt(sfile,'/','Unit',obj.info.unit);
        end
    end
    
    methods (Static)
        obj = load(ipfile);
        function obj = ex1
            [mlon,mlat] = meshgrid(1:5,11:16);
            rg(:,:,1) = mlon+mlat;
            rg(:,:,2) = mlon-mlat;
            rg(:,:,3) = -mlon+mlat;
            tt = generate_tt([2003,1],[2003,3]);
            obj = cLLZ(mlon,mlat,rg,tt);
            obj = obj.setinfo('unit','cm','name','ex1','history','create in the ex1 of cLLZ');
        end
    end
    
end

function one_line(mat, iline, Ndisp)
fprintf('  [%3d] ',iline)
fprintf('%.3f, ',mat(iline,1:Ndisp));
fprintf('..., ')
fprintf('%.3f, ',mat(iline,end-Ndisp+1:end));
fprintf('\n');
end
