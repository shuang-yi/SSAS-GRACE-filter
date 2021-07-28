classdef cSER
    % {{series; class}}
    properties
        tt
        y
        dy
        info % name, unit, history
    end
    
    methods
        % create {
        % 1. ser (.tt, .y, .dy, .name, .unit)
        % 2. tt,y,dy, name, unit
        % 3. cSER0, tt, y, dy
        function obj = cSER(varargin)
            
            
            if nargin == 3
                % 2. tt, y, dy
                [tt,y,dy] = deal(varargin{1:3});
                [n1,m1] = size(tt);
                [n2,m2] = size(y);
                if isnan(dy)
                    dy = nan(n2,m2);
                end
                [n3,m3] = size(dy);
                if ~all([n2,m2] == [n3,m3])
                    error('size of y (%d,%d) and dy (%d,%d) are not compatible',n2,m2,n3,m3);
                end
                if n1 == 1
                    if m1 ~= m2
                        error('size of tt (%d,%d) and y (%d,%d) are not compatible',n1,m1,n2,m2)
                    end
                    Nmonth = n2;
                    tt = tt';
                    y = y';
                    dy = dy';
                elseif m1 == 1
                    if n1 ~= n2
                        error('size of tt (%d,%d) and y (%d,%d) are not compatible',n1,m1,n2,m2)
                    end
                    Nmonth = m2;
                else
                    error('size of tt (%d,%d) is not supported, either dimension of it should be 1',n1,m2);
                end
                
                if Nmonth == 1
                    obj.tt = tt;
                    obj.y = y;
                    obj.dy = dy;
                    obj.info = datainfo;
                else
                    for ii = 1:Nmonth
                        obj(ii) = cSER(tt,y(:,ii),dy(:,ii));
                    end
                end
                
            elseif nargin == 1 && isstruct(varargin{1}) && isfield(varargin{1},'y') % 1. ser.tt,ser.y,ser.dy,ser.name
                ser0 = varargin{1};
                for ii = 1:numel(ser0)
                    if ~isfield(ser0,'dy')
                        ser0(ii).dy = nan(size(ser0(ii).y));
                    end
                    if ~isfield(ser0,'unit')
                        ser0(ii).unit = PAR.unit;
                    end
                    
                    obj(ii) = cSER(ser0(ii).tt,ser0(ii).y,ser0(ii).dy);
                end
            else
                error('unsupported input');
            end
            
        end
        % create }
        
        %% 1. query
        function Nmon = Nmonth(obj)
            for ii = 1:numel(obj)
                Nmon(ii) = numel(obj(ii).tt);
            end
        end
        
        function tran = trange(obj)
            for ii = 1:numel(obj)
                tran(ii,:) = trange_one(obj(ii));
            end
        end
        
        function obj = setinfo(obj,varargin)
            % SER.setinfo('name','xx','unit','xx','type','xx','history',{'xx','yy'});
            % SER.setinfo(struct('name','xx'));
            % SER.setinfo(SER2);
            %
            % multiple history should be stored in cell
            
            % {{cSER; set; info}}
            
            if isobject(varargin{1}) % input is an object
                a = varargin{1};
                for ii = 1:numel(obj)
                    rtmp = obj(ii).info;
                    
                    if numel(a) == 1  % asign one to multiple
                        obj(ii).info = rtmp.set(a.info);
                    elseif numel(a) == numel(obj) % asign multiple to multiple
                        obj(ii).info = rtmp.set(a(ii).info);
                    else
                        error('The number of infos is not capatible\n');
                    end
                    
                end
            else % input is {'name','xx'}
                a = struct(varargin{:});
                for ii = 1:numel(obj)
                    rtmp = obj(ii).info;
                    if numel(a) == 1  % asign one to multiple
                        obj(ii).info = rtmp.set(a);
                    elseif numel(a) == numel(obj) % asign multiple to multiple
                        obj(ii).info = rtmp.set(a(ii));
                    else
                        error('The number of infos is not capatible\n');
                    end
                    
                end
            end
            
            
        end
        
        function print(obj, iser)
            % {{print}}
            
            fprintf('===== info =====\n');
            
            if numel(obj)>1
                if nargin == 1
                    iser = 1;
                    fprintf('more the one input, only the %d is used\n',iser);
                end
                obj = obj(iser);
            end
            
            
            obj.info.print;
            
            
            fprintf(' - range of y \n');
            rg = obj.y;
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
            
            fprintf(' - total samples: %d\n',numel(obj.tt));
            fprintf('\n');
        end
        
        %% 2. reformat
        function varargout = reformat(obj,pattern)
            % {{reformat}}
            % [mat,sname] = reformat('matrix')
            
            if strcmp(pattern,'matrix')
                a = obj.Nmonth;
                if numel(unique(a))~= 1
                    disp(a);
                    error('cSER.reformat: the size of time series are not compatible\n');
                end
                
                mat(:,1) = obj(1).tt;
                ss{1} = 'time';
                for ii = 1:numel(obj)
                    mat(:,ii+1) = obj(ii).y;
                    mat_delta(:,ii) = obj(ii).dy;
                    ss{ii+1} = obj(ii).info.name;
                end
                varargout = {mat,mat_delta,ss};
                
            else
                error('unrecognized patter %s',pattern)
            end
        end
        
        %% setval
        function SER1 = setval(SER1,s_tran_seq,tran_seq,ser)
            % {{SER; setvalue}}
            % setval(SER1, 'tran', [2003,1,2003,12], 0)
            % setval(SER1, 'tran', [2003,1,2003,12], 1:12)
            % setval(SER1, 'seq', 3:15, 0)
            for ii = 1:numel(SER1)
                SER1(ii) = setone(SER1(ii),s_tran_seq,tran_seq,ser);
            end
        end
        
        %% resize
        function SER = extract(SER0,tran)
            % trange : 1. a single time ( [year,month] )
            %          2. time range ( [year1,month1,yaer2,month2] )
            for ii = 1:numel(SER0)
                SER(ii) = extract_one(SER0(ii),tran);
            end
        end
        
        function SER = extract_one(SER0,tran)
            tt = SER0.tt;
            
            if length(tran) == 2
                [minval,minloc] = min( abs( tt-generate_tt(tran) ) );
                tt = tt(minloc); y = SER0.y(minloc); dy = SER0.dy(minloc);
                SER = cSER(SER0,tt,y,dy);
                if minval ~= 0
                    fprintf('minloc= %d; min distances are %7.4f days\n',minloc,minval*365);
                end
                ss = sprintf('[%d,%d]',tran);
            elseif length(tran) == 4
                t1 = generate_tt(tran(1:2))-0.5/12; % first day of the month
                t2 = generate_tt(tran(3:4))+0.5/12; % end day of the month
                
                loc_t1 = sum(tt<t1)+1;
                loc_t2 = sum(tt<=t2);
                
                fprintf('loc_t1 = %d, loc_t2 = %d\n',loc_t1,loc_t2);
                
                ran = loc_t1:loc_t2;
                tt = tt(ran);
                y = SER0.y(ran);
                dy = SER0.dy(ran);
                SER = cSER(tt,y,dy);
                SER = SER.setinfo(SER0);
                ss = sprintf('[%d,%d, %d,%d]',tran);
            else
                error('length(trange) = %d',length(tran));
            end
            
            SER = SER.setinfo('history',sprintf('extract by %s',ss));
        end
        
        %% plot
        function plot(obj,varargin)
            PAR0 = struct('stype','o-','simple',0);
            PAR = var_initial(PAR0,varargin);
            
            if PAR.simple == 1 % only plot series
                for ii = 1:numel(obj)
                    plot(obj(ii).tt,obj(ii).y,PAR.stype);
                    hold on;
                end
                return;
            end
            
            for ii = 1:numel(obj)
                plot(obj(ii).tt,obj(ii).y,PAR.stype);
                ss{ii} = sprintf('%s',strstr(obj(ii).info.name,'title'));
                hold on;
                
                tran(ii,1:2) = obj(ii).tt([1,end]);
            end
            legend(ss,'location','best');
            if ~isempty(obj(1).info.unit)
                ylabel(obj(1).info.unit);
            end
            grid on;
            
            if obj(1).tt(1)>1800
                axis_year(floor(min(tran(:,1))),ceil(max(tran(:,2))));
            end
            
            hold off;
        end
        
        function eplot(obj,varargin)
            PAR0 = struct('iplot',1,'stype','o-');
            PAR = var_initial(PAR0,varargin);
            scolor = mycolors(numel(obj));
            
            for ii = 1:numel(obj)
                iplot = PAR.iplot;
                % single point:
                if numel(obj(ii).tt) == 1; iplot = 2;  end
                
                if iplot == 2 % errorbar
                    hh(ii) = errorbar(obj(ii).tt,obj(ii).y,obj(ii).dy,PAR.stype);
                elseif iplot == 1 % fill plot
                    ser = obj(ii);
                    if ~all(isnan(ser.dy))
                        h(1) = fill_two_line(ser.tt,ser.y+ser.dy,ser.tt,ser.y-ser.dy);
                        set(h(1),'facecolor',lighter(scolor(ii,:),0.5));
                        hold on;
                    end
                    hh(ii) = plot(ser.tt,ser.y,'color',scolor(ii,:),'linew',2);
                end
                ss{ii} = sprintf('%s (%s)',strstr(obj(ii).info.name,'title'),obj(ii).info.unit);
                hold on;
                tran(ii,1:2) = obj(ii).tt([1,end]);
            end
            legend(hh,ss,'location','best');
            
            if obj(1).tt(1)>1800
                axis_year(floor(min(tran(:,1))),ceil(max(tran(:,2))));
            end
            grid on;
            hold off;
        end
        
        %% calculation
        function obj3 = minus(obj1,obj2)
            % {{SER; minus}}
            
            if isnumeric(obj2) % numeric
                
                for ii = 1:numel(obj1)
                    if numel(obj2) == 1 % case 1: single value
                        a(ii) = cSER(obj1(ii),NaN,obj2,0);
                    else  % case 2: time series
                        if numel(obj1(ii).tt) ~= numel(obj2)
                            error('the numeric sequence should either be single value, or have the same size as the SER')
                        end
                        a(ii) = cSER(obj1(ii).tt,obj2(:),obj2(:)*0,'inherit',obj1(ii));
                    end
                end
                obj2 = a;
                
            elseif numel(obj2) == 1 % object
                obj2 = obj2(ones(numel(obj1)));
            end
            for ii = 1:numel(obj1)
                obj3(ii) = minus_one(obj1(ii),obj2(ii));
            end
        end
        
        %%
        function obj3 = minus_one(obj1,obj2)
            
            if (obj1.Nmonth ~= obj2.Nmonth) && obj2.Nmonth ~= 1
                error('The Mmonth of two inputs are %d and %d',obj1.Nmonth,obj2.Nmonth);
            end
            
            y = obj1.y - obj2.y;
            dy = sqrt(obj1.dy.^2 + obj2.dy.^2);
            % name = sprintf('%s minus %s',obj1.info.name,obj2.info.name);
            obj3 = cSER(obj1.tt,y,dy);
            obj3 = obj3.setinfo(obj1);
            % obj3 = obj3.setinfo('name',name,'history',sprintf('minus: %s',name));
        end
        
        function obj3 = plus(obj1,obj2)
            % {{SER; plus}}
            
            if isnumeric(obj2) % numeric
                
                for ii = 1:numel(obj1)
                    if numel(obj2) == 1 % case 1: single value
                        a(ii) = cSER(obj1(ii),NaN,obj2,0);
                    else  % case 2: time series
                        if numel(obj1(ii).tt) ~= numel(obj2)
                            error('the numeric sequence should either be single value, or have the same size as the SER')
                        end
                        a(ii) = cSER(obj1(ii).tt,obj2(:),obj2(:)*0,'inherit',obj1(ii));
                    end
                end
                obj2 = a;
                
            elseif numel(obj2) == 1 % object
                obj2 = obj2(ones(numel(obj1)));
            end
            for ii = 1:numel(obj1)
                obj3(ii) = plus_one(obj1(ii),obj2(ii));
            end
        end
        
        function obj3 = mtimes(obj1,a)
            % {{SER; plus}}
            
            if ~isnumeric(a) || numel(a)~= 1 % single numeric
                
                error('input should be single numeric');
                
            end
            for ii = 1:numel(obj1)
                obj3(ii) = mtimes_one(obj1(ii),a);
            end
        end
        
        function SER_mean=mean(obj,o_trange)
            % {{SH; mean}}
            if nargin == 1
                o_trange = [];
            end
            for ii = 1:numel(obj)
                SER_mean(ii) = mean_one(obj(ii),o_trange);
            end
        end
        
        function SER=demean(obj, o_trange)
            % {{mean; demean}}
            if nargin == 2
                SERm = obj.mean(o_trange);
                ss = sprintf('remove mean value during [%s]',num2str(o_trange));
            else
                SERm = obj.mean;
                ss = sprintf('remove mean value of the whole period');
            end
            
            SER = obj - SERm;
            SER = SER.setinfo('history',ss);
        end
        
    end
    methods (Static)
        function obj = ex1
            tt = generate_tt([2000,1],[2010,12]);
            y = [sin(2*pi*tt)+(tt-2000)*0.1,cos(2*pi*tt)-(tt-2000)*0.1];
            % y = [(1:132)',(2:133)'];
            dy = ones(size(y));
            name = {'sin','cos'};
            
            obj = cSER(tt,y,dy);
            obj = obj.setinfo('name',name,'unit','cm');
            obj.plot;
                       
        end
        
    end
end

function tran = trange_one(obj)
tt = obj.tt;
if isnan(tt(1))
    tran = NaN;
    return
end
t0 = time_transfer(tt(1),3);
t1 = time_transfer(tt(end),3);
tran = [t0(1:2), t1(1:2)];
end

function SER1 = setone(SER1,s_tran_seq,tran_seq,ser)
if strcmp(s_tran_seq,'tran')
    tran = tran_seq;
    tt = SER1.tt;
    if length(tran) == 2
        [~,minloc] = min( abs( tt-generate_tt(tran) ) );
        ran = minloc;
        
        ss = sprintf('set sample in [%d,%d] to be %s',tran, num2str(ser));
    elseif length(tran) == 4
        t1 = generate_tt(tran(1:2))-0.5/12; % first day of the month
        t2 = generate_tt(tran(3:4))+0.5/12; % end day of the month
        
        loc_t1 = sum(tt<t1)+1;
        loc_t2 = sum(tt<=t2);
        
        ran = loc_t1:loc_t2;
        
        ss = sprintf('ser samples within [%d,%d, %d,%d] to be %s',tran, num2str(ser));
        
    else
        error('length(trange) = %d',length(tran));
    end
elseif strcmp(s_tran_seq,'seq')
    ran = tran_seq;
    ss = sprintf('ser samples in %s to be %s',num2str(ran), num2str(ser));
end

if numel(ser) == 1
    ser = ones(size(ran))*ser;
end
if numel(ser) ~= numel(ran)
    error('setval: the size of input ser is %s, while the size of selected sampels is %s\n',num2str(size(ser)),num2str(size(ran)));
end
SER1.y(ran) = ser;

SER1 = SER1.setinfo('history',ss);

end

function obj3 = plus_one(obj1,obj2)

if (obj1.Nmonth ~= obj2.Nmonth) && obj2.Nmonth ~= 1
    error('The Mmonth of two inputs are %d and %d',obj1.Nmonth,obj2.Nmonth);
end

y = obj1.y + obj2.y;
dy = sqrt(obj1.dy.^2 + obj2.dy.^2);
% name = sprintf('%s plus %s',obj1.info.name,obj2.info.name);
obj3 = cSER(obj1.tt,y,dy);
obj3 = obj3.setinfo(obj1);
% obj3 = cSER(obj1.tt,y,dy,'inherit',obj1);
% obj3 = obj3.setinfo('name',name,'history',sprintf('plus: %s',name));
end

function obj3 = mtimes_one(obj1,a)

y = obj1.y *a;
dy = obj1.dy *a;
name = sprintf('%s times %g',obj1.info.name,a);
obj3 = cSER(obj1.tt,y,dy);
obj3 = obj3.setinfo(obj1);

obj3 = obj3.setinfo('history',sprintf('mtimes: %s',name));
end

function SER_mean = mean_one(obj,trange)
if ~isempty(trange)
    SER = obj.extract(trange);
else
    SER = obj;
    trange = [inf,inf];
end

tt = mean(SER.tt);
y = mean(SER.y);
dy = sqrt(sum(SER.dy.^2))/SER.Nmonth;

SER_mean = cSER(tt,y,dy);
SER_mean = SER_mean.setinfo(obj);
SER_mean = SER_mean.setinfo('history',sprintf('mean value over %s', mat2str(trange)));
end