classdef cSH
    % 1. create
    %  * SH
    %  * rC, rS, tt
    %  * sc, tt, 'sc'
    %  * cs, tt, 'cs'
    %  * lmcosi, tt, 'lmcosi'           % single month
    %  * vCS, LIST, tt, 'vCS'
    %  * Nmax, tt
    %
    % 2. query or return
    %  Nmonth, Nmax, trange
    %  optt, nm(n,m),degree(n)
    % 3. reformat
    %  [rC,rS,tt] = reformat('matrix')
    %  [vCS,LIST] = reformat('vCS',o_in0)
    %  [lmcosi] = reformat('lmcosi')
    % 4. setval
    %  setval(nlist,mlist,ser)
    % 5. resize
    %  resize(Nmax), extract(trange), remove_by_ym(tlist), remove_by_seq(seq_list)
    % 6. output
    %  optt
    %  nm(n,m)
    %  degree(n);
    %  output2file(opfile)
    % 7. plot
    %  plot_t('o-'); preview; kaula('o-');
    %  plot; plot_abs; plot_nm(n,m); plot_n; plot_m
    % 8. calculation
    %  minus, plus, times(a), mean(trange), demean(trange)
    %  equal(SH1,SH2)
    % 9. convert
    %  unit_convert(ifrom,ito)
    
    % {{SH; class}}
    properties
        C
        S
        tt
        info % unit, name, history
    end
    
    methods
        %% create {
        
        function obj = cSH(varargin)
            
            if nargin == 0
                msg = {'==1. how to create a cSH variable:'
                    'cSH(SH)'
                    'cSH(rC, rS, tt)'
                    'cSH(scs, tt, ''sc'') '
                    'cSH(css, tt, ''cs'') '
                    "cSH(lmcosi, tt, 'lmcosi') % single month"
                    "cSH(vCS, LIST, tt, 'vCS')"
                    "cSH(Nmax, tt)"};
                for ii = 1:numel(msg)
                    fprintf('%s\n',msg{ii});
                end
                
                fprintf('\n==2. all methods\n');
                methods(cSH(0,0));
            elseif nargin == 3
                if isnumeric(varargin{3}) % 1.1. rC, rS ,tt
                    [rC,rS,tt] = deal(varargin{1:3});
                    Nmon = size(rC,3);
                    if Nmon ~= numel(tt)
                        if isnan(tt)
                            tt = NaN*ones(Nmon,1);
                        else
                            error('input sizes are not consistent: C(%d), S(%d), tt(%d)',size(rC,3),size(rS,3),numel(tt));
                        end
                    end
                    if Nmon == 1
                        obj.C = rC;
                        obj.S = rS;
                        obj.tt = tt;
                        obj.info = datainfo;
                    else
                        for ii = 1:Nmon
                            obj(ii) = cSH(rC(:,:,ii),rS(:,:,ii),tt(ii));
                        end
                    end
                elseif strcmp(varargin{3},'sc')  % sc, tt, 'sc'
                    [scs,tt] = deal(varargin{1:2});
                    Nmax = size(scs,1)-1;
                    if size(scs,2) ~= 2*Nmax+1
                        error('the size of the input cannot be recognized to be SC: %d * %d\n',size(scs,1),size(scs,2));
                    end
                    for ii = 1:size(scs,3)
                        sc = scs(:,:,ii);
                        rS = zeros(Nmax+1);
                        rS(:,2:Nmax+1) = fliplr(sc(:,1:Nmax));
                        rC = sc(:,Nmax+1:2*Nmax+1);
                        
                        obj(ii) = cSH(rC,rS,tt(ii));
                    end
                elseif strcmp(varargin{3},'cs')  % sc, tt, 'cs'
                    [css,tt] = deal(varargin{1:2});
                    Nmax = size(css,1)-1;
                    if size(css,2) ~= Nmax+1
                        error('the size of the input cannot be recognized to be CS: %d * %d\n',size(css,1),size(css,2));
                    end
                    for ii = 1:size(css,3)
                        cs = css(:,:,ii);
                        rC = tril(cs);
                        rS = triu(cs,1);
                        rS = rS([end,1:end-1],:)';
                        
                        obj(ii) = cSH(rC,rS,tt(ii));
                    end
                elseif strcmp(varargin{3},'lmcosi')  % 1.3. lmcosi, tt, 'lmcosi'
                    lmcosi = varargin{1};
                    Nmax = max(lmcosi(:,1));
                    rC = zeros(Nmax+1);
                    rS = rC;
                    tt = varargin{2};
                    
                    Nc = size(lmcosi,1);
                    for ik = 1:Nc
                        in = lmcosi(ik,1);
                        im = lmcosi(ik,2);
                        rC(in+1,im+1) = lmcosi(ik,3);
                        rS(in+1,im+1) = lmcosi(ik,4);
                    end
                    
                    obj = cSH(rC,rS,tt);
                end
                
            elseif nargin == 1 && ((isstruct(varargin{1}) && isfield(varargin{1},'C')) || isa(obj,'cSH')) % 1. SH ['C','S','tt']
                SH = varargin{1};
                if numel(SH) == 1
                    obj = cSH(SH.C,SH.S,SH.tt);
                else
                    for ii = 1:numel(SH)
                        obj(ii) = cSH(SH(ii));
                        try % inherit the info, but it may fail
                            obj(ii) = obj(ii).setinfo(SH(ii).info);
                        catch
                            warning('info is not read in correctly');
                        end
                    end
                end
                
            elseif nargin == 4 && strcmp(varargin{4},'vCS')  % 3. vCS, LIST, tt, 'vCS'
                [vCS,LIST,tt] = deal(varargin{1:3});
                Nmax = max(LIST(:,2));
                
                [Ncs,Nmonth] = size(vCS);
                
                Ncs2 = size(LIST,1);
                if Ncs ~= Ncs2
                    disp(['length(vCS):',num2str(Ncs),';  ','size(LIST,1):',num2str(Ncs2)]);
                    error('vCS_to2D: LIST is not compatable with vCS');
                end
                
                rCs = NaN(Nmax+1,Nmax+1,Nmonth);
                rSs = rCs;
                
                for ii = 1:Ncs
                    in = LIST(ii,1);
                    im = LIST(ii,2);
                    if im >= 0 % C
                        if ~isnan(rCs(in+1,im+1,:))
                            mydisp('warning: duplicate inputs, should m < 0 in the LIST?\n');
                        end
                        rCs(in+1,im+1,:) = vCS(ii,:);
                    else  %S
                        if ~isnan(rSs(in+1,-im+1,:))
                            mydisp('warning: duplicate inputs\n');
                        end
                        rSs(in+1,-im+1,:) = vCS(ii,:);
                    end
                    
                end
                rCs(isnan(rCs)) = 0;
                rSs(isnan(rSs)) = 0;
                obj = cSH(rCs,rSs,tt);
                
            elseif nargin == 2 && numel(varargin{1}) == 1 % Nmax, tt
                Nmax = varargin{1};
                tt = varargin{2};
                rC = zeros(Nmax+1,Nmax+1,numel(tt));
                rS = rC;
                obj = cSH(rC,rS,tt);
            else
                error('unsupported input');
            end
            
        end
        %% create }
        
        % %query
        function obj = setinfo(obj,varargin)
            % SH.setinfo('name','xx','unit','xx','type','xx','history',{'xx','yy'});
            % SH.setinfo(struct('name','xx'));
            % SH.setinfo(SH2);
            %
            % multiple history should be stored in cell
            
            % {{cSH; set; info}}
            
            for ii = 1:numel(obj)
                rtmp = obj(ii).info;
                %                 if ~isa(rtmp,'datainfo') % in case the info type is the old version
                %                     rtmp = datainfo(rtmp);
                %                 end
                
                if isa(varargin{1},'cSH') % input is another cSH
                    
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
        
        function print(obj)
            % {{cSH; dispinfo; print}}
            
            
            fprintf('===== info =====\n');
            tt = obj.optt;
            if numel(obj)>1
                warning('more the one input, only the first is used\n');
                obj = obj(1);
            end
            
            obj.info.print;
            
            fprintf(' - Nmax: %d;   Nmonth: %d\n',obj.Nmax, numel(tt));
            
            fprintf(' - range of rC \n');
            
            [lmcosi] = obj.reformat('lmcosi');
            rg = lmcosi(:,3);
            
            fprintf('   (min) %g; (mean) %g; (median) %g; (max) %g;\n',...
                min(rg), mean(rg), median(rg), max(rg));
            
            fprintf(' - range of rS\n');
            rg = lmcosi(:,4);
            fprintf('   (min) %g; (mean) %g; (median) %g; (max) %g;\n',...
                min(rg), mean(rg), median(rg), max(rg));
            
            fprintf(' - trange:')
            if all(isnan(tt))
                fprintf(' %d * NaN\n',numel(tt));
            else
                tt1 = time_transfer(tt,3);
                fprintf(' from %d-%d-%d to %d-%d-%d\n',tt1(1,:),tt1(end,:));
            end
            fprintf('\n');
        end
        
        %% search
        function obj2=find(obj,varargin)
            % filter cSH by name/type/unit/history, contain mode
            
            % {{cSH; filter}}
            
            if mod(nargin,2) ~= 1
                mydisp('error: the number of input (%d) should be even\n',nargin-1);
            end
            
            for ii = 1:numel(obj)
                info(ii) = obj(ii).info;
            end
            
            ind = info.find(varargin{:});
            obj2 = obj(ind);
        end
        
        function Nmon = Nmonth(obj)
            Nmon = numel(obj);
        end
        
        function Nm = Nmax(obj)
            Nm = size(obj(1).C,1)-1;
        end
        
        function tran=trange(obj)
            % return the time range of the input
            % {{SH; tt; time; trange}}
            obj.tt = obj.optt;
            if isnan(tt(1))
                tran = NaN;
                return
            end
            t0 = time_transfer(tt(1),3);
            t1 = time_transfer(tt(end),3);
            tran = [t0(1:2), t1(1:2)];
        end
        
        %% reformat
        function varargout = reformat(obj,pattern,o_in0)
            % {{SH; reformat}}
            % [rC,rS,tt] = reformat('matrix')
            % [vCS,LIST] = reformat('vCS',o_in0)
            % [lmcosi] = reformat('lmcosi')
            if strcmp(pattern,'matrix')
                rC = zeros(obj.Nmax+1,obj.Nmax+1,obj.Nmonth);
                rS = rC;
                tt = zeros(obj.Nmonth,1);
                for ii = 1:numel(obj)
                    rC(:,:,ii) = obj(ii).C;
                    rS(:,:,ii) = obj(ii).S;
                    tt(ii) = obj(ii).tt;
                end
                varargout = {rC,rS,tt};
                
            elseif strcmp(pattern,'lmcosi')
                
                if obj.Nmonth > 1
                    warning('only support one month, consider ''vCS'' for multiple months');
                    obj = obj(1);
                end
                
                Nmax = obj.Nmax;
                
                Nc = (Nmax+1)*(Nmax+2)/2;
                lmcosi( Nc,4 ) =0;
                ik = 0;
                for in = 0:Nmax
                    for im = 0:in
                        ik = ik+1;
                        lmcosi(ik,:) = [in,im,obj.C(in+1,im+1),obj.S(in+1,im+1)];
                    end
                end
                varargout = {lmcosi};
            elseif strcmp(pattern,'vCS')
                if nargin == 2
                    o_in0 = 0;
                end
                [vCS,LIST]=SH_to1D(obj,o_in0);
                varargout = {vCS,LIST,obj.optt};
            elseif strcmp(pattern,'sc') % NS's format
                %{
    0      0      C(0,0) 0      0
    0      S(1,1) C(1,0) C(1,1) 0
    S(2,2) S(2,1) C(2,0) C(2,1) C(2,2)
                %}
                Nmax = obj.Nmax;
                Nmonth = obj.Nmonth;
                mat = zeros(Nmax+1,2*Nmax+1,Nmonth);
                for ii = 1:Nmonth
                    rC = obj(ii).C;
                    rS = obj(ii).S;
                    rS = rS(:,2:end);
                    mat(:,:,ii) = [fliplr(rS),rC];
                end
                varargout = {mat,obj.optt};
            elseif strcmp(pattern,'cs') % NS's format
                Nmax = obj.Nmax;
                Nmonth = obj.Nmonth;
                
                mat = zeros(Nmax+1,Nmax+1,Nmonth);
                for ii = 1:Nmonth
                    rC = obj(ii).C;
                    rS = obj(ii).S;
                    rS = rS';
                    rS = rS([2:end,1],:);
                    mat(:,:,ii) = tril(rC) + triu(rS,1);
                end
                varargout = {mat,obj.optt};
            else
                error('unrecognized patter %s',pattern)
            end
            
        end
        
        %% setval
        function SH2 = setval(SH1,nlist,mlist,ser)
            % {{SH; setvalue}}
            % SH2=sh_setval(SH1,1,-1, 2); % set all coefficients of (n,m)=(1,-1) to 2
            % SH2=sh_setval(SH1,'=2',[], 2); % set all coefficients of n=2 to 2
            % SH2=sh_setval(SH1,[],'>3 abs', 3); % set all coefficients of abs(m)>3 = 3
            % SH2 = sh_setval(SH1,1,1, [1,2,3]); % assign a series, only one element can be targeted in this case.
            %
            % See also: get_idx
            
            if nargin ~= 4
                error('four inputs are expected');
            end
            
            if ~ischar(nlist) && ~isempty(nlist) % change 1 to '=1'
                if isnumeric(nlist) && numel(nlist) == 1
                    nlist = sprintf('=%d',nlist);
                else
                    error('unsupported input\n');
                end
            end
            
            if ~ischar(mlist) && ~isempty(mlist) % change 1 to '=1'
                if isnumeric(mlist) && numel(mlist) == 1
                    mlist = sprintf('=%d',mlist);
                else
                    error('unsupported input\n');
                end
            end
            
            Nmax = SH1.Nmax;
            Nmonth = SH1.Nmonth;
            SH_idx = cSH.get_idx(Nmax,nlist,mlist);
            idx = SH_idx.reformat('vCS');
            
            [vCS,LIST,tt] = SH1.reformat('vCS');
            
            if numel(ser) > 1 && sum(idx) > 1
                mydisp('error: either the index targets only one coefficient, or the ser contain only one element\n')
            end
            
            if numel(ser) > 1
                if numel(ser) ~= Nmonth
                    mydisp('error: length of the ser is %d, while length of input SH is %d\n',numel(ser),Nmonth)
                end
            end
            
            vCS(idx==1,:) = ser;
            
            if numel(ser) <= 5
                s_ser = sprintf('size = %d, %s',numel(ser),mat2str(ser));
            else
                s_ser = sprintf('size = %d, %s, ...',numel(ser),mat2str(ser(1:5)));
            end
            
            if ~ischar(nlist)
                nlist = mat2str(nlist);
            end
            if ~ischar(mlist)
                mlist = mat2str(mlist);
            end
            
            history = sprintf('setval: n(%s),m(%s),series(%s)',...
                nlist,mlist,s_ser);
            
            SH2 = cSH(vCS,LIST,tt,'vCS');
            SH2 = SH2.setinfo(SH1);
            SH2 = SH2.setinfo('history',history);
        end
        
        %% resize
        function SH = resize(obj,Nmax2)
            % {{SH; resize; change}}
            if Nmax2 > obj.Nmax % expand the size
                aa = zeros(Nmax2+1);
                bb = aa;cc = aa;
                for ii = 1:obj.Nmonth
                    bb(1:obj.Nmax+1,1:obj.Nmax+1) = obj(ii).C;
                    cc(1:obj.Nmax+1,1:obj.Nmax+1) = obj(ii).S;
                    SH(ii) = cSH(bb,cc,obj(ii).tt);
                    SH(ii) = SH(ii).setinfo(obj(ii));
                end
                SH = SH.setinfo('history',sprintf('Resize to Nmax = %d',Nmax2));
            elseif Nmax2 < obj.Nmax % shrink the size
                for ii = 1:obj.Nmonth
                    SH(ii) = cSH(obj(ii).C(1:Nmax2+1,1:Nmax2+1),obj(ii).S(1:Nmax2+1,1:Nmax2+1),obj(ii).tt);
                    SH(ii) = SH(ii).setinfo(obj(ii));
                end
                SH = SH.setinfo('history',sprintf('Resize to Nmax = %d',Nmax2));
            else % size is not changed;
                SH= obj;
            end
        end
        
        function SH=extract(SH0,tran)
            % {{SH; extract; time; resize}}
            % SH0.tt, SH0.C, SH0.S
            % trange : 1. a single time ( [year,month] )
            %          2. time range ( [year1,month1,yaer2,month2] )
            
            tt = SH0.optt;
            
            if length(tran) == 2
                [minval,minloc] = min( abs( tt-generate_tt(tran) ) );
                SH = SH0(minloc);
                if minval ~= 0
                    fprintf('minloc= %d; min distances are %7.4f days\n',minloc,minval*365);
                end
            elseif length(tran) == 4
                t1 = generate_tt(tran(1:2))-0.5/12; % first day of the month
                t2 = generate_tt(tran(3:4))+0.5/12; % end day of the month
                
                loc_t1 = sum(tt<t1)+1;
                loc_t2 = sum(tt<=t2);
                
                fprintf('loc_t1 = %d, loc_t2 = %d\n',loc_t1,loc_t2);
                
                SH = SH0(loc_t1:loc_t2);
                
            else
                error('length(trange) = %d',length(tran));
            end
        end
        
        function SH=remove_by_ym(obj, tlist)
            % {{SH; remove; reshape}}
            % SH = remove(SH0, [2003,1; 2004,12] );
            
            Ntt = size(tlist,1);
            if size(tlist,2)~=2
                mydisp('error: the form of input is [y1,m1; y2, m2; ...]\n');
            end
            
            tts = obj.optt;
            iloc = [];
            for ii = 1:Ntt
                tt0 = time_transfer([tlist(ii,:),15],1);
                itmp = find(abs(tts - tt0) < 0.5/12 );
                if ~isempty(itmp)
                    iloc = [iloc, itmp];
                end
            end
            
            if ~isempty(iloc)
                fprintf('These sequence has been removed: ');
                fprintf('%d ',iloc);
                fprintf('\n');
                
                iseq0 = 1:obj.Nmonth;
                iseq1 = setdiff(iseq0, iloc);
                
                SH = obj(iseq1);
                SH = SH.setinfo('history',sprintf('remove_by_ym: %s',mat2str(tlist)));
            else
                SH = obj;
                fprintf('No sequence has been removed in cSH.remove\n');
            end
        end
        
        function SH=remove_by_seq(obj, iloc)
            % {{SH; remove; reshape}}
            % SH = SH.remove_by_seq([1,3,5]);
            
            if ~isempty(iloc)
                fprintf('These sequence has been removed: ');
                fprintf('%d ',iloc);
                fprintf('\n');
                
                iseq0 = 1:obj.Nmonth;
                iseq1 = setdiff(iseq0, iloc);
                
                SH = obj(iseq1);
                if numel(iloc) > 10
                    ss = sprintf('%s, ..., %s',mat2str(iloc(1:5)),mat2str(iloc(end-4:end)));
                else
                    ss = mat2str(iloc);
                end
                SH = SH.setinfo('history',sprintf('remove_by_seq: %s',ss));
            else
                fprintf('No sequence has been removed in cSH.remove\n');
            end
        end
        
        %% output
        
        function tt = optt(obj) % extract by time
            tt = [obj.tt]';
        end
        
        function ser = nm(obj,n,m) % extract by degree/order
            % {{SH; extract; degree; order}}
            ser = zeros(obj.Nmonth,numel(n)); % row: time; col: C, S;
            
            if numel(n) ~= numel(m)
                error('numel(n) = %d, but numel(m) = %d\n',numel(n),numel(m));
            end
            
            for ii = 1:obj.Nmonth
                for inm = 1:numel(n)
                    if m(inm) >= 0
                        ser(ii,inm) = obj(ii).C(n(inm)+1,m(inm)+1);
                    else
                        ser(ii,inm) = obj(ii).S(n(inm)+1,-m(inm)+1);
                    end
                end
            end
            
        end
        
        function ser = degree(obj,n) % extract the whole degree
            % {{degree; SH}}
            ser = zeros(obj.Nmonth,2*n+1); % row: time; col: Cn0,...,Cnn, Sn1,...,Snn;
            for ii = 0:n
                ser(:,ii+1) = obj.nm(n,ii);
                if ii > 0
                    ser(:,ii+1+n) = obj.nm(n,-ii);
                end
            end
        end
        
        function output2file(SH,op_file)
            
            fid=fopen(op_file,'w');
            if fid==-1
                error('error open file: %s',op_file);
            end
            
            if SH.Nmonth > 1
                fprintf('More than one month, only the first month will be written\n');
            end
            for in = 0:SH.Nmax
                for im = 0:in
                    fprintf(fid,'%5d%5d%23.15e%23.15e\n',in,im,SH(1).C(in+1,im+1),SH(1).S(in+1,im+1));
                end
            end
            
            fclose(fid);
        end
        
        function save(obj,opfile)
            % {{save; object; struct}}
            
            SH = classToStruct(obj);
            
            save(opfile,'SH');
            
        end
        
        function save_nc(obj,sfile)
            %{
write cSH in the format of netcdf
% {{netcdf; SH; write}}
            %}
            
            if nargin == 0
                fun_ex1;
                return;
            end
            % ----- end of head -----
            Nmax = obj.Nmax;
            Nmonth = obj.Nmonth;
            [rC,rS,tt] = obj.reformat('matrix');
            nccreate(sfile,'C',...
                'Dimensions', {'N',Nmax+1,'N',Nmax+1,'Ntime',Nmonth},...
                'FillValue','disable');
            ncwrite(sfile,'C',rC);
            ncwriteatt(sfile,'C','description','C(n,m,t)');
            
            nccreate(sfile,'S',...
                'Dimensions', {'N',Nmax+1,'N',Nmax+1,'Ntime',Nmonth},...
                'FillValue','disable');
            ncwrite(sfile,'S',rS);
            ncwriteatt(sfile,'S','description','S(n,m,t)');
            
            nccreate(sfile,'time',...
                'Dimensions', {'Ntime',Nmonth},...
                'FillValue','disable');
            ncwrite(sfile,'time',tt);
            ncwriteatt(sfile,'time','description','Time in decimal years');
            
            ncwriteatt(sfile,'/','creation_date',datestr(now));
            
        end
        
        
        %% plot
        function plot_t(obj,o_symbol)
            % {{SH; plot; time; tt}}
            if nargin == 1
                o_symbol = 'rs';
            end
            % tt=load_tt('csr_full');
            
            tt = obj.optt;
            
            if sum(isnan(tt)) > 0
                fprintf('There is NaN in the time: \n');
                fprintf('%.2f ',tt);
                error(' ');
            end
            yearlist = floor(tt);
            monthlist = (tt-yearlist) * 12;
            
            plot(yearlist,monthlist,o_symbol,'markersize',15,'markerfacecolor','none');
            for ii = min(yearlist):max(yearlist)
                count = sum(yearlist == ii);
                text(ii,12.5,sprintf('%d',count),'color','k','horizontalA','center');
            end
            
            % if nargin == 1
            set(gca,'xtick',min(yearlist(:)):max(yearlist(:)))
            axis([min(yearlist(:))-1,max(yearlist(:))+1,0,12])
            % end
            for ii = 1:12
                sticklabel{ii} = sprintf('%d',ii);
            end
            set(gca,'ytick',0.5:11.5,'yticklabel',sticklabel);
            set(gca,'ygrid','on')
            
            % output number
            Nmonth = numel(tt);
            for imon = 1:Nmonth
                text(yearlist(imon)+0.1,monthlist(imon),num2str(imon));
            end
        end
        
        function LLZ=preview(SH,varargin)
            % {{SH; plot; preview}}
            PAR0 = struct('to', 'ewh','resolution',1,...
                'iGauss', 00,'iPxMy',0,'iFan',0,'iddk',0,...
                'ifrom', 0,...
                'icran',[],...
                'iplot',1,'new',1);
            
            PAR = var_initial(PAR0,varargin);
            
            if PAR.new == 1
                myfigure;
            end
            
            % PAR = rmfield(PAR,'ii');
            SH = CS_smooth(SH(1),PAR);
            [LLZ]=LLZ_forward_ns(SH,'iplot',1,'to',PAR.to);
            if ~isempty(PAR.icran)
                mycaxis(PAR.icran);
            else
                mycaxis(LLZ);
            end
            
            if PAR.iGauss > 0
                fprintf('G%d is used for SH.preview\n',PAR.iGauss);
            end
            
        end % 300 km Gaussian, global, plot
        
        function varargout=kaula(SH,o_symb,inormalize)
            % {{SH; koala; plot}}
            
            if nargin < 2
                o_symb = NaN;
            end
            
            if nargin < 3
                inormalize = 0;
            end
            
            
            SH = SH(1);
            rC = SH.C;
            rS = SH.S;
            val = rC.^2+rS.^2;
            Nmax = SH.Nmax;
            nlist = (0:Nmax)'*2+1;
            val = sqrt( sum(val,2)./nlist ); % sum along with im
            
            % if rC is too small , set rC = rC*R_e
            if inormalize == 1 && max(val) < 1D-6
                val = val*6371D6; %mm
                ystr = 'mm of geoid';
            else
                ystr = 'absolution diff';
            end
            nn = 0:(size(rC,1)-1); nn = nn(:);
            
            
            if ~isnan(o_symb)
                h = semilogy(nn,val,o_symb);
                xlabel('degree');
                ylabel(ystr);
                grid on;
            end
            
            if nargout == 1
                varargout = {val};
            elseif nargout == 2
                varargout = {val,nn};
            elseif nargout == 3
                varargout = {val,nn,h};
            end
        end
        
        function [xx,yy,ccss1] = plot(SH1,iref)
            % {{SH; plot}}
            
            if SH1.Nmonth > 1
                fprintf('Only the first month will be plotted\n');
            end
            
            SH1 = SH1(1);
            Nmax = SH1.Nmax;
            cc = SH1.C;
            ss = SH1.S;
            ss = ss(:,2:end);
            ccss = [fliplr(ss),cc];
            ccss(ccss==0) = NaN;
            if nargin == 1
                iref = -ceil(log10(max(abs(ccss(~isnan(ccss))))));
            end
            ccss1 = (iref+log10(abs(ccss))).*sign(ccss);
            
            nn = 0:Nmax;
            mm = -Nmax:Nmax;
            [xx, yy] = meshgrid(mm,nn);
            
            
            mypcolor(xx,yy,ccss1);
            
            xlabel('Order (m)');
            ylabel('Degree (n)');
            title(sprintf('(%d+log10(value))*sign(value)',iref))
            ind = ~isnan(ccss1);
            vmin = floor(min(ccss1(ind)));
            vmax = ceil(max(ccss1(ind)));
            caxis([vmin,vmax]);
            grid on
            % mycpt('polar_inverse','icolormap',1);
            mycpt('polar','icolormap',1);
            
        end % log10(val) * sign(val)
        
        function [xx,yy,ccss1] = plot2(SH1)
            % {{SH; plot}}
            
            if SH1.Nmonth > 1
                fprintf('Only the first month will be plotted\n');
            end
            
            
            SH1 = SH1(1);
            Nmax = SH1.Nmax;
            cc = SH1.C;
            ss = SH1.S;
            ss = ss(:,2:end);
            ccss = [fliplr(ss),cc];
            ccss(ccss==0) = NaN;
            ccss1 = log10(abs(ccss));
            ccss2 = sign(ccss);
            nn = 0:Nmax;
            mm = -Nmax:Nmax;
            [xx, yy] = meshgrid(mm,nn);
            
            
            subplot(2,1,1)
            mypcolor(xx,yy,ccss1);
            xlabel('Order (m)');
            ylabel('Degree (n)');
            title('log10(value)')
            ind = ~isnan(ccss1);
            vmin = floor(min(ccss1(ind)));
            vmax = ceil(max(ccss1(ind)));
            caxis([vmin,vmax]);
            grid on
            
            subplot(2,1,2)
            mypcolor(xx,yy,ccss2);
            xlabel('Order (m)');
            ylabel('Degree (n)');
            title('sign(value)')
            grid on
            colormap jet
            
        end % log10(val) + sign(val)
        
        function [xx,yy,amp,phz] = plot_polar(SH1, imon, iplot)
            % {{SH; plot}}
            if nargin < 2
                imon = 1;
                if SH1.Nmonth > 1
                    fprintf('Only the first month will be plotted\n');
                end
            end
            
            if nargin < 3
                iplot = 1;
            end
            
            SH1 = SH1(imon);
            Nmax = SH1.Nmax;
            cc = SH1.C;
            ss = SH1.S;
            amp = sqrt(cc.^2+ss.^2);
            ind_nan = amp == 0;
            amp(ind_nan) = NaN;
            amp = log10(amp);
            phz = atan2d(ss,cc);
            phz(ind_nan) = NaN;
            
            nn = 0:Nmax;
            mm = 0:Nmax;
            [xx, yy] = meshgrid(mm,nn);
            
            if iplot == 1
                figure;
                subplot(1,2,1);
                mypcolor(xx,yy,amp);
                xlabel('Order (m)');
                ylabel('Degree (n)');
                title('Amplitude')
                caxis([-15,-9]);
                set(gca,'ydir','reverse');
                % ind = ~isnan(ccss1);
                % vmin = min(ccss1(ind));
                % vmax = max(ccss1(ind));
                % caxis([vmin,vmax]);
                
                subplot(1,2,2)
                mypcolor(xx,yy,phz);
                xlabel('Order (m)');
                ylabel('Degree (n)');
                title('Phase')
                caxis([-180,180]);
                set(gca,'ydir','reverse');
                
                colormap hsv
            end
        end
        
        function plot_abs(SH1, imon)
            % {{SH; plot}}
            if nargin == 1
                imon = 1;
                if SH1.Nmonth > 1
                    fprintf('Only the first month will be plotted\n');
                end
            end
            SH1 = SH1(imon);
            Nmax = SH1.Nmax;
            cc = SH1.C;
            ss = SH1.S;
            ss = ss(:,2:end);
            ccss = [fliplr(ss),cc];
            ccss(ccss==0) = NaN;
            ccss1 = abs(ccss);
            % ccss2 = sign(ccss);
            nn = 0:Nmax;
            mm = -Nmax:Nmax;
            [xx, yy] = meshgrid(mm,nn);
            
            mypcolor(xx,yy,ccss1);
            xlabel('Order (m)');
            ylabel('Degree (n)');
            title('Absolute value (white box indicates positive)')
            % ind = ~isnan(ccss1);
            % vmin = min(ccss1(ind));
            % vmax = max(ccss1(ind));
            % caxis([vmin,vmax]);
            
            hold on;
            ind = ccss > 0;
            % plot(xx(ind),yy(ind),'ws','markerfacecolor','w','markers',10);
            
            plot(xx(ind),yy(ind),'ws','markersize',5);
            
            hold off;
            grid on
            
            colormap jet
        end% absolute(val)
        
        function H = plot_nm(obj,n,m)
            % {{SH; plot; single}}
            ser = obj.nm(n,m);
            str = sprintf(' (%d,%d)',n,abs(m));
            
            if m < 0
                H = plot(obj.optt,ser,'o-');
                legend(['S',str]);
            else
                H = plot(obj.optt,ser(:,1),'o-');
                legend(['C',str]);
            end
            grid on;
        end
        
        function plot_n(obj,n)
            % {{SH; plot; degree}}
            if obj.Nmonth > 1
                fprintf('only the first month will be used\n');
                obj = obj(1);
            end
            xx = 0:n;
            plot(xx,obj.C(n+1,xx+1),'o-',xx,obj.S(n+1,xx+1),'o-');
            legend('C','S');
            xlabel('order (m)');
            grid on;
        end
        
        function plot_m(obj,m)
            % {{SH; plot; order}}
            if obj.Nmonth > 1
                fprintf('only the first month will be used\n');
                obj = obj(1);
            end
            xx = m:obj.Nmax;
            plot(xx,obj.C(xx+1,m+1),'o-',xx,obj.S(xx+1,m+1),'o-');
            legend('C','S');
            xlabel('degree (n)');
            grid on;
        end
        
        %% calculation
        function obj3 = minus(obj1,obj2)
            % {{SH; minus}}
            if obj1.Nmax ~= obj2.Nmax
                error('The Nmax of two inputs are %d and %d',obj1.Nmax,obj2.Nmax);
            end
            
            rC = zeros(obj1.Nmax+1,obj1.Nmax+1,obj1.Nmonth);
            rS = rC;
            
            if obj2.Nmonth ~= 1
                if obj1.Nmonth == obj2.Nmonth
                    for ii = 1:obj1.Nmonth
                        rC(:,:,ii) = obj1(ii).C - obj2(ii).C;
                        rS(:,:,ii) = obj1(ii).S - obj2(ii).S;
                    end
                    obj3 = cSH(rC,rS,obj1.optt);
                    obj3 = obj3.setinfo(obj1);
                    %         obj3 = obj3.setinfo(...
                    %     'history',sprintf('%s MINUS %s',inputname(1),inputname(2)));
                else
                    error('The sizes of two inputs are %d and %d',obj1.Nmonth,obj2.Nmonth);
                end
            else
                for ii = 1:obj1.Nmonth
                    rC(:,:,ii) = obj1(ii).C - obj2.C;
                    rS(:,:,ii) = obj1(ii).S - obj2.S;
                end
                obj3 = cSH(rC,rS,obj1.optt);
                obj3 = obj3.setinfo(obj1.info);
                %     obj3 = obj3.setinfo('history',sprintf('%s MINUS %s',inputname(1),inputname(2)));
            end
        end
        
        function obj3 = plus(obj1,obj2)
            if obj1.Nmax ~= obj2.Nmax
                error('The Nmax of two inputs are %d and %d',obj1.Nmax,obj2.Nmax);
            end
            
            rC = zeros(obj1.Nmax+1,obj1.Nmax+1,obj1.Nmonth);
            rS = rC;
            
            if obj2.Nmonth ~= 1
                if obj1.Nmonth == obj2.Nmonth
                    for ii = 1:obj1.Nmonth
                        rC(:,:,ii) = obj1(ii).C + obj2(ii).C;
                        rS(:,:,ii) = obj1(ii).S + obj2(ii).S;
                    end
                    obj3 = cSH(rC,rS,obj1.optt);
                    obj3 = obj3.setinfo(obj1);
                    %         obj3 = obj3.setinfo('history',sprintf('SH plus'));
                    
                else
                    error('The sizes of two inputs are %d and %d',obj1.Nmonth,obj2.Nmonth);
                end
            else
                for ii = 1:obj1.Nmonth
                    rC(:,:,ii) = obj1(ii).C + obj2.C;
                    rS(:,:,ii) = obj1(ii).S + obj2.S;
                end
                obj3 = cSH(rC,rS,obj1.optt);
                obj3 = obj3.setinfo(obj1);
                %     obj3 = obj3.setinfo('history',sprintf('SH plus'));
                
            end
        end
        
        function obj1 = times(obj1,a)
            % {{SH; times}}
            
            for ii = 1:numel(obj1)
                obj1(ii).C = obj1(ii).C * a;
                obj1(ii).S = obj1(ii).S * a;
            end
            
            obj1 = obj1.setinfo('history',sprintf('multiplied by %f',a));
        end
        
        function SH_mean=mean(obj,o_trange)
            % {{SH; mean}}
            if nargin == 2
                SH = obj.extract(o_trange);
            else
                SH = obj;
            end
            
            [rCs,rSs,tt]=SH.reformat('matrix');
            
            rCmean = mean(rCs,3);
            rSmean = mean(rSs,3);
            
            SH_mean = cSH(rCmean,rSmean,mean(tt));
        end
        
        function SH=demean(obj, o_trange)
            % {{SH; mean; demean}}
            if nargin == 2
                SHm = obj.mean(o_trange);
            else
                SHm = obj.mean;
                o_trange = obj.trange;
            end
            
            SH = obj - SHm;
            SH = SH.setinfo('history',sprintf('demean: %s',num2str(o_trange)));
        end
        
        %% compare
        function iequal = isequal(obj1,obj2)
            % {{SH; equal}}
            dSH = obj1 - obj2;
            vCS = dSH(1).reformat('vCS');
            
            max_diff = max(abs( vCS(:) ));
            
            if max_diff == 0
                iequal = 1;
            else
                mydisp('max diff = %g\n',max_diff);
                iequal = 0;
            end
        end
    end
    
    methods (Static)
%         [SH_idx]=get_idx(Nmax, nlist, mlist);
        function [SH_idx]=get_idx(Nmax, nlist, mlist)
            % return index by rules:
            % either nlist or mlist supports five types of input:
            % 1) '<20'
            % 2) '>10'
            % 3) '=5'
            % 4) combination of 1&2, '<20 >10', sequence does not matter
            % 5) combination of 3, '=1 =3 =5'
            % 6) empty or non-char inputs, then no filter is applied
            % 7) mlist has an extra input of 'abs', usually, m is negative in Snm
            %    and positive in Cnm, if 'abs' is set, the abs value of m will be used
            %
            % SH_idx: index in the form of sSH
            %
            % optional input:
            %   abs,
            %
            % See also: setval
            %
            % {{cSH; index}}
            
            if nargin == 0
                help(mfilename);
                fun_ex1;
                return;
            end
            
            % PAR0 = struct('abs_m',1);
            % PAR = var_initial(PAR0, varargin);
            
            % --------------------end of head-----------
            SH1 = cSH(Nmax,NaN);
            [~,LIST0] = SH1.reformat('vCS');
            
            if ischar(mlist) && contains(mlist,'abs')
                LIST = abs(LIST0);
            else
                LIST = LIST0;
            end
            
            if ischar(mlist)
                if contains(mlist,'=')
                    a = get_var(mlist,'=*','f',0); % two patterns are matched
                    idx1 = LIST(:,2)>Inf; % all false
                    if ~isempty(a)
                        for ii = 1:numel(a)
                            idx1 = LIST(:,2)==a(ii) | idx1;
                        end
                    end
                    idxm = idx1;
                else
                    a = get_var(mlist,'>*','f',0); % two patterns are matched
                    if ~isempty(a)
                        idx1 = LIST(:,2)>a;
                    else
                        idx1 = LIST(:,2)>-Inf; % all true
                    end
                    a = get_var(mlist,'<*','f',0); % two patterns are matched
                    if ~isempty(a)
                        idx2 = LIST(:,2)<a;
                    else
                        idx2 = LIST(:,2)<Inf; % all true
                    end
                    idxm = idx1 & idx2;
                    
                end
            else
                idxm = LIST(:,2)>-Inf; % all true
            end
            
            if ischar(nlist)
                if contains(nlist,'=')
                    a = get_var(nlist,'=*','f',0); % two patterns are matched
                    idx1 = LIST(:,1)>Inf; % all false
                    if ~isempty(a)
                        for ii = 1:numel(a)
                            idx1 = LIST(:,1)==a(ii) | idx1;
                        end
                    end
                    idxn = idx1;
                else
                    a = get_var(nlist,'>*','f',0); % two patterns are matched
                    if ~isempty(a)
                        idx1 = LIST(:,1)>a;
                    else
                        idx1 = LIST(:,1)>-Inf; % all true
                    end
                    a = get_var(nlist,'<*','f',0); % two patterns are matched
                    if ~isempty(a)
                        idx2 = LIST(:,1)<a;
                    else
                        idx2 = LIST(:,1)<Inf; % all true
                    end
                    idxn = idx1 & idx2;
                    
                end
            else
                idxn = LIST(:,1)>-Inf; % all true
            end
            
            idx = idxm & idxn;
            
            SH_idx = cSH(idx,LIST0,NaN,'vCS');
        end

%         obj = load(ipfile);
        
        function ex1() % 6 ways to create cSH
            
            % 1.
            SH1 = cSH(60,1:3);
            SH1 = SH1.setinfo('name','create_by_Nmax');
            
            % 2.
            SH0 = load_SH('csr06','trange',[2003,1,2003,12]);
            SH0 = cSH(SH0);
            lmcosi = SH0(1).reformat('lmcosi');
            SH2 = cSH(lmcosi,NaN,'lmcosi');
            SH2 = SH2.setinfo('name','create_by_lmcosi');
            SH0(1).isequal(SH2)
            
            % 3.
            [rCs,rSs,tt] = SH0.reformat('matrix');
            SH3 = cSH(rCs,rSs,tt);
            SH3 = SH3.setinfo('name','create_by_matrix');
            isequal(SH0,SH3)
            
            % 4.
            [vCS,LIST] = SH0.reformat('vCS');
            SH4 = cSH(vCS,LIST,[SH0.tt],'vCS');
            SH4 = SH4.setinfo('name','create_by_vCS');
            isequal(SH0,SH4)
            
            % 5.
            cs = SH0.reformat('cs');
            SH5 = cSH(cs,[SH0.tt],'cs');
            SH5 = SH5.setinfo(SH0(1));
            SH5 = SH5.setinfo('name','create_by_cs');
            isequal(SH0,SH5)
            
            % 6.
            cs = SH0.reformat('sc');
            SH6 = cSH(cs,[SH0.tt],'sc');
            SH6 = SH6.setinfo('name','create_by_sc');
            isequal(SH0,SH6)
            
        end
        function ex2() % setinfo, print, find, save
            
            a = cSH(60,1);
            a = a.setinfo('name','test1','type','GRACE_1','history','cSH, ex2');
            a.print;
            a(2) = cSH(60,2);
            a(2) = a(2).setinfo(struct('name','test2','type','GRACE_1'));
            a(3) = cSH(60,3);
            a(3) = a(3).setinfo(a(2));
            b = a.find('type','GRACE_1');
            b(3).print
            a.save('tmp.mat');
        end
        function obj = ex0
            obj = load_SH('csr06');
        end
    end
end

%%
% function for method 'vCS'
function [vCS,LIST]=SH_to1D(SHs,o_in0)

if nargin == 1
    o_in0 = 0;
end
%
% Ncs = (Nmax+1)^2 - (in0)^2;

[rCs,rSs]=SHs.reformat('matrix');


[vC,LIST1]=rC_2Dto1D(rCs,o_in0);
[vS,LIST2]=rS_2Dto1D(rSs,o_in0);

vCS = [vC;vS];
LIST = [LIST1;LIST2];

end

function [vC,list]=rC_2Dto1D(rCs,in0)

Nmax = size(rCs,1) - 1;
Nmonth = size(rCs,3);

Nc = (Nmax+1)*(Nmax+2)/2 - (in0)*(in0+1)/2;
vC( Nc, Nmonth ) =0;
list(Nc,2) = 0;
% list(:,1) = 1; % 1 stands for rC
ik = 0;
for in = in0:Nmax
    for im = 0:in
        ik = ik+1;
        vC(ik,:) = rCs(in+1,im+1,:);
        list(ik,1:2) = [in,im];
    end
end
end

function [vS,list]=rS_2Dto1D(rSs,in0)

if numel(rSs) == 1 % Nmax == 0;
    vS = [];
    list = [];
    return
end

if in0<1
    in0 = 1;
end

Nmax = size(rSs,1) - 1;
Nmonth = size(rSs,3);

Nc = (Nmax)*(Nmax+1)/2 - (in0-1)*(in0)/2;
vS( Nc, Nmonth ) =0;
list(Nc,2) = 0;
% list(:,1) = 2; % 2 stands for rS
ik = 0;
for in = in0:Nmax
    for im = 1:in
        ik = ik+1;
        vS(ik,:) = rSs(in+1,im+1,:);
        list(ik,1:2) = [in,-im];
    end
end

end

