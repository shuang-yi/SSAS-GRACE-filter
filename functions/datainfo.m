classdef datainfo
    % datainfo('name','dataname','unit','mm/yr','type','GRACE-FO','history','downloaded at
    %   xxx');
    % see examples:
    %  datainfo.ex1: create a simple datainfo data, and then reset the info
    %  datainfo.ex2: how to find data
    %
    % {{datainfo; class}}
    properties
        name
        unit
        type
        history
    end
    
    methods
        % create {
        function obj = datainfo(varargin)
            
            PAR0 = struct(...
                'unit',[],'name',[],'type',[],...
                'history',[]);
            
            [PAR] = var_initial2(PAR0,varargin);
            
            obj.name = PAR.name;
            obj.unit = PAR.unit;
            obj.type = PAR.type;
            
            %             obj.history{1} = sprintf('Created on %04d-%02d-%02d, %02d:%02d:%02.0f',clock);
            
            if ~isempty(PAR.history)
                obj = obj.set('history',PAR.history);
            else
                obj.history = {};
            end
        end
        % create }
        
        % multiple msgs should be stored in cell
        function obj = set(obj,varargin)
            % multiple history should be stored in cell
            
            % {{datainfo; set; info}}
            
            % add: append history, otherwise replace history
            
            PAR0 = struct('history',[],'name',[],'unit',[],'type',[],'add',1);
            
            if isa(varargin{1},'datainfo')
                varargin{1} = classToStruct(varargin{1});
            end
            
            PAR = var_initial(PAR0,varargin);
            
            if ~isempty(PAR.name)
                obj.name = PAR.name;
            end
            
            if ~isempty(PAR.unit)
                obj.unit = PAR.unit;
            end
            
            if ~isempty(PAR.type)
                obj.type = PAR.type;
            end
            
            if ~isempty(PAR.history)
                if iscell(PAR.history) % input is composed of cell data
                    if PAR.add == 1
                        
                        if isempty(obj.history)
                            n0 = 0;
                            obj.history = {}; % make sure it is a cell format
                        else
                            n0 = numel(obj.history);
                        end
                        N = numel(PAR.history);
                        obj.history(n0+1:n0+N) = PAR.history;
                        
                    else % replace
                        obj.history = PAR.history;
                    end
                else  % input is composed of string data
                    if PAR.add == 1
                        if isempty(obj.history)
                            n0 = 0;
                            obj.history = {}; % make sure it is a cell format
                        else
                            n0 = numel(obj.history);
                        end
                        obj.history{n0+1} = PAR.history;
                    else % replace
                        obj.history{1} = PAR.history;
                    end
                end
                obj.history = obj.history(:);
            end
            
        end
        
        % print the datainfo
        function print(obj)
            % {{datainfo; dispinfo; print}}
            
            msg = obj.history;
            if ~iscell(msg)
                msg = {msg};
            end
            
            fprintf('\nName: %s;   Unit: %s;   Type: %s\n',obj.name,obj.unit,obj.type);
            
            
            fprintf('\n--- history ---\n');
            for ii = 1:numel(msg)
                fprintf('  (%d) %s\n',ii, msg{ii});
            end
            fprintf('--- end of history ---\n\n');
            
            
        end
        
        % for name & unit, is true if it matches;
        % for history, is true if it contains
        function id = find(obj,varargin)
            
            % {{datainfo; filter; find}}
            
            
            PAR0 = struct('history',[],'name',[],'type',[],'unit',[]);
            PAR = var_initial(PAR0,varargin);
            
            Nobj = numel(obj);
            
            if ~isempty(PAR.name)
                ss = {obj.name};
                id1 = strcmp(ss, PAR.name);
            else
                id1 = true(1,Nobj);
            end
            
            if ~isempty(PAR.unit)
                ss = {obj.unit};
                id2 = strcmp(ss, PAR.unit);
            else
                id2 = true(1,Nobj);
            end
            
            if ~isempty(PAR.type)
                ss = {obj.type};
                id4 = strcmp(ss, PAR.type);
            else
                id4 = true(1,Nobj);
            end
            
            
            
            if ~isempty(PAR.history)
                for ii = 1:Nobj
                    ss{ii} = strjoin(obj(ii).history,',');
                end
                id3 = contains(ss,PAR.history);
            else
                id3 = true(1,Nobj);
            end
            
            id = id1&id2&id3&id4;
            
        end
    end
    
    methods (Static)
        function ex1 % create a data, and then reset the info
            % 1. create 
            obj = datainfo('name','data1','unit','m','type','type-a','history',{'an example','second line'});
            obj.print;
            
            % 2. reset the info
            s = struct('name','data2','type','type-b');
            obj = obj.set(s);
            obj.print;
        end
        
        function ex2 % example of find
            obj(1) = datainfo('name','data1','unit','m','history','type: type1');
            obj(2) = datainfo('name','data2','unit','mm','history','type: type1');
            obj(3) = datainfo('name','data3','unit','m','type','a','history','type: type2');
            
            fprintf('### The indices show the objectives matching the conditions  ###\n')
            id = obj.find('unit','m')
            id = obj.find('name','data2')
            id = obj.find('history','type1')
            id = obj.find('type','a')
            
        end
    end
    
end