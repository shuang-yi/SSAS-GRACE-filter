function [PAR,nvar] = var_initial(PAR0,varargin)
% initialization of input parameters
% === Usage 1:
% PAR0 = struct('iplot',0); 
% PAR = var_initial(PAR0, varargin);
%
% varargin could be {1, 2, 'a', 3}, then only PAR.a = 3 will be recognized,
%   and num_none_par = 2
%
% === Usage 2:
% PAR = var_initial(par0, par);
% par0: default value of each parameter, all parameters should be defined
% par: new value of each parameter, if omitted, the default value will be used.
% 
% Usage: 
% var_initial(struct('a',1, 'b',2, 'c','a;;b;;c'),struct('a',11, 'c', 'c'))
%   list of supported values are seperated by ';;', default is the first one
%
% var_initial(struct('c_a',1, 'c_name', []),struct('c_a',{1,2}, 'c_name', {'a','b'}))
%   variables starting with 'c_' will be expanded into individual settings
% See also: 

% {{function; var; input; parameter}}



if nargin == 0  
    fun_ex1;  
    return;
end
% ----- end of head -----

if nargin == 0
%     help(mfilename);
%     return;
% elseif nargin == 1
%     error('%s: input should >2',mfilename);
end

if nargin == 2
    if iscell(varargin{1})
        varargin = varargin{1};
    
    end
end

if isempty(varargin)
    PAR = var_ini(PAR0,{});
    return
end

if isstruct(varargin{1}) % input is PAR0, PAR
    PAR = varargin{1};
    for ii = 1:numel(PAR)
        PAR2(ii) = var_ini(PAR0,PAR(ii));
    end
    PAR = PAR2;
    
    return;
end

nvar = numel(varargin);

if mod(nvar,2) == 1
    error('The number of input is %d, should be even',nvar);
end

for ii = 2:2:nvar
    if iscell(varargin{ii})
        % a cell input will be expanded into matrix by MATALB automatically
        % , but usually we do not want it happens, unless the fieldname
        % begins with 'c_'
        name0 = varargin{ii-1};
        if ~(numel(name0) > 2 && strcmp(name0(1:2),'c_'))
            % , here protect it, if its name does not start with 'c_';
            varargin(ii) = {varargin(ii)};
        end
    end
end

PAR = struct(varargin{1:nvar});
for ii = 1:numel(PAR)
    PAR2(ii) = var_ini(PAR0,PAR(ii));
end
PAR = PAR2;
end

function par_new = var_ini(par0,par)
% parse structured inputs


par_new = par0;
sfield0 = fieldnames(par0);
% for ii = 1:numel(sfield0) %
%     if ~isfield(par,sfield{ii})
%         par.(sfield{ii}) = par0.(sfield{ii});
%     end

for ii = 1:numel(sfield0)
    if ischar(par0.(sfield0{ii})) && contains(par0.(sfield0{ii}),';;') % if it is a parameter with a list (separated by ;)
        slist = regexp(par0.(sfield0{ii}),';;','split');
        if isfield(par,sfield0{ii}) % it is defined
            ind = strcmp(slist,par.(sfield0{ii}));
            if sum(ind) > 1
                error(' More than one parameter value is in the list of %s', sfield0{ii});
            elseif sum(ind) == 0
                fprintf('%s (''%s'') can only have one value of:',sfield0{ii},par.(sfield0{ii}));
                fprintf(' ''%s'' ',slist{:}); fprintf('\n');
                error(' parameter value not in the list');
            else
                par_new.(sfield0{ii}) = slist{ind};
            end
        else % defaultly use the first one
            par_new.(sfield0{ii}) = slist{1};
        end
    else % non-list parameter
        if isfield(par,sfield0{ii}) %
            par_new.(sfield0{ii}) = par.(sfield0{ii});
        else % default value
            par_new.(sfield0{ii}) = par0.(sfield0{ii});
        end
    end
end

% make sure par and par0 share exactly the same fields {
if ~isempty(par)
    ik = 0;
    sfield2 = fieldnames(par);
    for ii = 1:numel(sfield2)
        if ~any(strcmp(sfield0,sfield2{ii}))
            fprintf('*error* ''%s'' cannot be found in PAR0\n',sfield2{ii});
            ik = 1;
        end
    end
    if ik == 1
        error('error in input');
    end
end
end

%%
function fun_ex1
PAR0 = struct('iplot',0,'stype','aa;;bb','var2','xx;;yy'); 

disp(PAR0);
PAR = var_initial(PAR0, 'iplot',1,'stype','bb');
disp(PAR)
PAR = var_initial(PAR0, {'iplot',1,'stype','bb'});
disp(PAR)
PAR = var_initial(PAR0, struct('iplot',1,'stype','bb'));
disp(PAR)
par2 = var_initial(struct('c_a',1, 'c_name', []),struct('c_a',{1,2}, 'c_name', {'a','b'}));
disp(par2(1));
disp(par2(2));
end