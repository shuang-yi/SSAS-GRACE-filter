function ss2=expand_format(ss,varargin)
% expand '%f@2 %d@3' to '%f %f %d %d %d'
% {{write; auto; expand}}
% See also

if nargin == 0
    help(mfilename);
    fun_ex1;
    return;
end

PAR0 = struct('ex',0);
PAR = var_initial(PAR0, varargin);

% --------------------end of head-----------

sss = regexp(ss,' ','split');
ss2 = [];
for ii = 1:numel(sss)
    if ~isempty(sss{ii})
        a = regexp(sss{ii},'@','split');
        
        if numel(a)==1
            ss2 = [ss2,a{1},' '];
        elseif numel(a) == 2
            itime = sscanf(a{2},'%d');
            if isempty(itime)
                error('cannot recognized a number from "%s"\n',a{2});
            end
            for ij = 1:itime
                ss2 = [ss2,a{1},' '];
            end
        else
            error('more than one @ is found in "%s"\n',sss{ii});
        end
    end
end
ss2 = ss2(1:end-1);
end

%% example

function fun_ex1()
expand_format('%f@2 %d@3 %*.2f@2')
end