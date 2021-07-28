function [val] = get_var(str,spattern0, stype, iwarn)
% {{input; var; read; parameter}}
% spattern 'a*b'
% stype: d, f, s
% iwarn: if 1, gives a warn if no values found.
% 
% example:
%  get_var('a2004_2005','','d') returns 2004 and 2005, spatterns can also
%    be '[a_]*'

if nargin == 0
    fun_ex1;
    return;
end

if strcmp(stype,'d')
    sformat = '(\d*)';
elseif strcmp(stype,'f')
    sformat = '([-0-9.]*)';
elseif strcmp(stype,'s')
    sformat = '(\w*)';
else
    mydisp('error: Only support d,f,s in the data type');
end

if nargin < 4 % default there is a warn
    iwarn = 1;
end

if isempty(regexp(spattern0,'\*','once')) % default, spattern0 is a prefix
    spattern0 = [spattern0,'*'];
end
spattern = strrep(spattern0, '*', sformat);
a = regexp(str,spattern,'tokens');
if isempty(a)
    if iwarn == 1
        mydisp('error: cannot find ''%s'' in str ''%s''',spattern,str);
    else
        val = [];
        return;
    end
end

if strcmp(stype,'d') || strcmp(stype,'f')
    for ii = 1:numel(a)
        val(ii) = sscanf(a{ii}{1},'%f');
    end
else
    val = a{1};
    if numel(val) == 1
        val = val{1};
    end
end

end

function fun_ex1()
str = 'a2004_b2005.1_cSTR1_dOTHER';
mydisp('input is "%s", four parameters are after a, b, c and d, respectively\n',str);

a = get_var(str,'[ab]*','f'); % two patterns are matched
c = get_var(str,'c*_','s');
d = get_var(str,'_d*','s');
mydisp('The parameters are: %d, %.1f, %s, %s\n',a,c,d);
end

