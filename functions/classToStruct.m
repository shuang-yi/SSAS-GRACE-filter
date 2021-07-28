function s = classToStruct(obj)
% {{class; object; convert; struct}}
if nargin == 0
    fun_ex1;
    return;
end
% ----- end of head -----

if ~isobject(obj)
    error('input should be an object\n');
end
props = properties(obj);
for ii = 1:numel(obj)
    for p = 1:numel(props)
        if isobject(obj(ii).(props{p})) % nested object value
            s(ii).(props{p})=classToStruct(obj(ii).(props{p}));
        else
            s(ii).(props{p})=obj(ii).(props{p});
        end
    end
end
end

%%
function fun_ex1
a = cSH(20,1);
a(2) = cSH(20,2);
a = a.setinfo('name','data1','history','a test in classToStruct');
b = classToStruct(a);
end