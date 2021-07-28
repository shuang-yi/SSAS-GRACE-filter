function m2=mlay(m1)
%{
force a matrix to be in the flat form (d1 <= d2)
%}
% {{matrix; transpose; lay}}

if nargin == 0  
    fun_ex1;  
    return;
end
% ----- end of head -----
N = ndims(m1);
if N ==1 
    m2 = m1;
elseif N == 2
    [i1,i2] = size(m1);
    if i1>i2
        m2 = m1.';
    else
        m2 = m1;
    end
else
    error('the input should only be 2D\n');
end

end
%%
function fun_ex1
mlay(1)
mlay(1:3)
mlay(11:14)
mlay(rand(2,4))
mlay(rand(4,2))
end