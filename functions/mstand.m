function m2=mstand(m1)
%{
force a matrix to be in the standing form (d1 >= d2)
%}
% {{matrix; transpose; stand}}

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
    if i1<i2
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
mstand(1)
mstand(1:3)
mstand(11:14)
mstand(rand(2,4))
mstand(rand(4,2))
end