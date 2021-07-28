function Y=matrix_lag(X,MM)
% {{SSA; lag; matrix}}
% X: time series
% MM: window width

if nargin == 0
    fun_ex1;
    return;
end

X = X(:);
N = numel(X);
L = N-MM+1;


Y = zeros(L,MM);
for m = 1:MM
    Y(:,m) = X((1:L) + m-1);    
end

end

function fun_ex1()
X = 1:12;
matrix_lag(X,3)
end