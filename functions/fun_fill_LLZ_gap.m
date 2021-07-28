function LLZ=fun_fill_LLZ_gap(LLZ,ran_or_ind_gap)
% mypath2('mypaper','P_G7_SSA/SSA-gap-filling/SSA_YI/GRACE')
%{

%}
% {{ssa; gap-filling; LLZ}}

KK =30;
M = 90;

% ----- end of head -----

if numel(ran_or_ind_gap) == 4 % input is ran
    rr = ran_or_ind_gap;
    ind_nan = LLZ.lon >=rr(1) & LLZ.lon<=rr(2) & LLZ.lat>=rr(3) & LLZ.lat<rr(4);
else % input is ind
    ind_nan = ran_or_ind_gap;
end

ind_nan = ind_nan';

ind_lat = sum(ind_nan,1)>0;

X = LLZ.rg';
X2 = X;
X2(ind_nan) = NaN;
X_nan_lat = X(:,ind_lat);

Y_nan = create_Y(X_nan_lat,M);

[u,~,~] =svds(Y_nan,KK);

EOF = u(:,1:KK);


X_todo = X2(:,ind_lat);
Ntodo = size(X_todo,2);
X_F = zeros(size(X_todo));
for ii = 1:Ntodo
    [X_F(:,ii)] = ssa_missing_m3_longitude(X_todo(:,ii), EOF, M, KK);
%     figure; plot(X_nan_lat(:,ii)); hold on; plot(X_F(:,ii));
end

X(:,ind_lat) = X_F;
LLZ.rg = X';

end

function Y = create_Y(X, M)

Nx = size(X,2); % number of channels

X2 = [X;X]; % repeat the matrix to recurrently sample the observations.

% X = X - mean(X);            % remove mean value
% X = X/std(X,1);             % normalize to standard deviation 1

N = size(X,1);
Y=zeros(N,M*Nx);

for ix = 1:Nx
    ik = (ix-1)*M;
    for m=1:M
        Y(:,m+ik) = X2((1:N)+m-1,ix);
    end
end
end

function [X2] = ssa_missing_m3_longitude(X, EOF, MM, KK)
% {{SSA; series; longitude}}

% X: time series
% MM: window length
% KK: number of RC

% RC(:,1:K)


if numel(KK)>1
    error('K should be a single value');
end

if numel(MM) > 1
    error('M should be a single value');
end

if KK > MM
    error('K = %d should < cons_M = %d\n',KK,MM);
end

X = mstand(X);
% XX0 = X;
N = size(X,1);

ind_nan = isnan(X);
% L = N - MM+1;
% L = N;


Q = zeros(N);
for m = 1:KK
    V = zeros(MM,N);
    for ii = 1:MM
%         V(ii,ii:(ii+MM-1)) = EOF(:,m);
        V(ii,:) = circshift(EOF(:,m)',[0,ii-1]);
    end
    %     PC2(:,m) = V*X;
    %     RC2(:,m) = (W*V'*V)*X;
    
    Q = Q+(V'*V)/MM;
end

X2 = X;
X2(ind_nan) = 0;

P = eye(N) - Q;
d = -P*X2;
G = P(:,ind_nan);
% G*x= d;
X0 = (G'*G)\G'*d;
X2(ind_nan) = X0;


end