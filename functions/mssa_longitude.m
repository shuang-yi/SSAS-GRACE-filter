function [LAMBDA, RC, EOF2, PC, ilist] = mssa_longitude(X, PAR)
% mypath2('share','Tools\SSA',1); % mssa_yi

% {{MSSA; LLZ; longitude}}

if nargin == 0
     fun_ex1;
     return;
end
% PAR0 = struct('M',90,'K',360,'only_pair',0,'corr',0.9,'ilat_list',[]);
% PAR = var_initial(PAR0, varargin);

M = PAR.M;
K = PAR.K;

% N = numel(t); % number of time epoches
% L = N-M+1;
[N,Nx] = size(X); % number of channels

K_try = M*Nx;
if K>K_try
    mydisp('K = %d, but should < %d\n',K,K_try);
    K = K_try;
end

if ~isfield(PAR,'ilat_list') || isempty(PAR.ilat_list)
    ilat_list = 1:Nx;
else
    ilat_list = PAR.ilat_list;
end
    
Y=create_Y(X(:,ilat_list), M);

% tic
[U,S,V] = svd(Y,'econ'); % Y = U*S*V';
% toc

% V2 = reshape(V,M,Nx,N);
% weight = max(abs(V2),[],1);
    
% [U,S,V] = svds(Y,K);
LAMBDA = diag(S.^2);
PC = U;
% EOF = V*S;

Y_all=create_Y(X, M);
EOF = Y_all'*U;

PC = PC(:,1:K);
EOF = EOF(:,1:K);
% EOF2 = zeros(N,M,Nx);
EOF2 = reshape(EOF,M,Nx,K);
for ii = 1:Nx
    EOF3(:,:,ii) = EOF2(:,ii,:);
end
EOF2 = EOF3;

%% Calculate reconstructed components RC
RC = zeros(N,K,Nx);

%     EOF2(:,1:K,ix) = EOF(ik+1:ik+M, :);

if PAR.only_pair == 1 % only return paired RC
    [ilist]=paired_PC(PC,PAR.corr,PAR.max_shift);
%     fprintf('Only_pair: %.2f of paired RC will be returned\n',numel(ilist)/K);
else % return all RCs
    ilist = 1:K;
end

% Y=create_Y(X, M, PAR.ind_nan);
% ind_nan = isnan(Y);

for m = ilist(:)' 
%     fprintf('m=%d\n',m);
    buf0 = PC(:,m) * EOF(:, m)'; % invert projection
%     buf0(ind_nan) = NaN;
    for ix = 1:Nx
        ik = (ix-1)*M;
        buf1 = buf0(:,ik+1:ik+M);
        
        buf2 = fun_shift_matrix(buf1);
        RC(:,m,ix) = mean_with_nan(buf2);
    end
end

end

function Y=create_Y(X, M, ind_nan)
if nargin < 3
    ind_nan = zeros(size(X))==1;
end
if numel(ind_nan) == 1 && (isnan(ind_nan) || isempty(ind_nan))
    ind_nan = zeros(size(X))==1;
end

X(ind_nan) = NaN;

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

function mat2 = fun_shift_matrix(mat)
[nx,ny] = size(mat);
mat2 = zeros(nx,ny);
for ij = 1:ny
    iloc = [nx-ij+2:nx,1:nx-ij+1];
    mat2(:,ij) = mat(iloc,ij);
end
end

function buf_mean=mean_with_nan(buf)
if sum(isnan(buf(:))) == 0
    buf_mean = mean(buf,2);
else
    buf_mean = zeros(size(buf,1),1);
    for ii = 1:size(buf,1)
        ind1 = ~isnan(buf(ii,:));
        buf_mean(ii) = mean(buf(ii,ind1));
    end
end
end

%%
function fun_ex1()

% Set general Parameters
M = 50;    % window length = embedding dimension
N = 200;   % length of generated time series
T = 20;    % period length of sine function
stdnoise = 1; % noise-to-signal ratio


% Create time series X
% First of all, we generate a time series, a sine function of length N with
% observational white noise

rng('default');

t = (1:N)';
X0 = sin(2*pi*t/T) + t*0.0;% + t.^2*0.0005;%+sin(2*pi*t/2/T);
noise = stdnoise*randn(size(X0));
X = X0 + noise;

K =10;
[lambda, RC ] = mssa_longitude(X, struct('M',M, 'K',K, 'only_pair',1,'corr',0.9,'max_shift',3));

figure;
subplot(2,1,1)
plot(lambda(1:min([20,numel(lambda)])),'o-');
title('eigenvalues LAMBDA');
subplot(2,1,2)
plot(t,X,'b-',t,sum(RC(:,1:2),2),'r-');
hold on;
plot(t,X0,'k--');
legend('Original','Reconstruction with RCs 1-2','true value');
end

function fun_ex2()
M = 30;    % window length of SSA
N = 200;   % length of generated time series
T = 22;    % period length of sine function
stdnoise = 1; % noise-to-signal ratio

rng('default');
%% Create time series X
% First of all, we generate two time series, a sine function
% of length N and the same function to the power of 3, both
% with observational white noise.

t = (1:N)';
X1 = sin(2*pi*t/T);     % sine function
X2 = cos(2*pi*t/T).^3;  % nonlinear transformation
X3 = sin(2*pi*t/T/2)+t/100;
X0 = [X1, X2, X3];
noise = stdnoise*randn(size(X0));  % Gaussian noise

X = X0+noise;

Nx = size(X,2);

K = 7;
[LAMBDA, RC, EOF, PC] = mssa_longitude(t,X,M,K);
myfigure(8,6)
subplot(2,3,[1,4])
plot(LAMBDA,'o-');
title('lambda');
subplot(2,3,[2,3]);
plot(EOF(:,1:K),'-','linew',2);
title('EOF');
subplot(2,3,[5,6]);
plot(PC(:,1:K),'-','linew',2);
title('PC');
myfigure(8,6);
isel = [1,2,3,4,6,7];
for ii = 1:Nx
    subplot(Nx,1,ii);
    plot(t,X(:,ii),t,sum(RC(:,isel,ii),2));
    hold on;
    plot(t,X0(:,ii),'k');
    hold off;
    if ii == 1
        legend('observation with noise',sprintf('MSSA, M=%d, K=%d',M,K),'True model');
    end
    title(sprintf('series %d',ii));
end

end