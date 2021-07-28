function [LAMBDA, RC, EOF2, PC, ilist] = mssa_longitude_missing(X, PAR)

% {{MSSA; LLZ; longitude}}


% PAR0 = struct('M',90,'K',360,'only_pair',0,'corr',0.9,'ilat_list',[]);
% PAR = var_initial(PAR0, varargin);

M = PAR.M;
K = PAR.K;

% N = numel(t); % number of time epoches
% L = N-M+1;
[N,Nx] = size(X); % number of channels

Nlist = 1:Nx;

K_try = M*Nx;
if K>K_try
    mydisp('K = %d, but should < %d\n',K,K_try);
    K = K_try;
end

ind_nan = PAR.ind_nan;
ind_nan_col = sum(ind_nan,1)>0;

ithreshold = 1D-2;

if sum(ind_nan(:)) > 0
    N_iter = 10000; % max number of iteration  
    iK = 1; 
else
    N_iter = 1;
    iK = K;   
end

X0 = X;
X1 = zeros(size(X));
X_save = zeros(size(X));

X(ind_nan) = 0;

isave = 1;
iKsave = 0;
iend = 0;
for iter = 1:N_iter
    if iend == 1
        break;
    end
%     tic
    v_diff = rms(X(ind_nan) - X1(ind_nan))/rms(X(ind_nan));
    X(ind_nan) = X1(ind_nan);
%     ttt = toc;
%     fprintf('time1=%f\n',ttt);
    if v_diff < ithreshold
        if mod(iK,2) == 0
            [ipaired]=paired_PC(PC(:,1:2),0.9,round(N/10));
            if isempty(ipaired)
                X_save = X_save+squeeze(RC(:,1,:));
                X = X-squeeze(RC(:,1,:));
                isave = isave+1;
                iKsave = iKsave+1;
            else
                X_save = X_save+squeeze(sum(RC(:,1:2,:),2));
                X = X-squeeze(sum(RC(:,1:2,:),2));
                isave = isave+2;
                iKsave = iKsave+2;
            end
            iK = 1;
            if iKsave >= K
                X = X0;
                X(ind_nan) = X_save(ind_nan);
                iK = iKsave;
                iend = 1;
            end
        else
            iK = iK+1;
%             fprintf('ii=%d, ik=%d\n',iter,iK);
        end
    end
    
%     ttt = toc;
%     fprintf('time2=%f\n',ttt);
    
    %%
    Y=create_Y(X, M);
    
    % tic
    [U,S,V] = svd(Y,'econ'); % Y = U*S*V';
    % toc
    
    % V2 = reshape(V,M,Nx,N);
    % weight = max(abs(V2),[],1);
    
    % [U,S,V] = svds(Y,K);
    
    %     PC = U;
    %     EOF = V*S;
    
    PC = U(:,1:iK);
    EOF = V(:,1:iK)*S(1:iK,1:iK);
    
%     ttt = toc;
%     fprintf('time3=%f\n',ttt);

    
    %% Calculate reconstructed components RC
    
    ilist = 1:iK;
    
    RC = zeros(N,numel(ilist),Nx);
    for m = ilist(:)'
        %     fprintf('m=%d\n',m);
        buf0 = PC(:,m) * EOF(:, m)'; % invert projection
        %     buf0(ind_nan) = NaN;
        if iend == 1
            ixlist = 1:Nx;
        else % only return these in NaN
            ixlist = Nlist(ind_nan_col);
        end
        
        for ix = ixlist
            ij = (ix-1)*M;
            buf1 = buf0(:,ij+1:ij+M);
            
            buf2 = fun_shift_matrix(buf1);
            %  RC(:,m,ix) = mean_with_nan(buf2);
            RC(:,m,ix) = mean(buf2,2);
        end
        
    end
    
%     ttt = toc;
%     fprintf('time4=%f\n',ttt);
    
    X1 = squeeze(sum(RC,2));
end

%%
LAMBDA = diag(S.^2);

EOF2 = reshape(EOF,M,Nx,iK);
for ii = 1:Nx
    EOF3(:,:,ii) = EOF2(:,ii,:);
end
EOF2 = EOF3;

if PAR.only_pair == 1 % only return paired RC
    [ilist]=paired_PC(PC,PAR.corr,PAR.max_shift);
end

end

function Y = create_Y(X, M, ind_nan)
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
