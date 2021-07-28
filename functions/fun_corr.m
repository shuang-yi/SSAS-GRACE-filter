function [vcorr,Ctoep]=fun_corr(x, maxlag, varargin)
%
%
% See also

if nargin == 0
    help(mfilename);
    fun_ex1;
    return;
end

PAR0 = struct('ex',0,'type',1);
PAR = var_initial(PAR0, varargin);

% --------------------end of head-----------

if PAR.type == 3 % covariance, mannual
    Y = matrix_lag(x,maxlag+1);
    Ctoep = Y'*Y;
    vcorr = NaN;
    return;
end

x = x(:);
Np = numel(x);
i0 = 1;
i1 = Np;
for ii = 0:maxlag
    x1 = x(i0:i1-ii);
    x2 = x(i0+ii:i1);
    ind = ~isnan(x1+x2);
    if sum(ind)>1
        if PAR.type == 1 % covariance
            vcorr(ii+1) = sum(x1(ind).*x2(ind))/sum(ind);
        elseif PAR.type == 2 % correlation
            vcorr(ii+1) = corr(x1(ind),x2(ind));
        end
    else
        vcorr(ii+1) = NaN;
    end
end



Ctoep=toeplitz(vcorr);

end

%% example

function fun_ex1()
x = sin(linspace(0,pi*2,100));
x = x+randn(size(x))*1;
maxlag = numel(x)/2;

xx = xcorr(x,maxlag,'unbiased');
result1=toeplitz(xx(maxlag+1:end));

x(4:30) = NaN;
[vcorr,Ctoep]=fun_corr(x, maxlag);
% result2(1) = sum(x.^2);

myfigure('f','name','fun_corr');
subplot(1,2,1)
imagesc(result1); colorbar;
title('toeplitz');
subplot(1,2,2)
imagesc(Ctoep); colorbar;
title('corr');
colormap jet
end

function fun_ex2()
x = sin(linspace(0,pi*2,100));
x = x+randn(size(x))*1;
maxlag = 10;

result1 = xcorr(x,maxlag,'unbiased');

result2=fun_corr(x, maxlag);

end