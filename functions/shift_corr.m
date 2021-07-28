function [vcorr,x2,opt_shift,max_corr,min_corr]=shift_corr(x,y,max_shift)
if nargin == 0
    fun_ex2;
    return;
end
x = x(:);
y = y(:);
ishift = (-max_shift:max_shift)';
vcorr = zeros(2*max_shift+1,1);
for ii = 1:numel(ishift)
    x2 = shift_series(x,ishift(ii));
    vcorr(ii) = corr(x2,y);
end

% return x2 that gives the max correlation
max_corr = max(vcorr);
min_corr = min(vcorr);
iloc = vcorr==max_corr;
opt_shift = ishift(iloc);

x2 = shift_series(x,ishift(iloc));

end

function x2=shift_series(x,ii)
% ii: positive, right shift
%     negative, left shift
ii = -ii;
x = x(:);
nx = numel(x);
ii = mod(ii,nx);
x3 = repmat(x,2,1);
iloc = ((1+ii):(nx+ii));
x2 = x3(iloc);
end

function fun_ex1()
x = 1:10;
shift_series(x,-1)
end

function fun_ex2()
x = sin(linspace(0,4*pi,200));
y = cos(linspace(0,4*pi,200));
a = load('tmp.mat');
x = a.PC(:,10);
y = a.PC(:,11);
max_shift = 10;
[v,x2,opt_shift] = shift_corr(x,y,max_shift);
xshift = -max_shift:max_shift;
myfigure('f');
subplot(2,1,1)
plot(xshift,v,'o-');
grid on;
title(sprintf('max corr: %.3f; ishift: %d',max(v),opt_shift));

subplot(2,1,2)
plot(x,'r-');
hold on;
plot(x2,'c-');
plot(y,'b-');
hold off;
legend('x',sprintf('x2,ishft=%d',opt_shift),'y');
end
