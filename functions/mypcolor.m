function g1 = mypcolor(x,y,z)
% {{plot}} {{2D}}
%
% mypcolor( x, y, z)
% mypcolor( z )

x = squeeze(x);
[n,m] = size(x);

ilabel = 0;
if nargin == 1
   z = x;
   [x,y] = meshgrid(1:m,1:n);
   ilabel = 1;
elseif nargin == 3
    y = squeeze(y);
    z = squeeze(z);
end

if n == 1
    dy = 1;
else
    dy = y(2,1) - y(1,1);
end

if m == 1
    dx = 1;
else
    dx = x(1,2) - x(1,1);
end

if dx==0 && dy==0
    x = x';
    y = y';
    z = z';
    dx = x(1,2) - x(1,1);
    dy = y(2,1) - y(1,1);
    [n,m] = size(x);
end
% x0 = min(x(1,:));
% x1 = max(x(1,:));
% y0 = min(y(:,1));
% y1 = max(y(:,1));
% xx = ( x0-dx/2 ):dx:(x1+dx/2);
% yy = ( y0-dy/2 ):dy:(y1+dy/2);
xx = linspace( x(1,1)-dx/2, x(1,m)+dx/2, m+1 );
yy = linspace( y(1,1)-dy/2, y(n,1)+dy/2, n+1 );
rtmp =  mean(z(:));
z(n+1,:) = rtmp;
z(:,m+1) = rtmp;
g1 = pcolor(xx,yy,z*1);
set(g1,'linestyle','none');
colorbar;
colormap jet
if ilabel == 1 % the location in a matrix
    set(gca,'ydir','reverse');
    xlabel('y');ylabel('x');
end
set(gca,'layer','top')
end