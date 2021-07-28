function [p, dp, ddp] = plm(l, m, thetaRAD)

% PLM Fully normalized associated Legendre functions for a selected order M
%
% HOW: 
%    p            = plm(l, thetaRAD)			- assumes M=0
%    p            = plm(l, m, thetaRAD)
%    [p, dp]      = plm(l, m, thetaRAD)
%    [p, dp, ddp] = plm(l, m, thetaRAD)
%
% IN:
%    l ........ degree (vector). Integer, but not necessarily monotonic.
%               For l < m a vector of zeros will be returned.
%    m ........ order (scalar). If absent, m = 0 is assumed.
%    thetaRAD . co-latitude [rad] (vector)
%
% OUT:
%    p ........ Matrix with Legendre functions. The matrix has length(thetaRAD) 
%               rows and length(l) columns.
%    dp ....... Matrix with first derivative of Legendre functions. The matrix 
%               has length(thetaRAD) rows and length(l) columns.
%    ddp ...... Matrix with second derivative of Legendre functions. The matrix 
%               has length(thetaRAD) rows and length(l) columns.
%
% SEE ALSO:
%    LEGPOL, YLM, IPLM
%
% REMARKS:
% *  Previous versions calculated the derivatives towards the latitude, 
%    i. e. dP/d\phi are calculated. Please check your code in order to get 
%    the derivatives correctly towards the co-latitude!
%      ->  dP/d\thetaRAD      =   -dP/d\phi
%      ->  d^2P/d\thetaRAD^2  =  d^2P/d\phi^2
% *  For the rotation bundle, negative order m  is considered at the moment 
%    by the relation:
%       P[l,-m] = (-1)^m * P[l,+m]
%    Further changes in normalization might be necessary for the complex
%    versions of spherical harmonics.
%
% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    Matthias WEIGELT (MW), DoGE, UofC
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2015-09-03: MR, remove the input driven autorotating behaviour of output matrices
%    2015-02-10: MA, negative order m for rotation purposes
%    2013-01-29: MA, comments
%    2013-01-23: MA, input argument thetaRAD [deg -> radian]
%    2008-04-04: MW, extension for second derivative
%    2004-11-24: MW, speed up calculation
%    2004-08-13: MW, extension for first derivative
%    1998-07-13: NS, Pmm non-recursive anymore
%    1997-06-09: NS, help text brushed up
%    1994-08-08: NS, initial version
% -------------------------------------------------------------------------
% license:
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the  Free  Software  Foundation; either version 3 of the License, or
%    (at your option) any later version.
%  
%    This  program is distributed in the hope that it will be useful, but 
%    WITHOUT   ANY   WARRANTY;  without  even  the  implied  warranty  of 
%    MERCHANTABILITY  or  FITNESS  FOR  A  PARTICULAR  PURPOSE.  See  the
%    GNU General Public License for more details.
%  
%    You  should  have  received a copy of the GNU General Public License
%    along with Octave; see the file COPYING.  
%    If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------

% Some input checking.
if nargin == 2
   thetaRAD = m;
   m  = 0;
end
if min(size(l)) ~= 1,  error('Degree l must be vector (or scalar)'), end
if any(rem(l,1) ~= 0), error('Vector l contains non-integers.'), end
if max(size(m)) ~= 1,  error('Order m must be scalar.'), end
if rem(m,1)     ~= 0,  error('Order m must be integer.'), end

% negative orders:
msign = 1; if m<0; msign = (-1)^m;end; m= abs(m);


% Preliminaries.
lcol = size(l,2);
trow = size(thetaRAD,1);
lmax = max(l);
if lmax < m, 
    p = zeros(length(thetaRAD),length(l));
    dp = zeros(length(thetaRAD),length(l));
    ddp = zeros(length(thetaRAD),length(l));
    return
end
n    = length(thetaRAD);				% number of latitudes
t    = thetaRAD(:);
x    = cos(t);
y    = sin(t);
lvec = l(:)';					% l can be used now as running index.

if min(t) < -1e-14 || (max(t) - pi) > 1e-14
    warning('Is the co-latitude ''thetaRAD'' given in radian?')
end

% Recursive computation of the temporary matrix ptmp, containing the Legendre
% functions in its columns, with progressing degree l. The last column of
% ptmp will contain zeros, which is useful for assignments when l < m.
ptmp  = zeros(n,lmax-m+2);
if nargout >= 2,                  % first derivative needs also P_{n,m+1} and P_{n,m-1}
    ptmp_m1 = zeros(n,lmax-m+3); 
    ptmp_p1 = zeros(n,lmax-m+1); 
    dptmp   = zeros(n,lmax-m+2); 
end
if nargout == 3,                  % second derivative needs also dP_{n,m+1} and dP_{n,m-1}
    dptmp_m1 = zeros(n,lmax-m+3); 
    dptmp_p1 = zeros(n,lmax-m+1); 
    ptmp_m2  = zeros(n,lmax-m+4); % but these first derivative need dP_{n,m+2} and dP_{n,m-2}
    ptmp_p2  = zeros(n,lmax-m); 
    ddptmp   = zeros(n,lmax-m+2); 
end

%--------------------------------------------------------------------
% sectorial recursion: PM (non-recursive, though)
%--------------------------------------------------------------------
ptmp(:,1) = secrecur(m,y);
if nargout >= 2,  % frist derivative needs preceding and subsequent element
    if m > 0,    ptmp_m1(:,1) = secrecur(m-1,y); end     % preceding elements
    if m < lmax, ptmp_p1(:,1) = secrecur(m+1,y); end     % subsequent elemtens
end
if nargout == 3,  % second derivative needs P_{n,m+2} and P_{n,m-2} as well
    if m > 1,      ptmp_m2(:,1) = secrecur(m-2,y); end   % preceding elements
    if m < lmax-1, ptmp_p2(:,1) = secrecur(m+2,y); end   % subsequent elemtens
end

%--------------------------------------------------------------------
% l-recursion: P
%--------------------------------------------------------------------
ptmp = lrecur(ptmp,x,m,lmax);
if nargout >= 2,  % frist derivative needs preceding and subsequent element
    if m > 0,    ptmp_m1 = lrecur(ptmp_m1,x,m-1,lmax); end    % preceding elements
    if m < lmax, ptmp_p1 = lrecur(ptmp_p1,x,m+1,lmax); end    % subsequent elemtens
end
if nargout == 3,  % second derivative needs P_{n,m+2} and P_{n,m-2} as well
    if m > 1,      ptmp_m2 = lrecur(ptmp_m2,x,m-2,lmax); end  % preceding elements
    if m < lmax-1, ptmp_p2 = lrecur(ptmp_p2,x,m+2,lmax); end  % subsequent elemtens
end

%--------------------------------------------------------------------
% now compute the derivatives 
%--------------------------------------------------------------------
if nargout >= 2,                                   % first derivative
    dptmp = derivALF(dptmp,ptmp_m1,ptmp_p1,m,lmax);
end
if nargout == 3,                                   % second derivative
    if m > 0,    dptmp_m1 = derivALF(dptmp_m1,ptmp_m2,ptmp,m-1,lmax); end
    if m < lmax, dptmp_p1 = derivALF(dptmp_p1,ptmp,ptmp_p2,m+1,lmax); end
    ddptmp = derivALF(ddptmp,dptmp_m1,dptmp_p1,m,lmax);
end

%--------------------------------------------------------------------
% The Legendre functions have been computed. What remains to be done, is to
% extract the proper columns from ptmp, corresponding to the vector lvec. 
% If l or thetaRAD is scalar the output matrix p reduces to a vector. It should
% have the shape of respectively thetaRAD or l in that case.
%--------------------------------------------------------------------
lind       = (lvec < m);			        % index into l < m
pcol       = lvec - m + 1;			        % index into columns of ptmp
pcol(lind) = (lmax-m+2)*ones(sum(lind),1);	% Now l < m points to last col.
p          = msign * ptmp(:,pcol);			% proper column extraction + handeling of negative order
if nargout >= 2, 
    dp =  msign * dptmp(:,pcol);            % proper column extraction + handeling of negative order
end   
if nargout == 3, 
    ddp = msign * ddptmp(:,pcol);           % proper column extraction + handeling of negative order
end   


% *************************************************************************
%                            INLINE FUNCTIONS
% *************************************************************************
% function for the sectorial recursion, non-recursive though
function out = secrecur(m,y)
if m == 0
   fac = 1;
else
   mm  = 2*(1:m);
   fac = sqrt(2*prod((mm+1)./mm));
end
out = fac*y.^m;             % The 1st column of ptmp.
% -------------------------------------------------------------------------

% function for the l-recursion
function in = lrecur(in,x,m,lmax)
for l = m+1:lmax
   col   = l - m + 1;			% points to the next column of ptmp
   root1 = realsqrt( (2*l+1)*(2*l-1)/((l-m)*(l+m)) ) ;
   root2 = realsqrt( (2*l+1)*(l+m-1)*(l-m-1) / ( (2*l-3)*(l-m)*(l+m) ) );

   % recursion
   if l == m+1
       in(:,col) = root1 *x.*in(:,col-1);
   else
       in(:,col) = root1 *x.*in(:,col-1) - root2 *in(:,col-2);
   end
end
% -------------------------------------------------------------------------

% function to calculate the derivate
function in = derivALF(in,miin,plin,m,lmax)

l = (m:lmax+1);
if m == 0
    in(:,1) = 0;
    if lmax > m, in(:,2:end) = bsxfun(@times,plin,-realsqrt((l(2:end)+1).*l(2:end)./2)); end %   (-ones(n,1)*realsqrt((l(2:end)+1).*l(2:end)./2)).*plin; end
elseif m == 1
    in(:,1) = miin(:,2);
    if lmax > m, in(:,2:end) =  bsxfun(@times,miin(:,3:end),realsqrt((l(2:end)+1).*l(2:end)./2)) - ...
            0.5.*bsxfun(@times,plin,realsqrt((l(2:end)-1).*(l(2:end)+2))); end
elseif m == lmax
    in(:,1) = realsqrt(m./2).*miin(:,2);
else
    in(:,1) = realsqrt(m./2).*miin(:,2);
    if lmax > m, in(:,2:end) = 0.5.*bsxfun(@times,miin(:,3:end),realsqrt((l(2:end)+m).*(l(2:end)-m+1))) - ...
            0.5.*bsxfun(@times,plin,realsqrt((l(2:end)-m).*(l(2:end)+m+1))); end    
end







