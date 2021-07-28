function [p,plmplus,plmmin] = iplm(l,m,theRAD,dt)

% IPLM Integrals of the fully normalized associated Legendre functions
% over blocks for a selected order M. 
%
% HOW: p = iplm(l,m,theRAD)		- assumes dt = theRAD(2)-theRAD(1)
%      p = iplm(l,m,theRAD,dt)
%
% IN:
%    l ........ degree (vector). Integer, but not necessarily monotonic.
%               For l < m a vector of zeros will be returned.
%    m ........ order (scalar)
%    theRAD ... co-latitude [rad] (vector)
%    dt ....... integration block-size [rad] (scalar). Default: dt = theRAD(2)-theRAD(1)
%
% OUT: 
%    p ........ Matrix with integrated Legendre functions.
%               Functions are integrated from theRAD(i)-dt/2 till theRAD(i)+dt/2.
%               The matrix has length(TH) rows and length(L) columns, unless L 
%               or TH is scalar. Then the output vector follows the shape of 
%               respectively L or TH. 
% 
% USES:
%    PLM
%
% REMARKS:
%    The blocks at the pole might become too large under circumstances.
%    This is not treated separately, i.e. unwanted output may appear.
%    In case TH is scalar, dt will be 1 (arbitrarily).

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Dimitris TSOULIS (DT), IAPG, TU-Munich  
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    Markus ANTONI (MA), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2012-01-23: MR, input of plm in radian
%    1999-02-??: NS, brush-up (layout, help text, inactive lines,...)
%                    output variable (column extraction)
%                    Plm's BEFORE for-loop -> large speed-up
%                    variable redefinition (e.g. loop variable l)
%    1998-12-??: DT, initial version 
%----------------------------------------------------------------------------
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
% ----------------------------------------------------------------------------

% diagnostics and preliminaries
narginchk(3, 4);

if nargin < 4
   if max(size(theRAD)) == 1
      dt = pi/180;
   else
      dt = theRAD(2) - theRAD(1);
   end
end
if min(size(l)) ~= 1,  error('Degree L must be vector (or scalar)'), end
if any(rem(l,1) ~= 0), error('Vector L contains non-integers.'), end
if max(size(m)) ~= 1,  error('Order M must be scalar.'), end
if rem(m,1) ~=0,       error('Order M must be integer.'), end
if max(size(dt)) ~= 1, error('Block size DT must be scalar.'), end
if dt == 0,            error('DT cannot be zero'); end

% init.
lcol   = size(l, 2); 
trow   = size(theRAD, 1);
n      = length(theRAD);				% number of latitudes
theRAD = theRAD(:);     
if min(theRAD) < 0 || max(theRAD) > pi
    warning('Is the co-latitude ''theta'' given in radian?')
end
lmax = max(l);
mfix = m;				   	% m can be used now as running index.
lvec = l(:)';			   		% l can be used now as running index.
l    = mfix:lmax;

% Initialization of cosine, sine and Plm functions
stplus  = sin(theRAD+dt/2);
stmin   = sin(theRAD-dt/2);
ctplus  = cos(theRAD+dt/2);
ctmin   = cos(theRAD-dt/2);
plmplus = ones(n,lmax+1);
plmmin  = ones(n,lmax+1);
plmplus(:,l+1) = plm(l,mfix,(theRAD+dt/2));	% tesserals
plmmin(:,l+1)  = plm(l,mfix,(theRAD-dt/2));
if mfix > 0
   m   = 1:mfix;
   mm  = 2*m;
   fac = sqrt(2*cumprod((mm+1)./mm));
   [mgr,stp] = meshgrid(m,stplus);
   [fgr,stm] = meshgrid(fac,stmin);
   plmplus(:,m+1)  = fgr.*(stp.^mgr);		% sectorials
   plmmin(:,m+1)   = fgr.*(stm.^mgr);
end


% Initialization of the temporary matrix ptmp, containing the integrated
% functions in its columns, with progressing degree l. The last column of
% ptmp will contain zeros, which is useful for assignments when l < m.
ptmp   = zeros(n,lmax+2);				% L+2 columns !!!
ptmp00 = cos(theRAD - dt/2) - ctplus;
ptmp11 = sqrt(3) / 2 * (dt - ctplus.*stplus + ctmin.*stmin);
ptmp10 = sqrt(3) / 2 * (stplus.^2 - stmin.^2);
ptmp(:,1) = ptmp00;


% Compute first the integrals of order m == 0
if mfix == 0
   
   ptmp(:,2) = ptmp10;
   for l = 2:lmax					% loop over the degree l  
      rootnm  = sqrt( (2*l+1) * (2*l-1) / l^2 );
      root1nm = sqrt( (2*l-1) * (2*l-3) / (l-1)^2 );
      ptmp(:,l+1) = rootnm / (l+1) * ( (l-2) * ptmp(:,l-1) / root1nm + ...
                    stplus.^2 .* plmplus(:,l) - stmin.^2 .* plmmin(:,l) );
   end    
   
else
% Compute the integrals of order m > 0

% First we compute the diagonal element IPmm (lmax == mfix)
   ptmp(:,2) = ptmp11;   
   for l = 2:mfix					% loop over the degree l
      rootmm  = sqrt( (2*l+1) / (2*l) );
      root1mm = sqrt( (2*l-1) / (2*l-2));
      if l == 2
         root1mm = sqrt(3);
      end
      ptmp(:,l+1) = rootmm / (l+1) * ( l * root1mm * ptmp(:,l-1) - ... 
                ( ctplus.*plmplus(:,l+1) - ctmin.*plmmin(:,l+1) ) / rootmm );
   end
      
% the arbitrary element IPlm ( computed only when lmax > mfix)
   
   if lmax > mfix
      
% first we do the element IPlm, for which l - m = 1
      l = mfix + 1;
      rootnm = sqrt( (2*l+1) * (2*l-1) / (l+mfix) / (l-mfix) );
      ptmp(:,l+1) = rootnm / (l+1) * (stplus.^2 .* plmplus(:,l) - ... 
                    stmin.^2 .* plmmin(:,l) );
       
% now we do the rest
      for l = mfix+2:lmax                   % loop over the degree l
         rootnm  = sqrt( (2*l+1) * (2*l-1) / (l+mfix) / (l-mfix) );
         root1nm = sqrt( (2*l-1) * (2*l-3) / (l-1+mfix) / (l-1-mfix) );
         ptmp(:,l+1) = rootnm / (l+1) * ( (l-2) * ptmp(:,l-1) / root1nm + ...
                       stplus.^2 .* plmplus(:,l) - stmin.^2 .* plmmin(:,l) );
      end
      
   end
   
end


% The integrated functions have been computed. What remains to be done, is to
% extract the proper columns from ptmp, corresponding to the vector lvec. 
% If l or theta is scalar the output matrix p reduces to a vector. It should
% have the shape of respectively theta or l in that case.

% p     = zeros(n, length(lvec));		% size declaration.
lind  = find(lvec < mfix);			% index into l < m
pcol  = lvec + 1;				    % index into columns of ptmp
pcol(lind) = (lmax+2) * ones(size(lind));	% Now l < m points to last col.
p     = ptmp(:, pcol);			    % proper column extraction 

if max(size(lvec)) == 1   && min(size(theRAD)) == 1 && (trow == 1), p = p'; end
if max(size(theRAD)) == 1 && min(size(lvec)) == 1   && (lcol == 1), p = p'; end
