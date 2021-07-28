function cs = gsha(f, method, grd, lmax)

% GSHA global spherical harmonic analysis
%
% IN:
%    f ....... global field of size (lmax+1)*2*lmax or lmax*2*lmax
%    method .. string argument, defining the analysis method:
%              - 'ls' ..... least squares
%              - 'wls' .... weighted least squares 
%              - 'aq' ..... approximate quadrature 
%              - 'fnm' .... first neumann method
%              - 'snm' .... second neumann method
%              - 'mean' ... block mean values (use of integrated Plm)
%    grd .... optional string argument, defining the grid:
%              - 'pole', 'mesh' ...... (default if lmax+1), equi-angular (lmax+1)*2*lmax, 
%                                      including poles and Greenwich meridian.
%              - 'block', 'cell' ..... (default if lmax), equi-angular block midpoints lmax*2lmax
%              - 'neumann', 'gauss' .. Gauss-Neumann grid (lmax+1)*2*lmax
%    lmax .... maximum degree of development
%
% OUT:
%    cs ...... Clm, Slm in |C\S| format
%
% USES:
%    plm, iplm, neumann, sc2cs 
%
% SEE ALSO:
%    GSHS
%
% REMARKS:
%    TBD - Zlm-functions option
%        - eigengrav, GRS80
%        - When 'pole' grid, m = 1 yields singular Plm-matrix!

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Dimitris TSOULIS (DT), IAPG, TU-Munich
%    Nico SNEEUW (NS), IAGP, TU Munich
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    Matthias ROTH (MR), GI, Uni Stuttgart 
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2016-01-22: MR, change waitbar --> twaitbar
%    2015-07-06: MR, remove typo "ischar", thanks to Uli Mayer for pointing 
%                    them out
%    2014-01-15: MR, revise help text, beautify code
%    2012-01-23: MA, input of plm in radian
%    1999-02-01: NS, brush up (help text, layout, removal of unused commands, ...)
%                    restructuring of 'mean' method over 1st and 2nd analysis
%                    'mean' method as quadrature instead of LS in 2nd step 
%    1998-11-??: DT, initial version
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

%--------------------------------------------------------------------------
% diagnostics and preliminaries
%--------------------------------------------------------------------------
narginchk(2, 4); 

%--------------------------------------------------------------------------
% Grid definition
%--------------------------------------------------------------------------
[rows, cols] = size(f);

if cols == 2 * rows                     % 'block' | 'cell' grid 
   if nargin < 4, lmax = rows; end      % default
   if nargin < 3, grd = 'block'; end	% default
   if ~strcmp(grd, 'block') && ~strcmp(grd, 'cell') 
      error('Your GRID variable should be either block or cell')
   end
   
   n     = rows;
   dt    = 180 / n;
   theta = (dt/2:dt:180)';
   lam   = (dt/2:dt:360);               % dt = dlam

elseif cols == 2 * rows - 2             % 'pole' | 'mesh' | 'neumann'
   
   if nargin < 4, lmax = rows-1; end	% default
   if nargin < 3, grd = 'pole'; end    % default
  
   n     = rows - 1;
   dt    = 180 / n;
   if strcmp(grd,'pole') || strcmp(grd,'mesh') 
      theta = (0:dt:180)';
      lam   = (0:dt:360-dt);
   elseif strcmp(grd,'neumann') || strcmp(grd,'gauss') 
      [gw,gx] = neumann(n+1);
      theta = acos(flipud(gx))*180/pi;
      lam   = (0:dt:360-dt);
   else
      error('The wrong type of GRID')
   end
   
else
   error('Invalid size of matrix F')
end	

theRAD = theta * pi/180;
if min(theRAD)<0 || max(theRAD)>pi
    warning('Is the co-latitude ''thetaRAD'' given in radian?')
end

%--------------------------------------------------------------------------
% further diagnostics
%--------------------------------------------------------------------------
if ~ischar(grd) || ~ischar(method) % deprecated: ~isstr(grd) || ~isstr(method)
    error('GRID and METHOD must be strings')
end
if strcmp(method, 'snm') && ~(strcmp(grd, 'neumann') || strcmp(grd, 'gauss'))
	error('2nd Neumann method ONLY on a ''neumann''/''gauss'' GRID')
end
if strcmp(method, 'mean') && ~(strcmp(grd, 'block') || strcmp(grd, 'cell'))
    error('Block mean method ONLY on a ''block''/''cell'' GRID')
end
if lmax > n
    error('Maximum degree of development is higher than number of rows of input.')
end


%--------------------------------------------------------------------------
% Init.
%--------------------------------------------------------------------------
L   = n;
% a   = zeros(rows, L+1);
% b   = zeros(rows, L+1);
clm = zeros(L+1, L+1);
slm = zeros(L+1, L+1);

%--------------------------------------------------------------------------
% 1st step analysis: Am(theta) & Bm(theta)
%--------------------------------------------------------------------------
m = 0:L;
c = cos(lam' * m * pi/180);
s = sin(lam' * m * pi/180);

% preserving the orthogonality (except for 'mean' case)
% we distinguish between 'block' and 'pole' type grids (in lambda)

if strcmp(grd,'block') || strcmp(grd,'cell')
    if strcmp(method, 'mean')
        dl = dt;                    % longitude block size
        c(:, 1) = dl/360 * ones(2 * n, 1);	% ICm for m = 0, including 1/(1+dm0)/pi
        m       = 1:L;
        ms      = 2 ./ m .* sin(m * dl/2 * pi/180) / pi;   
        c(:, 2:L+1) = c(:,2:L+1) .* ms(ones(2*n,1),:);	% ICm
        s(:, 2:L+1) = s(:,2:L+1) .* ms(ones(2*n,1),:);	% ISm
    else
        c = c/L; s = s/L;
        c(:, 1)   = c(:,1)/2;       % / (1 + dm0 - dmL)
        s(:, L+1) = s(:,L+1)/2;     % / (1 - dm0 + dmL)
        c(:, L+1) = zeros(2*n, 1);	% CLL unestimable
        s(:, 1)   = zeros(2*n, 1);	% Sl0    -"-
    end
else                                % 'pole' | 'mesh' | 'neumann'
    c = c/L; s = s/L;
    c(:,[1 L+1]) = c(:,[1 L+1])/2;	% / (1 + dm0 + dmL)
    s(:,[1 L+1]) = zeros(2*n,2);	% / (1 - dm0 - dmL), Sl0 & SLL unestimable                 
end

a = f * c;
b = f * s;

%--------------------------------------------------------------------------
% 2nd step analysis: Clm & Slm
%--------------------------------------------------------------------------
% hwb = twaitbar('init', [], 'GSHA Percentage of orders m ready ');  % initialize the waitbar

switch lower(method)
  case 'ls'            % Least squares solution
    for m = 0:L
        p  = plm(m:L, m, theRAD);
        ai = a(:, m+1);
        bi = b(:, m+1);
        clm(m+1:L+1, m+1) = p \ ai;
        slm(m+1:L+1, m+1) = p \ bi;  
%         hwb = twaitbar((m + 1) / (L + 1), hwb); % update the waitbar
    end
  case 'wls'        % Weighted Least Squares
    si = sin(theRAD);
    si = 2 * si / sum(si);
    for m = 0:L
        p   = plm(m:L, m, theRAD);
        ai  = a(:, m+1);
        bi  = b(:, m+1);
        d   = 1:length(theRAD);
        pts = p' * sparse(d,d,si);
        clm(m+1:L+1, m+1) = (pts * p) \ pts * ai;
        slm(m+1:L+1, m+1) = (pts * p) \ pts * bi;  
%         hwb = twaitbar((m + 1) / (L + 1), hwb); % update the waitbar
    end
  case 'aq'         % Approximate Quadrature
    si = sin(theRAD);
    si = 2 * si / sum(si);
    for m = 0:L
        p  = plm(m:L, m, theRAD);
        ai = a(:, m+1);
        bi = b(:, m+1);
        clm(m+1:L+1, m+1) = (1 + (m==0))/4 * p' * (si.*ai);
        slm(m+1:L+1, m+1) = (1 + (m==0))/4 * p' * (si.*bi);   
%         hwb = twaitbar((m + 1) / (L + 1), hwb); % update the waitbar
    end
  case 'fnm'        % 1st Neumann method (exact up to L/2)
    w = neumann(cos(theRAD));
    for m = 0:L
        p  = plm(m:L, m, theRAD);
        ai = a(:, m+1);
        bi = b(:, m+1);
        clm(m+1:L+1, m+1) = (1 + (m==0))/4 * p' * (w.*ai);
        slm(m+1:L+1, m+1) = (1 + (m==0))/4 * p' * (w.*bi);
%         hwb = twaitbar((m + 1) / (L + 1), hwb); % update the waitbar
    end
  case 'snm'        % 2nd Neumann method (exact)
    % weigths determined above already
    for m = 0:L
        p  = plm(m:L, m, theRAD);
        ai = a(:, m+1);
        bi = b(:, m+1);
        clm(m+1:L+1, m+1) = (1 + (m==0))/4 * p' * (gw.*ai);
        slm(m+1:L+1, m+1) = (1 + (m==0))/4 * p' * (gw.*bi);
%         hwb = twaitbar((m + 1) / (L + 1), hwb); % update the waitbar
    end
  case 'mean'       % block mean values
    for m = 0:L
        p  = iplm(m:L, m, theRAD); % integrated Legendre
        ai = a(:, m+1);
        bi = b(:, m+1);
        clm(m+1:L+1, m+1) = (1 + (m==0))/4 * p' * ai;
        slm(m+1:L+1, m+1) = (1 + (m==0))/4 * p' * bi;   
%         hwb = twaitbar((m + 1) / (L + 1), hwb); % update the waitbar
    end
  otherwise
    error('Choose a valid METHOD')
end

% twaitbar('close', hwb); % finalize the waitbar

%--------------------------------------------------------------------------
% Write the coefficients Clm & Slm in |C\S| format
%--------------------------------------------------------------------------
slm = fliplr(slm);
cs  = sc2cs([slm(:, 1:L) clm]);
cs  = cs(1:lmax+1, 1:lmax+1);
