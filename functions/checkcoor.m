function [lamRAD ,phiRAD, r] = checkcoor(lamRAD, phiRAD, r, ae, checkdim)

% A small script for checking spherical coordiantes (longitude/latiude/r)
% Are the vectors of the same dimension?
% Are the angles given in degree?
% Fill a scalar radius by the same elements?
%
% The script is executed within the SH synthesis routines like 
% gradpshs.m, gshsptw.m, nablaPot.m, pshs.m and tensphs.m

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%   2013-02-13: MR, change function names, brush up comments
%   2013-02-06: MA, script -> function
%   2013-01-29: MA, comments
%   2013-01-23: MA, initial version
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

if isempty(r)
    r = ae;
end

% check if r, lam and phi are vectors
if ~isvector(lamRAD), error('longitude ''lam'' must be scalar or vectorial.'); end
if ~isvector(phiRAD), error('latitude ''phi'' must be scalar or vectorial.');  end
if ~isvector(r),      error('radius ''r'' must be scalar or vectorial.');      end

% standing vectors
r    = r(:);
phiRAD = phiRAD(:);
lamRAD  = lamRAD(:);

% check units
if max(abs(phiRAD)) > pi/2
    warning('Is the latitude ''phi'' given in radian?')
end
if max(abs(lamRAD)) > 2*pi
    warning('Is the longitude ''lam'' given in radian?')
end
if min(r) < ae
    warning('radius ''r'' is smaller than the Earth radius')
end




switch checkdim
case 'pointwise'
    
    
    if isscalar(r) == 1
        r = r*ones(size(phiRAD));
    end

    if numel(phiRAD) ~= numel(lamRAD),    error('''phi'' and ''lam'' must have the same length.'); end
    if numel(phiRAD) ~= numel(r),         error('''phi'',''r'' and ''lam'' must have the same length.'); end
otherwise
end


