function uc = upwcon(degree,height,const)

% UPWCON(degree,height,const) returns the upward continuation (R/r)^l, in which
%
% uc = upwcon(degree,height)
% uc = upwcon(degree,height,const)
%
% INPUT
% degree - Spherical harmonic degree, must be integer [scalar/vector]
% height - Height above mean Earth radius [m] [scalar/vector]
% const  - Reference ellipsoid constants: Gravitational constant (GM) and 
%          semi-major axis (ae) given as [GM ae]. See uberall/constants for
%          units. (optional)                        - def: GRS80 constants
%
% OUTPUT
% uc    - Upward continuation terms.
% 
% USES: 
%    uberall/lying, uberall/standing, uberall/constants
%
% REMARKS:
%    If both degree and height are vectors, degree will be(come) a row vector
%    and height a column vector.
%    load necessary constants

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    Matthias WEIGELT (MW), GI, Uni Stuttgart
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    Balaji DEVARAJU (BD), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%   2014-03-04: BD, brush up help text
%   2013-02-13: MR, change function names, brush up comments
%   2013-01-29: MA, comments
%   1994-04-22: NS, initial version
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

if nargin < 3, const = []; end

% load necessary constants
if isempty(const)
    constants;
elseif numel(const) == 2
    ae = const(2);
else
    eval([const ';']);
end


if nargin < 2, height = 0; end
if length(height) > 1
   degree = lying(degree);
   height = standing(height);
   rr = ae*ones(size(height)) ./ (ae+height);
   uc = (rr * ones(size(degree))) .^ (ones(size(height)) * degree);
else    
   rr = ae / (ae+height);
   uc = rr.^degree;
end
