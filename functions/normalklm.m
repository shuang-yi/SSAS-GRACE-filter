function nklm = normalklm(lmax,typ)

% NORMALKLM returns an ellipsoidal normal field
% consisting of normalized -Jn, n=0,2,4,6,8
%
% IN:
%    lmax ....... maximum degree
%    typ ........ either 'wgs84' (equipotential ellipsoid), default,
%                        'grs80',
%                 or     'he' (hydrostatic equilibrium ellipsoid)
% OUT:
%    nklm ....... normal field in CS-format (sparse)
%
% REMARKS:
%    .J2,J4 values for hydrostatic equilibrium ellipsoid from Lambeck (1988)
%    "Geophysical Geodesy", p.18

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    Matthias WEIGELT (MW), DoGE, U Calgary
%    Balaji DEVARAJU (BD), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2014-03-09: BD,  brushed up help text
%    xxxx-xx-xx: MW,  added 'wgs84'
%    2000-03-24: NS,  brush up help text
%                   - add ellipsoid in hydrostatical equilibrium 
%                   - remove all dependencies
%    1996-10-30: NS, output will be a sparse matrix.
%    1995-02-27: NS, initial version
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

% checks and defaults
narginchk(1,2);

if nargin < 2,       typ = 'wgs84'; end     % default type of ellipsoid
if numel(lmax) > 1,  error('LMAX should be scalar.'),   end
if rem(lmax,1) ~= 0, error('LMAX should be integer.'),  end
if lmax < 0,         error('LMAX should be positive.'), end

% initialize 
switch lower(typ)
    case 'wgs84'
        J2     =  1.08262982131e-3;     % earth's dyn. form factor (= -C20 unnormalized)
        J4     = -2.37091120053e-6;     % -C40 unnormalized
        J6     =  6.08346498882e-9;     % -C60 unnormalized
        J8     = -1.42681087920e-11;    % -C80 unnormalized
        jcoefs = [1 -J2 -J4 -J6 -J8]';
        l      = (0:2:min(lmax,8))';
    case 'grs80'
        J2     =  1.08263e-3;          % earth's dyn. form factor (= -C20 unnormalized)
        J4     = -2.37091222e-6;     % -C40 unnormalized
        J6     =  6.08347e-9;        % -C60 unnormalized
        J8     = -1.427e-11;         % -C80 unnormalized
        jcoefs = [1 -J2 -J4 -J6 -J8]';
        l      = (0:2:min(lmax,8))';
    case {'he','hydro'}
        J2     = 1.072618e-3;		% earth's dyn. form factor (= -C20 unnormalized)
        J4     = 0.2992e-5;     	% -C40 unnormalized
        jcoefs = [1 -J2 -J4]';
        l      = (0:2:min(lmax,4))';   
    otherwise
        error(['Unknown type of ellipsoid: ' typ])
end

% calculate
coefs = jcoefs(1:length(l))./sqrt(2*l+1);
nklm  = sparse(l+1,ones(size(l)),coefs,lmax+1,lmax+1);

