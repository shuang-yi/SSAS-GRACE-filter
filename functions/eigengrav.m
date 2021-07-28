function tf = eigengrav(lmax, fstr, h, const)

% EIGENGRAV(lmax, fstr, h, type) returns the isotropic spectral transfer
% (or: eigenvalues) of several gravity related quantities. 
% Upward continuation may be included.
%    
% IN:
%    lmax .. spherical harmonic degree                             [n x 1]
%    fstr .. string, denoting the functional under consideration:
%            'none', 
%            'geoid',
%            'dg', 'gravity' ... gravity anomaly,
%            'potential', 
%            'tr' .............. gravity disturbance, 
%            'trr' ............. (d^2/dr^2)
%            'slope' ........... size of surface gradient, 
%            'water' ........... equivalent water thickness, 
%            'smd' ............. surface mass density.
%            When fstr is left empty, a menu pops up.
%    h ..... (default: 0), height above Earth mean radius [m].   
%    type .. type of constants loaded for the calculations
%            1. If left empty constants of the GRS80 ellipsoid will be used
%            2. A 2 element array with GM and semi-major axis of the ellipsoid
%               can be given [GM ae]
%            3. A file with the pathname as a string can be given where GM and
%               ae must have the same variable names.
%
% OUT:
%    tf .... transfer. Size and shape equal to lmax. Units are respectively 
%            [none], [m], [mGal], [mGal], [E], [m^2/s^2], [rad], [m], [kg/m^2].
%                                                                  [n x 1]
%
% USES:
%    upwcon, lovenr, uberall/constants, uberall/isint
%
% SEE ALSO:
%    UPWCON

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAGP, TU Munich
%    Matthias WEIGELT (MW), DoGE, UofC
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    Balaji DEVARAJU (BD), IfE, LUH
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2015-05-22: MR, replace if ... elseif by switch, replace isstr by
%                    ischar
%    2014-11-14: BD, changes in testing of the input constants
%                    brushed up help text
%    2014-01-14: MR, brush up help text, beautify code
%    2013-02-13: MR, change function names, brush up comments
%    2013-01-30: MA, comments/removing of smoothing option
%    2007-02-15: MW, adding water equivalent thickness and surface mass density
%    1996-09-27: NS, - menu-driven input to get FSTR.
%                    - slope option added.
%    1995-05-03: NS, initial version
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

% defaults
if nargin < 4, const= []; end
if nargin < 3, h    = 0;  end
if nargin < 2, fstr = ''; end
if isempty(h), h    = 0;  end

% general checks
if nargin < 1,        error('Gimme more'),           end
if ~isscalar(h),      error('H should be scalar.'),  end
if ~ischar(fstr),     error('FSTR must be string.'), end 
if any(~isint(lmax)), error('L must be integer.'),   end
if ~isvector(lmax),   error('L must be vector.'),    end
if min(lmax) < 0,     error('Negative L occurs.'),   end

% Get the proper functional string
fstr = lower(deblank(fstr));
if isempty(fstr)
   fstrmat = char('none', 'geoid', 'dg', 'tr', 'trr', 'potential', 'slope', 'water');
   fchoice = menu('CHOOSE GRAVITY FUNCTIONAL',...
               'None', 'Geoid', 'Gravity Anomaly',...
               'Gravity Disturbance', 'Radial Grav. Gradient',...
               'Potential', 'Slope', 'Equivalent Water Thickness');
   fstr    = deblank(fstrmat(fchoice, :));
end

% load necessary constants
if isempty(const)
    constants;
elseif numel(const) == 2
    GM = const(1);
    ae = const(2);
elseif ischar(const)
    eval([const ';']);
else
    error('Please verify the input constants.')
end

% treat eigenvalue part and dimensioning
r = ae + h;
switch fstr
    case 'none'
        tf = ones(size(lmax));                      % []
    case 'geoid'
        tf = ones(size(lmax)) * r;                  % [m]
    case 'potential'
        tf = ones(size(lmax)) * GM/r;		    	% [m^2/s^2]
    case {'gravity', 'dg'}
        tf = (lmax-1) * GM/r/r * 1e5;		    	% [mGal]
    case 'tr'
        tf = (lmax+1) * GM/r/r * 1e5;		     	% [mGal]
    case 'gr'
        tf = (lmax+1) * GM/r/r * 1e8;		     	% [microGal]
    case 'trr'
        tf = (lmax+1) .* (lmax+2) * GM/r/r/r * 1e9;	% [E]
    case 'slope'
        tf = sqrt(lmax .* (lmax+1));			    % [rad]
    case 'water'
        tf = (2.*lmax+1) ./ (1+lovenr(lmax)) * r * 5.517 / 3;  % [m]
    case 'ewh'
        tf = (2.*lmax+1) ./ (1+lovenr(lmax)) * r * 5.517 / 3 * 1D2;  % [cm]
    case 'smd'
        tf = (2.*lmax+1) ./ (1+lovenr(lmax)) * r * 5517 / 3;   % [kg/m^2]
    otherwise
        error('Requested functional FSTR not available.')
end

% treat upward continuation part: when r = R + h is used in the dimensioning
% it is equal for all functionals
if h > 0
   tf = upwcon(lmax, h, const) .* tf;			% add upw. cont. to transfer
end


