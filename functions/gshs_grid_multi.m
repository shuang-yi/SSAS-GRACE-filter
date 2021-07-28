function f = gshs_grid_multi(field, lamRAD, phiRAD, a_E, varargin)

% GSHSAG calculates a global spherical harmonic synthesis for any grid 
% defined by lam and phi (both vectors). The radius must be scalar.
%
% f = gshsag(field, lamRAD, phiRAD, a_E)
%
% IN:
%    field ....... gravity field in |c\s| or /s|c\ format
%    lamRAD ...... longitude [rad]                                   [n, 1]
%    phiRAD ...... latitude  [rad]                                   [m, 1]
%    a_E ......... semi major axis of the Earth rotational ellipsoid [1, 1]
% OPTIONAL:
%    'height' .... (default: 0), height [m]                        [scalar]
%    'max_lm' .... maximum degree/order (default: determined from field)
%                                                                  [scalar]
%    'quant' ..... optional argument, defining the field quantity: [string]
%                  - 'potential' ... (default), potential [m^2/s^2], needs 'GM'
%                  - 'tr' .......... gravity disturbance [mGal], needs 'GM'
%                  - 'trr' ......... 2nd rad. derivative [E], needs 'GM'
%                  - 'none' ........ coefficients define the output
%                  - 'geoid' ....... geoid height [m] 
%                  - 'dg', 'gravity' gravity anomaly [mGal], needs 'GM'
%                  - 'slope' ....... slope [arcsec]
%                  - 'water' ....... equivalent water thickness [m] 
%                  - 'smd' ......... surface mass density
%    'sub_WGS84' . (default: true) determines whether WGS84 is subtracted.
%    'GM' ........ geocentric gravitational constant GM
%    'legendre' .. (default: 'plm') Legendre functions algorithm:
%                  - 'plm' ... unstable (for d/o > ~1800) PLM
%                  - 'mex' ... fast, X-number stabilized LEGENDRE_MEX
%    'waitbar' ... if set, shows a waitbar (default: false)
%    'curvature' . (default: false) consider curvature of the chosen quantity [bool]
%
% OUT:
%    f ....... field quantity                                       [n, m]
%
% EXAMPLE: see SHBUNDLE/example/example_gshs_grid.m
%
% USES:
%    vec2cs, normalklm, eigengrav, plm, Legendre\_mex, cs2sc, sc2cs, 
%    UBERALL/checkcoor, UBERALL/twaitbar
%
% SEE ALSO:
%    GSHS_, GSHS_PTW, REGIONALSHS, GSHA

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Matthias WEIGELT (MW), DoGE, UofC
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2018-11-27: MA, extra warning if a reference field is subtracted
%                    (request of external users; avoidable via getopt.m)
%    2015-05-22: MR, reprogram parameter interface (hence, rename function), 
%                    revise help text and code
%    2014-01-15: MR, revise help text, beautify code
%    2013-02-13: MR, change function names, brush up comments
%    2013-01-30: MA, comments/removing of smoothing option
%    2013-01-23: MA, input of plm in radian
%    2012-03-06: MW, add curvature calculations
%    2009-03-04: MW, change input co-latitude -> latitude
%    2005-10-25: MW, initial version
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
% INPUT CHECK and PREPARATION
%--------------------------------------------------------------------------


%% define defaults and check optional parameters
defaultparams = {'height', 0;
                 'max_lm', inf;
                 'quant', 'potential';
                 'sub_wgs84', true;
                 'gm', nan;
                 'legendre', 'plm';
                 'waitbar', false;
                 'curvature', false}; 
params = getopt(defaultparams, false, varargin);  

% check Legendre function algorithm
plm_func = false;
switch params.legendre
    case 'mex'
        if exist('Legendre_mex', 'file') == 3; % check if compiled Legendre_mex exists
            plm_func = @Legendre_mex;
        else
            warning('Legendre_mex is not compiled! Falling back to plm...');
            plm_func = @plm;
        end
    case 'plm'
        plm_func = @plm;
    otherwise
        error('Unknown Legendre function.');
end

% field size determination, rearrange field and subtract reference field
[row, col, Nmonth] = size(field);

if col == row              % field in |C\S|-format 
    field = cs2sc(field);  % convert to /S|C\-format
elseif col ~= 2 * row - 1
   error('Input "field" not in cs or sc format');
% else: field is already in /S|C\-format
end

if params.max_lm > (row - 1)     % desired max_lm is bigger than what the field provides
    params.max_lm = row - 1;     % adjust max_lm 
elseif params.max_lm < (row - 1) % if max_lm is smaller than what the field provides
    error('error here, this function needs to be added')
    field = field(1:(params.max_lm + 1), (row - params.max_lm):(row + params.max_lm)); % adjust field
    row = size(field, 1); % update row
% else: everything is ok
end

if params.sub_wgs84
    for ii = 1:Nmonth
        field(:,:,ii) = field(:,:,ii) - cs2sc(full(normalklm(params.max_lm, 'wgs84')));
    end
    if params.display_warning == true
        warning('A reference field (WGS84) is removed from your coefficients')
    end
end

% check if r, lamRAD and phiRAD are vectors
[lamRAD, phiRAD, ~] = checkcoor(lamRAD, phiRAD, 0, 0, 'grid');
theRAD              = (pi/2 - phiRAD);
if ~isscalar(params.height)
    error('''height'' must be scalar.');
end

% prepare l, m, often used variables
m       = (0:params.max_lm);
l       = m';
mlam    = (lamRAD * m)';           
len_the = length(theRAD);

% apply transfer function
transf  = eigengrav(l, params.quant, params.height, [params.gm, a_E]);
for ii = 1:Nmonth
    field(:,:,ii)   = diag(transf(:)) * field(:,:,ii);
end
 
%----------------------------------------------------------------------------
% CALCULATION
%----------------------------------------------------------------------------
if params.waitbar
    WAIT = twaitbar('init', [], 'gshs_grid: calculating');
end;

if params.curvature
    TAp  = []; TAp(len_the, row) = 0; % faster than TAp = zeros(len_the, row)
    TBp  = []; TBp(len_the, row) = 0;
    TAl  = []; TAl(len_the, row) = 0;
    TBl  = []; TBl(len_the, row) = 0;
    TApp = []; TApp(len_the, row) = 0; 
    TBpp = []; TBpp(len_the, row) = 0;
    TAll = []; TAll(len_the, row) = 0; 
    TBll = []; TBll(len_the, row) = 0;
    TApl = []; TApl(len_the, row) = 0; 
    TBpl = []; TBpl(len_the, row) = 0;
    
    % unwrapping to avoid if ... else for order m = 0 (i.e. Snm = 0)
    Cnm = field(:, row);              % get Cnm coefficients for order 0
    [P, dP, ddP] = plm(l, 0, theRAD); % calc fully normalized Legendre Polynoms
    TAp(:, 1)  = -dP * Cnm; % all variables multiplied with Snm are 0, of course  
    TBl(:, 1)  =  -P * Cnm;
    TApp(:, 1) = ddP * Cnm;
    TAll(:, 1) =  -P * Cnm;
    TBpl(:, 1) =  dP * Cnm;
    for m = 1:(row - 1)
        Cnm = field(:, row + m);          % get Cnm coefficients for order m
        Snm = field(:, row - m);          % get Snm coefficients for order m
        [P, dP, ddP] = plm_func(l, m, theRAD); % calc fully normalized Legendre Polynoms
        TAp(:, m+1)  = -dP * Cnm;  TBp(:, m+1)  = -dP * Snm;
        TAl(:, m+1)  =   P * Snm;  TBl(:, m+1)  =  -P * Cnm;
        TApp(:, m+1) = ddP * Cnm;  TBpp(:, m+1) = ddP * Snm;
        TAll(:, m+1) =  -P * Cnm;  TBll(:, m+1) =  -P * Snm;
        TApl(:, m+1) = -dP * Snm;  TBpl(:, m+1) =  dP * Cnm;
        if params.waitbar; WAIT = twaitbar(m / (row - 1), WAIT); end;
    end
    
    % sum the partial derivatives
    cosmlam = cos(mlam);
    sinmlam = sin(mlam);
    fp  = TAp  * cosmlam + TBp  * sinmlam; fp2 = fp.^2;
    fl  = TAl  * cosmlam + TBl  * sinmlam; fl2 = fl.^2;
    fpp = TApp * cosmlam + TBpp * sinmlam;
    fll = TAll * cosmlam + TBll * sinmlam;
    fpl = TApl * cosmlam + TBpl * sinmlam;
    
    % now do the final summation
    f = ((1 + fp2) .* fll - 2 .* fp .* fl .* fpl + (1 + fl2) .* fpp) ./ (2 .* (realsqrt(1 + fp2 + fl2)).^3);
    
else
    TA = []; TA(len_the, row, Nmonth) = 0; % faster than TA = zeros(length(idx), row);
    TB = []; TB(len_the, row, Nmonth) = 0;   

    % unwrapping to avoid if ... else for order m = 0
    Cnm = field(:, row, :);                % get Cnm coefficients for order 0
    for ii = 1:Nmonth
        TA(:, 1, ii) = plm(l, 0, theRAD) * Cnm(:,:,ii); % for m = 0: Snm = 0, hence TB also = 0
    end
    for m = 1:row-1
        Cnm = field(:, row + m, :);     % get Cnm coefficients for order m
        Snm = field(:, row - m, :);     % get Snm coefficients for order m
        P = plm_func(l, m, theRAD);  % calc fully normalized Legendre Polynoms
        % ------------------------------------------------------------------
        for ii = 1:Nmonth
            TA(:, m + 1, ii) = P * Cnm(:,:,ii);
            TB(:, m + 1, ii) = P * Snm(:,:,ii);
        end
        if params.waitbar; WAIT = twaitbar(m / (row - 1), WAIT); end;
    end
    
    % now do the final summation
    for ii = 1:Nmonth
        f(:,:,ii) = TA(:,:,ii) * cos(mlam) + TB(:,:,ii) * sin(mlam);
    end
    f = squeeze(f);
end

if params.waitbar; twaitbar('close', WAIT); end;

