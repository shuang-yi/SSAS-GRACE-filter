function [varargout] = Legendre_mex(varargin) 

% Computes fully normalized Legendre functions of the first kind and their
% first and second derivatives (tested for Lmax <= 10000).
%
% For a degree/order Lmax > ~1800 the addition theorem is breached. To
% avoid the breach, this function follows Fukushima's concept of X-numbers 
% to avoid a breach of the addition theorem for higher degree/order of lmax
% (see: Fukushima T (2012), "Numerical computation of spherical harmonics
% of arbitrary degree and order by extending exponent of floating point
% numbers", Journal of Geodesy, Vol. 86, pp. 271-285).
% The function switches automatically to the stable (but slightly slower)
% X-number method for Lmax >= 1200 (just to be on the safe side).
%
% For speed, the function was implemented as mex-file in C. 
% Compile it by calling:
%
%    "mex Legendre_mex.c legendre_od.c x_numbers.c"
%
% (Prior to calling, install a Matlab-supported C-compiler).
%
% The function merges/replaces the following functions: 
%   Legendre_mex, Legendre_x_mex, legendreP,
%   plm, legendreP, legendreP_Xnum.
% The latter three functions will be kept to ensure that people without
% C-compiler can calculate Legendre functions.
%
% This function supports (nearly) the former PLM calling scheme 
% (PARAMETERS MUST BE IN THE ORDER DESCRIBED HERE!):
%
% PLM mode:
% IN:
%    l .......... scalar or lying vector: degrees of Legendre functions,
%                 [integer],
%    m .......... (optional) scalar: order of Legendre functions; if not
%                 given, m = 0 is assumed, [integer],
%    thetaRAD ... lying vector of co-latitudes [rad].
%    'plm' ...... (optional) use Nico's PLM style
%    'xnum'/'noxnum' ... (optional) force use of X-numbers/switch off
%                 X-numbers. If left out, the function determines itself if
%                 X-numbers are necessary.
% OUT:
%    nL ......... a matrix of size (length(thetaRAD) x length(l)),
%                 every column corresponds to a value of thetaRAD, per column
%                 the order of coefficients corresponds to the order of the l 
%                 and thetaRAD vectors.
%    dnL ........ like nL, just for the first derivatives
%    ddnL ....... like nL, just for the second derivatives
%               
% SPEED_OD mode:
% IN:
%    Lmax ....... maximum degree/order of Legendre functions, [integer],
%    thetaRAD ... (lying vector) of co-latitudes in [rad]
%    'speed_od' ..(necessary) to use the optimized computation algorithm
%    'xnum'/'noxnum' ... (optional) like above in PLM mode.
%
% OUT: 
%    nL ......... a matrix of size ((Lmax + 1) * Lmax / 2 + 1) times 
%                 length(thetaRAD), every column corresponds to a value of
%                 phi, per column the order of coefficients is: 
%                 00, 10, 20, 30, ... 11, 21, 31, ... 22, 32, ... 
%    dnL ........ like nL, just for the first derivatives
%    ddnL ....... like nL, just for the second derivatives
%
% EXAMPLES:
%    nL = Legendre_mex(0:10, [20 25] * pi / 180);
%    [nL, dnL] = Legendre_mex(0:10, 5, [20 25] * pi / 180);
%
%    nL = Legendre_mex(200, [10 20 30 40] * pi / 180, 'speed_od');
%    [nL, dnL] = Legendre_mex(200, [10 20 30 40] * pi / 180, 'speed_od');
%    [nL, dnL, ddnL] = Legendre_mex(200, [10 20 30 40] * pi / 180, 'speed_od');
%
% SEE ALSO:
%    plm, legendreP, legendreP_Xnum, legendreIndex

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2015-07-17: MR, major revision, rearrange code to other files,
%                    reprogram parameter interface to also force/not allow
%                    use of X-numbers
%    2014-10-08: MR, major revision, merge with other Legendre functions 
%    2013-11-21: MR, add license information
%    2013-11-19: MR, add calculation of dnL and ddnL
%    2013-11-18: MR, remove mistake in filename
%    2013-07-02: MR, speed-up, add help text
%    2013-06-17: MR, initial version
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

if exist('Legendre_mex', 'file') ~= 3
    error('Please, compile first by calling function ''mex_compile.m''!');
end
