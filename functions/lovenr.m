function kn = lovenr(n)

% LOVENR gives the LOVE number of the elastic earth for a certain degree n
%
% IN:
%   n ........ spherical harmonic degree (up to 200)
%
% OUT:
%   kn ....... LOVE number of degree n
%
% REMARKS:
%   The elastic LOVE numbers are taken from the paper by WAHR et al., 
%   "Time variability of the earth's gravity field: hydrological and 
%   oceanic effects and their possible detection using GRACE",  
%   JGR, Vol. 103, No. B12, p 30205-30229, 1998

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Matthias ROTH (MR), GIS, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2002-04-16: MR, initial version
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


% data points for LOVE numbers as provided by Wahr for selected degrees
% where missing degrees can be interpolated linearly
l  = [0   1     2     3     4     5    6    7    8    9   10   12   15   20   30   40   50   70  100  150  200]';
kl = [0 270 -3030 -1940 -1320 -1040 -890 -810 -760 -720 -690 -640 -580 -510 -400 -330 -270 -200 -140 -100  -70]'./1e4;

% here we go
kn = interp1(l,kl,n,'linear',0);
