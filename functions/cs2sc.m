function [sc, rows, cols] = cs2sc(field, backval)

% CS2SC(FIELD,backval) converts the square (L+1)x(L+1) matrix FIELD, containing
% spherical harmonics coefficients in |C\S| storage format into a 
% rectangular (L+1)x(2L+1) matrix in  /S|C\format.
%
% IN:
%    field .... the square (L+1)x(L+1) matrix FIELD , containing
%               spherical harmonics coefficients in |C\S| storage format
%    backval .. optional input and describes the matrix entries, 
%               where m > l. Default is 0!
% OUT: 
%    sc ....... rectangular (L+1)x(2L+1) matrix in  /S|C\format
%    rows ..... rows of rectangular (L+1)x(2L+1) matrix  in  /S|C\format
%    cols ..... columns of rectangular (L+1)x(2L+1) matrix  in  /S|C\format

% ----------------------------------------------------------------------------
% project: SHBundle 
% ----------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    Matthias WEIGELT (MW), DoGE, UofC
%    <bundle@gis.uni-stuttgart.de>
% ----------------------------------------------------------------------------
% revision history:
%    2010-07-01: MW, file is now insensitive to SC-format input, 
%                added output rows and cols
%    1999-07-01: MW, V5 adaptation, eliminiation of TRAPSTRIP
%    1994-07-22: NS, initial version
% ----------------------------------------------------------------------------
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
% --------------------------------------------------------------------------

if nargin == 1, backval = 0; end
[rows,cols, Nmonth] = size(field);

if (rows~=cols) && (cols~=2*rows-1)
    error('Input neither in cs nor in sc format');
elseif cols==2*rows-1
    sc = field;
else
    sc = zeros(rows,2*rows-1,Nmonth);
    for ii = 1:Nmonth
        lmax = rows-1;
        c    = tril(field(:,:,ii));
        s    = rot90(triu(field(:,:,ii),1),-1);
        mask = backval*ones(lmax+1,2*lmax+1);
        a    = fliplr(triu(mask(:,1:cols-1)));
        b    = triu(mask(:,cols:2*lmax+1),1);
        sc(:,:,ii)   = [a b] + [s(:,2:lmax+1) c];
    end
end
