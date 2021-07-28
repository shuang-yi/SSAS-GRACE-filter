function [cs, rows, cols] = sc2cs(field)

% SC2CS(FIELD) converts the rectangular (L+1)x(2L+1) matrix FIELD, containing
% spherical harmonics coefficients in /S|C\ storage format into a 
% square (L+1)x(L+1) matrix in |C\S| format.
%
% IN:
%    field .... the rectangular (L+1)x(2L+1) matrix FIELD, containing
%               spherical harmonics coefficients in /S|C\ storage format
%
% OUT: 
%    cs ....... square (L+1)x(L+1) matrix in |C\S| format
%    rows ..... the number of rows for matrix FIELD
%    cols ..... the number of columns for matrix FIELD

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    Matthias WEIGELT (MW), DoGE, UofC
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    20??-07-10: MW, - file is now insensitive to CS-format input
%                    - output rows and cols added
%    1994-07-22: NS, initial version
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

[rows,cols] = size(field);
if (rows~=cols) && (cols~=2*rows-1)
    error('Input neither in cs nor in sc format');
elseif cols == rows
    cs = field;
else
    lmax = rows-1;
    c    = field(:,lmax+1:2*lmax+1);
    s    = [zeros(lmax+1,1) field(:,1:lmax)];
    cs   = tril(c) + triu(rot90(s),1);
end
