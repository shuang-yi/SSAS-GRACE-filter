function iia = isint(a,tol)

% ISINT checks whether the elements of a matrix are integer. 
% Returns a mask, which can be used as a single Boolean, in the sense that 
% IIA is true when IIA = ones(size(A)). Otherwise IIA is false.
%
%    IIA = ISINT(A)         - strong form
%    IIA = ISINT(A,tol)     - weak form
%
% SEE ALSO:
%    FLUSHINT

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    1997-01-15: NS, definitions changed (Functionality still the same).
%    1994-05-02: NS, initial version
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

if nargin == 1
   iia = ~rem(a,1);
else
   iia = abs(round(a)-a) < tol;
end
