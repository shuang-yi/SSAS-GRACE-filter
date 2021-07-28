function f = ispec(a,b)

% ISPEC(A,B) returns the function F from the spectra A and B.
%
%IN:
%    a ...... cosine coefficients 
%    b ...... sine coefficients          
%             
%             a and b are defined by:
%             f(t) = a_0 + SUM_(i=1)^n2 a_i*cos(iwt) + b_i*sin(iwt)
%   
%             with w = ground-frequency and n2 half the number of samples (+1).
%             Note that no factor 2 appears in front of the sum.
% 
% OUT:
%    F = ISPEC(A,B) considers A the cosine- and B the sine-spectrum.
%    F = ISPEC(S) assumes S = [A B].
%    If  A and B are matrices, Fourier operations are columnwise.
% 
% USES: 
%    spec
%
% SEE ALSO:
%    SPEC, FFT

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    1994-06-29: NS, initial version
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
   if min(size(a)) == 2 && size(a, 2) == 2
      a = a(:);				        % Put A upright.
   end
   m = size(a, 2);
   if rem(m, 2) ~= 0 
      error('If one input argument: number of columns must be even.') 
   end
   b = a(:, m/2+1:m);				% Split into A and B spectra.
   a = a(:, 1:m/2);
elseif nargin == 2
   if sum(size(a) == size(b)) ~= 2
      error('Sizes of A and B don''t match!')
   end 
   if min(size(a)) == 1
      a = a(:);				        % Put vectors A and B upright.
      b = b(:);
   end
else
   error('Check your input!')
end
[n2, ~] = size(a);

a(1,:) = a(1,:) * 2;
if abs(b(n2, :)) < 1e-10 			% case where n must be even
   n   = 2*n2 - 2;				
   a(n2, :) = a(n2, :) * 2;			% simulate 100% aliasing
   fs  = (a - 1i * b) / 2;	
   fs  = [fs; conj(fs(n2-1:-1:2, :))] * max(n,1);% 'max' required for n2 = 1
else			                    % case where n must be odd
   n   = 2*n2 - 1;				
   fs  = (a - 1i*b) / 2;	
   fs  = [fs; conj(fs(n2:-1:2, :))] * n;
end

f = real(ifft(fs));
