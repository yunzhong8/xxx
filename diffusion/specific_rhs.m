function f = specific_rhs(x,y,nel)
%unit_rhs   unit RHS forcing function
%   f = specific_rhs(x,y,nel)
%   input
%          x          x coordinate vector
%          y          y coordinate vector 
%          nel        number of elements  
%   IFISS function: DJS; 28 February 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 
% f=ones(nel,1);

% Define the RHS forcing function
f = 32 * pi^2 * sin(4 * pi * x) .* sin(4 * pi * y);
% Ensure the output is a column vector for compatibility
f = f(:);
return