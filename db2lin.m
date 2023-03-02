function varargout = db2lin(varargin)
% DESCRIPTION y = db2lin(x)
%  Converts dB to linear.
%  Takes any number of arguments.  -Inf will map to 0.
% INPUT
%  x1  -- Real, matrix
%  x2  -- Real, matrix
%  x3  -- Real, ... 
% OUTPUT
%  y1 -- As 1st input, but coverted.
%  y2 -- As 2nd input, but coverted.
%  y3 -- As 3rd input, but ...
% INPUT
%  Any number of any real matrices.
% OUTPUT
%  y --  The input arguments converted.  
% TRY
%  db2lin(3), db2lin([0 -inf]), db2lin([0 0], -3+zeros(2,2,2)) 
% SEE ALSO
%  lin2db

for i = 1:nargin
 varargout{i} = 10 .^ (varargin{i}/10);
end
