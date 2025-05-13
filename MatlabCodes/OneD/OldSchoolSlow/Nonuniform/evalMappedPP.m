function out = evalMappedPP( pp, width, offset, x )
% EVALMAPPEDPP - Evaluates the supplied piecewise polynomial
%   scaled by width and shifted by offset for the real world 
%   coordinates of x
%   values beyond the left-most and right-most breakpoints
%   are assumed to be zero.

[brk, c, l, d, k] = unmkpp( pp );

minVal = brk(1);
maxVal = brk(size(brk,2)) ;

%x = x - offset;
x = offset - x;
x = x / width;
out = ppval( pp, x ) .* (x >= minVal & x <= maxVal);