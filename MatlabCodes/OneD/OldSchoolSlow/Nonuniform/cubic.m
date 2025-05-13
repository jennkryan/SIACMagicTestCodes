function out = cubic( x )
%out = x .* x .* x;
out = -1/3.*x.*(2.*x-1).*(x-1);