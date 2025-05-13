function out = periodicEval( u, boundaries, x )
%periodicEval    Evaluates the Legendre-based piecewise-polynomial
%                 expression at x, in a periodic fashion
%                 Evaluates for arbitrary interval size
%
%   periodicEval( u, boundaries, x)
%       u, the 2D matrix representing the modes of each element
%       boundaries, a list of increasing real values definining the
%       elements of which u is defined
%       x, the value in world coordinates to evaluate
%
%   Make the following assumptions:
%       1 - the input function (u) is piecewise polynomial in legendre
%           basis
%       2 - The domain of u is [ boundaries, boundaries( last ) ]
%       3 - The function u is periodic with period domainWidth
%       4 - Assume the intervals are defined on [lower, upper)


domainMin = boundaries(1);
domainMax = boundaries( size( boundaries, 2) );
domainWidth = domainMax - domainMin;
numIntervals = size( boundaries, 2) - 1;

evalX = [];
% For each x -
%   Map x into the real domain
%   Find the interval that x lies on
%   map x into the [-1, 1] legendre domain on that interval
%   return the pair, (x, interval)
for i = 1:size(x,2)
    tempX = x(i);
    % map into real
    while (tempX < domainMin )
        tempX = tempX + domainWidth;
    end
    while ( tempX > domainMax )
        tempX = tempX - domainWidth;
    end
    % find interval
    interval = 1;
    for interval = 2:numIntervals+1
        if ( boundaries( interval ) > tempX )
            break;
        end
    end
    % interval is [interval - 1, interval )
    % map to legendre domain
    intervalWidth = boundaries(interval) - boundaries( interval - 1 );
    tempX = (tempX - boundaries( interval - 1)) / intervalWidth * 2 - 1;
    % return pair
    evalX(:,i) = [ tempX; interval - 1 ];
end

% evaluate values
out = zeros(1, size(x,2) );
modes = size(u,2);
for L = 1:modes
    for v = 1:size(x,2)
        out(v) = out(v) + u(evalX(2,v), L) * sqrt(L-0.5)*JacobiPoly(L-1, evalX(1,v), 0.0, 0.0);
    end
end
   
