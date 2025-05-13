function out = nonUniformPost( elements, boundaries, BSorder, RS, kernelScale )
% elements = approximation modes (size(numElements,order))
% boundaries = element boundaries size(numElements+1)
% evalPoints = number of evaluation points per element
% ell = B-Spline order (smoothness+2)
% kernelScaling = scaling for the convolution kernel.

% polynomials are all of degree k
%   map k+1 eval points (quadrature zeros) to each interval
%   that means there will be (k+1) * numElements = N points
%   output is a 2 X N array, first row contains x values in world coords
%   second row contains the corresponding y values
%%%%%%%%%%%
% set up criteria
[numElements, order] = size(elements);
kapprox = order -1;
[z, w] = JacobiGZW(order, 0.0, 0.0);
evalPoints = 2*order;
zEval = zeros(evalPoints,1);
wEval = zeros(evalPoints,1);
zEval(1:order) = 0.5.*(z-1);
zEval(order+1:evalPoints) = 0.5.*(z+1);
NumP = size(zEval, 1);
%%%%%%%%%%%

%%%%%%%%%%%
% set x values -- this could alternatively be supplied externally
sampleCount = numElements * evalPoints; %order;
out = zeros(2, sampleCount );
offset = 1;
last = offset + evalPoints-1;
for i = 1:numElements
    h = boundaries( i + 1) - boundaries(i);
    out(1, offset:last) = h * 0.5 * (zEval' + 1) + boundaries( i );
    offset = last + 1;
    last = offset + evalPoints-1;
end
%%%%%%%%%%%
%%%%%%%%%%%
% USE THIS IF WE WANT TO USE THE LARGEST ELEMENT IN THE MESH FOR THE KERNEL
% SCALING.
if kernelScale ==0
    for n=1:numElements
        kernelScale = max(kernelScale,boundaries( n + 1 ) - boundaries( n ));
        elementWidth(n) = boundaries(n+1) - boundaries(n);
    end
end
%%%%%%%%%%%
% set up kernel statistics -- bspline for evaluation, etc.
bs = getBSplinePP( BSorder );     % this is ORDER of the bspline -- i.e. order 1 = degree 0 (constant bspline)

%%%%%%%%%%%
% NUMBER OF B-SPLINES
if RS == 0
    rspline = min(order+BSorder,2*order-1);
else
    rspline = RS;
end
if 2*floor(0.5*rspline) == rspline
    rspline = rspline+1;
end
fprintf('Using %i B-splines of order %i\n',rspline,BSorder)
% convolution
%  NUMBER OF GAUSS POINTS FOR EVALUATION OF POLYNOMIAL
kgl = ceil(0.5*(kapprox+BSorder+1));
[z, w] = JacobiGZW(kgl, 0.0, 0.0);  % need enough quad points to support k + BSOrder degree

%%%%%%%%%%%
%  Kernel stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
% Using symmetric kernel ==> lambda = 0
lambda=0;
%coeff = KernelCoeffs( rspline, BSorder-1, lambda)
coeff = getKernelCoeffsMatrix( rspline, BSorder-1);
kernelWidth = rspline + BSorder-1;
halfWidth = kernelWidth * 0.5;



compactSupport = [lambda - halfWidth, lambda + halfWidth] * kernelScale;
kernelBreaks = (lambda - halfWidth:lambda + halfWidth) * kernelScale;


outIndex = 1;   % where to write the x value to in the output
for n = 1:numElements
    elementWidth = boundaries( n + 1 ) - boundaries( n );
    for xIndex = 1:evalPoints
        % find all the breakpoints for integration
        offset = out(1, (n - 1) * evalPoints + xIndex);  % get the breakpoint
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        localCompactSupport = offset + compactSupport ;   % scale compact support to current evaluation point
        localKernelBreaks = offset + kernelBreaks ;  % kernel breaks wrt current evaluation point.
        uhatBreaks = funcBreakPoints_periodic( boundaries, localCompactSupport );  % in world coordinates (outside of interval definition)
        breaks = [localKernelBreaks uhatBreaks];
        breaks = sort(breaks);

        % integrate product of convolved kernel and uhat
        intVal = 0;
        for i = 1:(size(breaks,2) - 1)
            a = breaks( i );
            b = breaks( i + 1 );
            %   map quadrature points
            quadZ = (z' + 1) * 0.5 * (b - a) + a;
            % evaluate uhat and kernel at quadrature points
            uhatVal = periodicEval( elements, boundaries, quadZ );
            %kernelVal = evalKernelPP( bs, rspline, BSorder-1, lambda, elementWidth, offset, coeff, quadZ);
            kernelVal = evalKernelPP( bs, rspline, BSorder-1, lambda, kernelScale, offset, coeff, quadZ);
            %   perform quadrature
            intVal = intVal + (b-a)*0.5*dot( w, kernelVal .* uhatVal );   
        end
        %out(2,outIndex) = intVal / elementWidth;  % rescale the integral by the element width
        out(2,outIndex) = intVal / kernelScale;
        outIndex = outIndex + 1;
    end
end
    