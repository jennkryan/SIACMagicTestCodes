function out = evalKernelPP( bspline, rspline, deg, lambda, scale, offset, coeff, x )

out = x * 0;
for gamma = 0:rspline-1
    xnode = -0.5*(rspline-1) + gamma + lambda; 
    cIndex = gamma + 1;
    out = out + coeff(cIndex) * evalMappedPP( bspline, scale, offset - (xnode * scale), x );
end
