function coeff = getKernelCoeffsMatrix( numBS, bSplineDegree )

Amat = zeros(numBS,numBS);
bvec = zeros(1,numBS);
bvec(1) = 1;

BSorder = bSplineDegree+1;
BSint = zeros(1,numBS);
bspline = getBSplinePP(BSorder );

for n=0:numBS-1
    sumint = 0;
    for j=0:bSplineDegree
        sumint = sumint + (-1)^(bSplineDegree+j)*nchoosek(bSplineDegree,j) ...
            *((1-0.5*BSorder+j)^(n+BSorder)-(-0.5*BSorder+j)^(n+BSorder));
    end
    BSint(n+1) = factorial(n)/factorial(n+BSorder)*sumint;
end
    
RS = ceil(0.5*(numBS-1));

for m=0:numBS-1
    for gam = -RS:RS
        sumcoeff = 0;
        for n = 0:m
            sumcoeff = sumcoeff + (-1)^m*nchoosek(m,n)*gam^(m-n)*BSint(n+1);
        end
        Amat(m+1,gam+RS+1) = (-1)^m*sumcoeff;
    end
end

AImat = inv(Amat);
coeff = AImat*bvec';

end