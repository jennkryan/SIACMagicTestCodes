
function [out] = KernelCoeff(numBS,orderBS)

% grabs the kernel coefficients for the zero-th order B-spline

 A=zeros(numBS,numBS);
 A(1,1:numBS) = 1;
 for mrow = 1:numBS-1
     for mcol = 1:numBS
         gam = mcol - floor((numBS-1)/2)-1;
         if orderBS == 1
             A(mrow+1,mcol) = ((gam+0.5)^(mrow+1)-(gam-0.5)^(mrow+1))/(mrow+1);
         else
             jsum = 0;
             for j = 0:orderBS-1
                    jsum = jsum + (-1)^(j+orderBS-1)*nchoosek(orderBS-1,j)*((2*j+2-orderBS+2*gam)^(mrow+orderBS)-(2*j-orderBS+2*gam)^(mrow+orderBS));
             end
             A(mrow+1,mcol) = factorial(mrow)/factorial(mrow+orderBS)*0.5^(mrow+orderBS)*jsum;
         end
     end
 end

b = zeros(numBS,1);
b(1)=1.0;

%call the lu_factor function LU = linalg.lu_factor(A)
kcoeff = A\b;

out = kcoeff;

end

    
    
    