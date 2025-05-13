%JK Ryan
%January 2025

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Bpsi] = BStest(n, knots,x)

    
Bomega = @(x,ell) (x-knots(1:n-ell+2))./(knots(ell:n+1)-knots(1:n-ell+2));
for m=1:n
    if m==1
        if n==1
            B1 = @(x) heaviside(x-knots(1))-heaviside(x-knots(2));
        else
            B1 = @(x) heaviside(x-knots(1:n))-heaviside(x-knots(2:n+1));
        end
        Bfnew = B1(x);
    else
        if m==2
            weights = Bomega(x,m);
            Bm = weights(1:n-1).*Bfnew(1:n-1)+(1-weights(2:n)).*Bfnew(2:n);
        else
            weights = Bomega(x,m);
            Bm = weights(1:n-m+1).*Bm(1:n-m+1)+(1-weights(2:n+2-m)).*Bm(2:n+2-m);
        end
        Bfnew = Bm;
    end
end
Bpsi = Bfnew;
    
end
