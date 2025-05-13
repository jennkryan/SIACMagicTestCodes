
function CC = SIACcoeff(moments,smoothness)

numcoeff = moments + 1; % number of function translates.
%Hscale = input('Input scaling parameter (H*dx, usually tied to mesh, default is 1):  ');
Hscale = 1;
dissipation = smoothness+2;
SIACcoeff = zeros(numcoeff,1);
Cmatrix = zeros(numcoeff);
Cmatrix(1,1:numcoeff) = 1;
xrepro = zeros(numcoeff,1); % polynomial reproduction.
xrepro(1,1) = 1;

for mrow=1:numcoeff-1
    for mcol = 1:numcoeff
        gamma = mcol-floor((numcoeff-1)/2)-1;
        if dissipation == 1
            Cmatrix(mrow+1,mcol) = ((gamma+0.5)^(mrow+1)-(gamma-0.5)^(mrow+1))/(mrow+1);
        else
            nsum = 0;
            for n=0:dissipation-1
                nsum = nsum + (-1)^(n+dissipation-1)*nchoosek(dissipation-1,n)*((2*n+2-dissipation+2*gamma)^(mrow+dissipation)-(2*n-dissipation+2*gamma)^(mrow+dissipation));
            end
            Cmatrix(mrow+1,mcol) = factorial(mrow)/factorial(mrow+dissipation)*0.5^(mrow+dissipation)*nsum;
        end
    end
end
SIACcoeff = Cmatrix\xrepro; %SIACcoeff are what is implemented for the physical space filter.

CC = SIACcoeff;
% numcoeff2 = floor((numcoeff-1)/2) +1; 
% CC = zeros(numcoeff2,1);
% for k=1:numcoeff2
%     CC(k) = SIACcoeff(numcoeff2+1-k);
% end
% CC are what is implemented for the Fourier SIAC filter.

return