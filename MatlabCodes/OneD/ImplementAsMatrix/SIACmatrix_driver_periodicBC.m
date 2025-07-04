
%Takes FD data (POINT VALUES) and calculates filtered solution
%In order to do so, it pretends it is a nodal DG solution with the same
%XNodes in each elemXent

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%!!!! COMMENT OUT THESE LINES IF READING IN DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Nx = 16;
% domain
xleft = -1;
xright = 1;
xlength = xright - xleft;
dx = (xright-xleft)/Nx;
xmesh = xleft:dx:xright; % elemXent interfaces -- row vector

% set approximation degree -- cannot choose if reading in data.
degX = 2;
DGorderX = degX+1;

fexact = @(x) sin(pi.*x);

DOFx = DGorderX*Nx;
dXall = xlength/DOFx;

% Local XNodes
%[XNodes,w]=JacobiGZW(DGorder,1.0,1.0);
XNodes = (-1+(2.*((1:DGorderX)-1))/(DGorderX-1))'; %column vector
% Global XNodes
Xall = xmesh(1:Nx)+(0.5*dx).*(XNodes+1); % all data points
Unode = fexact(Xall)+randn(size(Xall)).*dx^(DGorderX);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%!!!! UNCOMMENT THESE LINES IF READING IN DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Here is where you read in data (data locaction, data values), but we will create data.
% Normally, you would read in y, and the function.  From the function, you would have

% read in Xall = data location
% read in Unode = data values -- assumes columns are elemXents, rows are
% XNodes.
% DOFx = size(Xall);    % Total number of given data points
% dXall = Xall(2)-Xall(1); % distance between data points
% DGinfo = divisors(Ny);
% %numDiv = length(DGinfo);
% %possible divisors
% idof = ceil(length(DGinfo)/2); % number of possibilities
% % DGinfo(1:floor(length(DGinfo)/2)) are possible polynomial degrees;
% % DGinfo(ceil(length(DGinfo)/2):end) are possible number of elemXents
% % Pairs: p=DOF(1), Nx=DOF(end); p=DOF(2), Nx=DOF(end-1),...
% % Nx = DOF(numDiv+1-idof);
% % DGorder = DOF(idof);
% % THIS ASSUMPTION CAN BE REVISED.
% % we take the highest degree polynomial out of all such that the support of
% % the full kernel is less than the size of the domain.
% klength = 3*(DGinfo(1:ceil(0.5*length(DGinfo)));
% for pair = 1:floor(length(DGinfo)/2)
%      if klength(pair) < DGinfo(length(DGinfo)+1-pair)
%          NX = DGinfo(length(DGinfo)+1-pair)
%          DGorderX = DGinfo(pair)
%      else
%          break
%      end
% end

% xleft = Xall(1);    % Determine left boundary
% xright = Xall(end);   % Determine right boundary
% dx = (xright-xleft)/Nx; % mesh size
% xmesh = xleft:dx:xright;
% Assuming using same XNodes in every elemXent
% XNodes(1:DGorderX) = (2/dx).*(Xall(1:DGorderX)-Xall(1))-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Global nodal mesh matrix
xglobal = reshape(Xall,[DGorderX,Nx]); % reshape XNodes so that each column gives all the points associated to an elemXent.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%!!!! UNCOMMENT THESE LINES IF READING IN DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%degX = DGorderX-1; %uncomment this line if reading

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Apxbasis = zeros(DGorderX,DGorderX); %(modes,XNodes)
for m=1:DGorderX
    Apxbasis(m,:) = legendreP(m-1,XNodes);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No need to uncomment here -- THIS TO CHECK POST-PROCESSOR

% calculate modes via quadrature
% [zq,wq] = JacobiGZW(DGorder,0.0,0.0);
% qbase = zeros(DGorder,DGorder);
% qImass = zeros(DGorder,1);
% for m=1:DGorder
%     qbase(m,:) = legendreP(m-1,zq);
%     qImass(m) = 0.5*(2*m-1);
% end
% uhat = zeros(DGorder,Nx);
% for j = 1:Nx
%     xq = xmesh(j)+0.5*dx*(zq+1);
%     fq = fexact(xq);
%     fqw = fq.*wq;
%     for m=1:DGorder
%         uhat(m,j) = qbase(m,:)*fqw*qImass(m);
%     end
% end
% 
% % create nodal solution
% Unode = zeros(size(uhat));
% for j=1:Nx
%     for k=1:DGorder 
%         Unode(k,j) = Apxbasis(:,k)'*uhat(:,j);
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct Matrix for converting XNodes to modes
NtoMx = zeros(DGorderX);
for k=0:DGorderX-1 % Loop over Lagrange basis
    basis = @(x) prod(x-XNodes(1:k))./prod(XNodes(k+1)-XNodes(1:k))*prod(x-XNodes(k+2:DGorderX))./prod(XNodes(k+1)-XNodes(k+2:DGorderX));
    for m=0:DGorderX-1 % Loop over Legendre basis
        NtoMx(m+1,k+1) = integral(@(x) basis(x).*legendreP(m,x),-1,1,'ArrayValued',true)*(0.5*(2*m+1));
    end
end
uhat = zeros(DGorderX,Nx);
for j=1:Nx
    uhat(:,j) = NtoMx*Unode(:,j);
end


uapprox = zeros(DGorderX,Nx);
for j=1:Nx %elemXents
    for k=1:DGorderX %XNodes
        uapprox(k,j) = Apxbasis(:,k)'*uhat(:,j);
    end
end

subplot(1,2,1)
plot(xglobal(:),fexact(xglobal(:)),'k','LineWidth',2,'DisplayName','Target')
hold on
plot(xglobal(:),uapprox(:),'b--','LineWidth',2,'DisplayName','Input Data')

subplot(1,2,2)
semilogy(xglobal(:),abs(fexact(xglobal(:))- uapprox(:)),'b--','LineWidth',2,'DisplayName','Input Data')
hold on

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % STUFF FOR POSTPROCESSOR
% % kernel width is moments + BSorderX
% % full filter: smoothness = deg-1, moments = 2*deg'
filterparam = 0;
while strcmp(filterparam,'Y') == 0 && strcmp(filterparam,'N') == 0
    filterparam = input('Do you want to choose the filter parameters? (Y/N)  ',"s");
    if strcmp(filterparam,'Y') == 1
        Xsmoothness = input('Input smoothness/dissipation/number of continuous derivatives (>= -1): ');
        Xmoments = input('Input number of moments to satisfy (>=0, should be even. 2 moments preserves mean and variance):  ');
    elseif strcmp(filterparam,'N') == 1
        Xsmoothness = degX-1;
        Xmoments = 2*degX; 
    else
        disp(['ERROR: Please only answer Y if yes or N if no.'])
    end
end
% % Ensure an even number of moments.
if mod(Xmoments,2) ~= 0
    Xmoments = Xmoments + 1;
end
X = Xmoments + 1;
BSorderX = Xsmoothness + 2;
BSknotsX = (-BSorderX/2:BSorderX/2)';
BSSupportX = [floor(BSknotsX(1)),ceil(BSknotsX(BSorderX+1))];
BSlenX = (BSSupportX(2)-BSSupportX(1))+1;




% 
% % Grab the symmetric SIAC weights
CCX = SIACcoeff(Xmoments,Xsmoothness);

% disp(['Using ',num2str(Xmoments+1),' B-Splines of order ',num2str(Xsmoothness+2),'.'])
% disp('SIAC Kernel coefficients')
% disp(CCX)
% 
% % Calculate integrals of B-Splines
[BSIntX] = GrabIntegrals(Xsmoothness, XNodes);

kernellengthX = 2*ceil((Xmoments + BSorderX)/2) + 1;
halfkerX = ceil((Xmoments + BSorderX)/2);

% % % % SIAC eveluated at FD points
SIACmatrixX = zeros(DGorderX,kernellengthX,DGorderX);
% % % % Post-processing matrix
for k = 1:DGorderX % loop over XNodes
    for igam = 1:Xmoments+1 % basis functions
        SIACmatrixX(:,igam:igam+BSlenX-1,k) = SIACmatrixX(:,igam:igam+BSlenX-1,k)+CCX(igam).*BSIntX(:,1:BSlenX,k);
    end
end
% 
% % Once we have the SIACmatrixX, we can apply it to our modal coefficients:
ustar = zeros(DGorderX,Nx);
for elemX = 1:Nx % assume we don't have a periodic solution so that we only post-process the interior elemXents
    if elemX < halfkerX+1
        RIdx = elemX + halfkerX;
        LIdx = kernellengthX-RIdx;
        VX(:,(kernellengthX+1-RIdx):kernellengthX) = uhat(:,1:RIdx);
        VX(:,1:LIdx) = uhat(:,(Nx+1-LIdx):Nx);
    elseif elemX > Nx-halfkerX
        RIdx = Nx+1-(elemX-halfkerX);
        LIdx = kernellengthX-RIdx;
        VX(:,1:RIdx) = uhat(:,(elemX-halfkerX):Nx);
        VX(:,RIdx+1:kernellengthX) =  uhat(:,1:LIdx);
    else
        VX = uhat(:,elemX-halfkerX:elemX+halfkerX);
    end
    for k = 1:DGorderX
        ustar(k,elemX) = sum(SIACmatrixX(:,:,k).*VX,'all');
    end
end

subplot(1,2,1)
hold on
plot(xglobal(:),ustar(:),':','LineWidth',3,'DisplayName','Filtered data')
legend
xlabel('x',FontSize=12)
ylabel('Data values',FontSize=12)
title('Input and Filtered data',FontSize=12,FontWeight='bold')

subplot(1,2,2)
semilogy(xglobal(:),abs(fexact(xglobal(:)) - ustar(:)),':','LineWidth',3,'DisplayName','Filtered data')
legend

xlabel('x',FontSize=12)
ylabel('log(|Error compared to exact|)',FontSize=12)
title('log(|Error|)',FontSize=12,FontWeight='bold')

disp(['Using N = ',num2str(Nx),' elemXents of order ',num2str(DGorderX)])
disp(['DG error L-inf error:        ',num2str(max(abs(fexact(xglobal(:))- uapprox(:))))])
disp(['Filtered error L-inf error:  ',num2str(max(abs(fexact(xglobal(:)) - ustar(:))))])

figure(2)
plot(xglobal,uapprox-ustar,'LineWidth',2)
xlabel('X')
ylabel('log(|Data - Filtered Data|)',FontSize=12)
title('Difference: Data - Filtered Data',FontSize=12,FontWeight='bold')