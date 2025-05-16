
%Takes FD data (POINT VALUES) and calculates filtered solution
%In order to do so, it pretends it is a nodal DG solution with the same
%nodes in each element

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%!!!! COMMENT OUT THESE LINES IF READING IN DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Nx = 16, Ny = 10; 
% domain
xleft = -1;
xright = 1;
xlength = xright - xleft;
dx = (xright-xleft)/Nx;
xmesh = xleft:dx:xright; % element interfaces -- row vector

% set approximation degree -- cannot choose if reading in data.
deg = 3;
DGorder = deg+1;

fexact = @(x) sin(pi.*x);

Ny = DGorder*Nx;
dy = dx/DGorder;
% Local nodes
%[Nodes,w]=JacobiGZW(DGorder,1.0,1.0);
Nodes = (-1+(2.*((1:DGorder)-1))/(DGorder-1))'; %column vector
% Global nodes
y = xmesh(1:Nx)+(0.5*dx).*(Nodes+1); % all data points
Unode = fexact(y)+randn(size(y)).*dx^(DGorder);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%!!!! UNCOMMENT THESE LINES IF READING IN DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Here is where you read in data (data locaction, data values), but we will create data.
% Normally, you would read in y, and the function.  From the function, you would have

% read in y = data location
% read in Unode = data values -- assumes columns are elements, rows are
% nodes.
% Ny = size(y);    % Total number of given data points
% dy = y(2)-y(1); % distance between data points
% DOF = divisors(Ny);
% %numDiv = length(DOF);
% %possible divisors
% idof = ceil(length(DOF)/2); % number of possibilities
% % DOF(1:floor(length(DOF)/2)) are possible polynomial degrees;
% % DOF(ceil(length(DOF)/2):end) are possible number of elements
% % Pairs: p=DOF(1), Nx=DOF(end); p=DOF(2), Nx=DOF(end-1),...
% % Nx = DOF(numDiv+1-idof);
% % DGorder = DOF(idof);
% % THIS ASSUMPTION CAN BE REVISED.
% % we take the highest degree polynomial out of all such that the support of
% % the full kernel is less than the size of the domain.
% klength = 3*(DOF(1:ceil(0.5*length(DOF)));
% for pair = 1:floor(length(DOF)/2)
%      if klength(pair) < DOF(length(DOF)+1-pair)
%          NX = DOF(length(DOF)+1-pair)
%          DGorder = DOF(pair)
%      else
%          break
%      end
% end

% xleft = y(1);    % Determine left boundary
% xright = y(end);   % Determine right boundary
% dx = (xright-xleft)/Nx; % mesh size
% xmesh = xleft:dx:xright;
% Assuming using same nodes in every element
% Nodes(1:DGorder) = (2/dx).*(y(1:DGorder)-y(1))-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Global nodal mesh matrix
xglobal = reshape(y,[deg+1,Nx]); % reshape nodes so that each column gives all the points associated to an element.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%!!!! UNCOMMENT THESE LINES IF READING IN DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%deg = DGorder-1; %uncomment this line if reading

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Apxbasis = zeros(DGorder,DGorder); %(modes,nodes)
for m=1:DGorder
    Apxbasis(m,:) = legendreP(m-1,Nodes);
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

% Construct Matrix for converting nodes to modes
NtoM = zeros(DGorder);
for k=0:DGorder-1 % Loop over Lagrange basis
    basis = @(x) prod(x-Nodes(1:k))./prod(Nodes(k+1)-Nodes(1:k))*prod(x-Nodes(k+2:DGorder))./prod(Nodes(k+1)-Nodes(k+2:DGorder));
    for m=0:DGorder-1 % Loop over Legendre basis
        NtoM(m+1,k+1) = integral(@(x) basis(x).*legendreP(m,x),-1,1,'ArrayValued',true)*(0.5*(2*m+1));
    end
end
uhat = zeros(DGorder,Nx);
for j=1:Nx
    uhat(:,j) = NtoM*Unode(:,j);
end


uapprox = zeros(DGorder,Nx);
for j=1:Nx %elements
    for k=1:DGorder %nodes
        uapprox(k,j) = Apxbasis(:,k)'*uhat(:,j);
    end
end

subplot(1,2,1)
plot(y(:),fexact(y(:)),'k','LineWidth',2,'DisplayName','Target')
hold on
% plot(y(:),Unode(:),'--','LineWidth',2,'DisplayName','Given data')
% hold on
plot(y(:),uapprox(:),'b--','LineWidth',2,'DisplayName','Input Data')

subplot(1,2,2)
% plot(y(:),fexact(y(:))- Unode(:),'--','LineWidth',2,'DisplayName','Given data')
% hold on
plot(y(:),fexact(y(:))- uapprox(:),'b--','LineWidth',2,'DisplayName','Input Data')
hold on

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % STUFF FOR POSTPROCESSOR
% % kernel width is moments + BSorder
% 

% % full filter: smoothness = deg-1, moments = 2*deg'
smoothness = 2;%deg-1; %input('Filter properties: Input desired smoothness):  ');
moments = 2*deg; %input('Filter properties: Input number of moments:  ');
numsplines = moments + 1;
BSorder = smoothness + 2;
BSknots = (-BSorder/2:BSorder/2)';
BSsupport = [floor(BSknots(1)),ceil(BSknots(BSorder+1))];
BSlen = (BSsupport(2)-BSsupport(1))+1;



% % Ensure an even number of moments.
% % if mod(moments,2) ~= 0
% %     moments = moments + 1;
% % end
% 
% % Grab the symmetric SIAC weights
CC = SIACcoeff(moments,smoothness);

% disp(['Using ',num2str(moments+1),' B-Splines of order ',num2str(smoothness+2),'.'])
% disp('SIAC Kernel coefficients')
% disp(CC)
% 
% % Calculate integrals of B-Splines
[BSInt] = GrabIntegrals(smoothness, Nodes);

kernellength = 2*ceil((moments + BSorder)/2) + 1;
halfker = ceil((moments + BSorder)/2);
% 
% % need to figure out how to efficiently sum integrals.
% % In elements i = j-(n+r)/2:j have sum_gamma=1:i+(n+r)/2+1
% % c_gamma*I(i+(n+r)/2+2-gamma)
% % % % SIAC eveluated at FD points
SIACmatrix = zeros(DGorder,kernellength,DGorder);
% % % % Post-processing matrix
for k = 1:DGorder % loop over nodes
    for igam = 1:moments+1
        SIACmatrix(:,igam:igam+BSlen-1,k) = SIACmatrix(:,igam:igam+BSlen-1,k)+CC(igam).*BSInt(:,1:BSlen,k);
    end
%     disp(['Node ',num2str(k),':   SIACmatrix'])
%     disp(SIACmatrix)
end
% 
% % Once we have the SIACmatrix, we can apply it to our modal coefficients:
ustar = zeros(DGorder,Nx-(2*halfker));
for elem = (halfker+1):(Nx-(halfker)) % assume we don't have a periodic solution so that we only post-process the interior elements
    for k=1:DGorder %loop over nodes
        V = uhat(:,elem-halfker:elem+halfker);
        ustar(k,elem-halfker) = sum(SIACmatrix(:,:,k).*V,'all');
        if elem == (halfker+1) & k==1
            disp(size(V))
            disp(size(SIACmatrix(:,:,k).*V))
        end
        % In 2D ustar(kx_node,ky_node,j_x,j_y) =
        % sum(SIACmatrix(:,:,kx).*V(kx,:,j_x,j_y),'all')*sum(SIACmatrix(:,:,ky).*V(:,ky,j_x,j_y),'all')
    end
end
z = y(:,(halfker+1):(Nx-(halfker)));

subplot(1,2,1)
% plot(y(:),fexact(y(:)),'k','LineWidth',2,'DisplayName','Exact')
% hold on
% plot(y(:),Unode(:),'--','LineWidth',2,'DisplayName','Given data')
hold on
plot(z(:),ustar(:),':','LineWidth',3,'DisplayName','Filtered data')
legend
xlabel('x',FontSize=12)
ylabel('Data values',FontSize=12)
title('Input and Filtered data',FontSize=12,FontWeight='bold')

subplot(1,2,2)
% plot(y(:),fexact(y(:))- Unode(:),'--','LineWidth',2,'DisplayName','Given data')
% hold on
% plot(y(:),fexact(y(:))- uapprox(:),'b--','LineWidth',2,'DisplayName','Approximation')
% hold on
plot(z(:),fexact(z(:)) - ustar(:),':','LineWidth',3,'DisplayName','Filtered data')
legend

xlabel('x',FontSize=12)
ylabel('log(|Error compared to exact|)',FontSize=12)
title('log(|Error|)',FontSize=12,FontWeight='bold')

disp(['Using N = ',num2str(Nx),' elements of order ',num2str(DGorder)])
disp(['DG error L-inf error:        ',num2str(max(abs(fexact(y(:))- uapprox(:))))])
disp(['Filtered error L-inf error:  ',num2str(max(abs(fexact(z(:)) - ustar(:))))])

figure(2)
plot(z,uapprox(:,halfker+1:Nx-halfker)-ustar,'LineWidth',2)
xlabel('X')
ylabel('log(|Data - Filtered Data|)',FontSize=12)
title('Difference: Data - Filtered Data',FontSize=12,FontWeight='bold')