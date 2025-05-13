clear all 
close all

%func = @one;
%func = @linear;
%func = @quadratic;
%func = @cubic;
func = @oneSine;
%func = @(x) exp(-200*(x-.5).^2);
%func = @(x) abs(cos((pi).*x+10^(-2)));


% this WILL break if you don't have enough elements to support the degree
% of the kernel.  If you receive an index out of bounds error, increase the
% number of elements
N = 20%input('Input number of elements:  ');  %number of elements
M = 3%input('Input polynomial degree (NOT ORDER):  ');   %polynomial degree

%evalPoints = input('Input number of evaluation points per element:  ');
evalPoints = 2*(M+1);

[z,w] = JacobiGZW(M+1,0.0,0.0);  %get quadrature points and weights
% six quadrature points at which to evaluate the funcitons -- non
% superconvergent
%[zEval, wEval] = JacobiGZW(evalPoints, 0.0, 0.0);
zEval = zeros(2*(M+1),1);
wEval = zeros(2*(M+1),1);
zEval(1:M+1) = 0.5.*(z-1);
zEval(M+2:2*(M+1)) = 0.5.*(z+1);
NumP = size(zEval, 1);


for i = 0:M
  Lvalue(i+1,:) = JacobiPoly(i,z,0.0,0.0)';
  %LvalueMass(i+1) = dot(w,Lvalue(i+1,:).*Lvalue(i+1,:));
end
for i = 0:M
    Lvalue(i+1,:) = sqrt(i+0.5)*Lvalue(i+1,:);
end
    

%!!!!!!!!!!!!!!! IMPORTANT!!!!!!!!!!!!!!!!!!!!!!
%ell=M+1; % order of B-spline for filtering  For superconvergence should be M+1
    ellp2 = M-1%input('Input smoothness required (Usually = degree - 1.  For this implementation, usually >=0).  0 = continuous:  ');
    ell = ellp2+2;

%!!!!!!!!!!!!!!! IMPORTANT!!!!!!!!!!!!!!!!!!!!!!
RS = 2*M+1;%input('Number B-splines/moments+1 (0 = automatic):  ');

MeshMax = 4;

for mesh = 1:MeshMax

    %%%%% UNCOMMENT THESE TWO LINES FOR UNIFORM INTERVALS
    h = 1.0/N;
    x = 0:h:1;
    %%%%% UNCOMMENT THIS SINGLE LINE FOR NON_UNIFORM INTERVAL
    % x = nonUniformIntervals( N, 1.0, 0 ); %element definitions

    for nel = 1:N
        h = x( nel + 1 ) - x( nel );
        zmap = 0.5*h*(z+1.0) + x(nel);  %affine mapping of [-1,1] quadrature points on to element
        f = feval( func, zmap);           %function to be mapped evaluated at mapped quadrature points
        for i = 0:M
            %uhat(nel,i+1) =  dot(w,f'.*Lvalue(i+1,:))/LvalueMass(i+1);
            uhat(nel,i+1) =  dot(w,f'.*Lvalue(i+1,:));
        end
    end

    save('modes.txt','uhat','-ascii')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PLOT TRUE DATA AND APPROXIMATION DATA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    count = 1;
    for nel = 1:N
        h = x( nel + 1 ) - x( nel );
        zmap = 0.5*h*(zEval'+1.0) + x(nel);
        f = feval( func, zmap);            %function to be mapped evaluated at mapped plotting points
        xvec(count:count+NumP-1) = zmap;
        ftrue(count:count+NumP-1) = f;
    
        %reconstruct polynomial as linear combination of basis functions
        utemp = zeros(size(zEval));
        for i = 0:M
            utemp = utemp + uhat(nel,i+1)*sqrt(i+0.5)*JacobiPoly(i, zEval, 0.0 ,0.0); %Lvalue(i+1,:)';
        end
        fapprox(count:count+NumP-1) = utemp;
        count = count + NumP;
    end

    save('fapprox.txt','fapprox','-ascii');
% plot 
    error = max(abs(ftrue-fapprox));
    %close all; 
    figure(1)
%    subplot(1,2,1)
    if mesh == MeshMax
        plot(xvec,ftrue,'-b','LineWidth',2); 
    end
    hold on
    grid on
    plot(xvec,fapprox,'LineStyle','-.','LineWidth',3);hold on;
%     intPoints = feval( func, x );
%     plot(x, intPoints, 'o');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% POST PROCESS AND PLOT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %tic;        % time it
    % set kernelScale = 0 in order to use the largest element.
    kernelScale = 0;
    %kernelScale = 2*h;
    %kernelScale = 0.5*h;
    %kernelScale = h/(M+1.);
    filtered = nonUniformPost(uhat, x, ell, RS, kernelScale );
    %toc;

    figure(1)
    plot(filtered(1,:), filtered(2,:),'LineStyle',':','LineWidth',3);
    if mesh == MeshMax
        plot(xvec,ftrue,'-b'); 
    end
    hold on
    grid on

    

    fapproxError = abs(ftrue-fapprox)+ 1.0e-14 ;
    filterError = abs(ftrue - filtered(2,:))+ 1.0e-14 ;

    figure(2);
    subplot(1,2,1)
    semilogy(xvec, fapproxError,'LineWidth',2);
    if mesh == MeshMax
        ylabel('log_1_0( error )','FontSize',12,'FontWeight','bold');
        xlabel('x, world space','FontSize',12,'FontWeight','bold');
        title('Error of approximation','FontSize',14);
    end
    hold on
    subplot(1,2,2)
    semilogy(xvec, filterError,'LineWidth',2);
    if mesh == MeshMax
        ylabel('log_1_0( error )','FontSize',12,'FontWeight','bold');
        xlabel('x, world space','FontSize',12,'FontWeight','bold');
        title('Error of post-processed approximation','FontSize',14);
    end
    hold on

    ppapprox = filtered(2,:);
    save('ppapprox.txt','ppapprox','-ascii')
    
    DGerror = max(abs(fapproxError));
    SIACerr = max(abs(filterError));
    if mesh ==1
        s1 = sprintf('p=%2i,   ell = %2i\n',M,ell);
        s2 = sprintf('      N     DG L-Inf   Order        SIAC L-inf  Order\n');
        s3 = sprintf('    %3i %13.3d           %13.3d\n',N,DGerror,SIACerr);
        disp([s1,s2,s3])
    else
        DGorder = log(OldDGerr/DGerror)/log(2);
        SIACorder = log(OldSIACerr/SIACerr)/log(2);
        s4 = sprintf('    %3i %13.3d  %f %13.3d   %f\n',N,DGerror,DGorder,SIACerr,SIACorder);
        disp([s4])
    end
    
    N = 2*N;
    OldDGerr = DGerror;
    OldSIACerr = SIACerr;
end