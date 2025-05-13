function [BSInt] = GrabIntegrals(smoothness,Nodes)
%GrabIntegrals(moments,smoothness,Nodes)
% calculates integrals of the B-splines 
% if BSorder is even, then there are 2*BSorder of them;
% otherwise, if BSorder is odd, there are (2*BSorder+1).
% BSInt(cell,deg,node)
% Assumes data converted into Legendre basis

    
    % stuff for DG polynomial
    DGorder = length(Nodes);
    degDG = DGorder -1;

    % stuff for B-spline
    BSorder = smoothness + 2;
    degBS = BSorder -1;
    BSknots = (-BSorder/2:BSorder/2)';
    BSsupport = [floor(BSknots(1)),ceil(BSknots(BSorder+1))];
    BSlen = (BSsupport(2)-BSsupport(1))+1;

    % Calculate Bspline integrals
    % BInt(mode basis, cell, node)
    BIntL=zeros(DGorder,BSlen,DGorder);
    BIntR=zeros(DGorder,BSlen,DGorder);
    for k=1:DGorder % Loop over evaluation points.  In this case, they are the NODES.
        zeta = Nodes(k);
        xicell = zeta-sign(zeta)*mod(BSorder,2);
        for i=BSsupport(1):BSsupport(2) %Loop over the elements in the support of the B-Spline
            j = i-BSsupport(1)+1; % set array index
            for m=1:DGorder %Loop over basis functions -- Legendre Polynomials
                BIntL(m,j,k) = 0.5*integral(@(x) BStest(BSorder,BSknots,0.5*(zeta-x)-i)*legendreP(m-1,x),-1,xicell,'ArrayValued',true);
                BIntR(m,j,k) = 0.5*integral(@(x) BStest(BSorder,BSknots,0.5*(zeta-x)-i)*legendreP(m-1,x),xicell,1,'ArrayValued',true); 
            end
        end
    end

    BSInt=zeros(DGorder,BSlen,DGorder);
    for k=1:DGorder
        BSInt(:,:,k) = BIntL(:,:,k)+BIntR(:,:,k);
    end

end % function end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%