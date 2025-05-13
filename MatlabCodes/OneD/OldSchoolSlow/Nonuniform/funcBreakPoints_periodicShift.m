function breaks = funcBreakPoints_periodicShift( boundaries, support, shift )
% Returns a list (in world coords) of the periodic breakpoints of the
% function as defined by the boundaries and the region of support
% doesn't include the actual support limits -- strictly internal breaks
%
% assumes that the size of the support is less than that of the domain

    breaks = [];
    
    Nelem = size(boundaries, 2) - 2*shift;
    lower = support(1);
    upper = support(2);
    supportWidth = upper - lower;
    
    domainMin = boundaries(1 + shift);
    domainMax = boundaries( size( boundaries, 2) - shift );
    domainWidth = domainMax - domainMin;
    
    % this works because I'm assuming that the support will always fit into the domain
    doubleBoundary = [ boundaries(shift+1:shift+Nelem+1) boundaries(2:size(boundaries,2)) + domainWidth ];
    
    % shift support's lower end into the domain
    offset = 0;     % the number of support widths I had to shift the support to get it to start in the domain
    while ( lower < domainMin )
        offset = offset + 1;
        lower = lower + domainWidth ;
        upper = upper + domainWidth ;
    end
    while ( lower > domainMax )
        offset = offset - 1;
        lower = lower - domainWidth ;
        upper = upper - domainWidth ;
    end
    
    boundaryIndex = 1;
    while ( doubleBoundary( boundaryIndex ) <= lower )
        boundaryIndex = boundaryIndex + 1;
    end
    
    brkIndex = 1;
    while ( doubleBoundary( boundaryIndex ) < upper )
        breaks( brkIndex ) = doubleBoundary( boundaryIndex );
        boundaryIndex = boundaryIndex + 1;
        brkIndex = brkIndex + 1;
    end
    
    breaks = breaks - domainWidth * offset;
    