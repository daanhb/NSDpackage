function [h, dhdp, w] = nearlyFinitePathFix(nfParamPoint, nfCritPoint, criticalStartPoints, pathPower, IC, G, freq, N, COVsingularities, RelTol)

    Pend=nfParamPoint*2^(-1/pathPower); %truncate the path far from the critical point

    %make the first bit of the path:
    [p, Wstart] = pathQuadFinite( Pend, COVsingularities, freq, N );
    [~, Hstart] = ode45(@(t,y) NSDpathODE(t,y,pathPower-1,G, IC, false), p, IC, odeset('RelTol',RelTol) );      
    
    %now go from here to the nearby critical point, in a straight line:
    L=abs(Hstart(end,1)-nfCritPoint);
    [zEnd, Wend] = quad_gauss(N, 0, L);
    Hend = Hstart(end,1) + zEnd*(nfCritPoint-Hstart(end,1))/L;
    [pEnd, Wend] = pathQuadFinite( L, COVsingularities, freq, N );
    Hend = Hstart(end,1) + pEnd*(nfCritPoint-Hstart(end,1))/L;
    %piece this wonky path together:
    h = [Hstart(:,1); Hend];
    w = [Wstart; Wend];%.*exp(1i*freq*G{1}(Hend))
    dhdp = [Hstart(:,2); ones(length(Wend),1)];
    
end