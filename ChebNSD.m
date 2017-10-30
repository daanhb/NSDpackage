function [ X, W ] = ChebNSD( a,b,freq,N,g,Dg )
%Where important information about the phase g(x) is unknown, this code
%will approximate the information using Chebfun.
%here g and it's derivative Dg are passed as anonymous functions


%ASSUMES STATIONARY POINTS LIE IN INTERIOR, AND NO SINGULARITIES
    ErrTol=1E-3;

%first determine stationary points
    ChebDg=chebfun(Dg,[a b], 'splitting', 'on');
    stationaryPoints=unique(roots(ChebDg));
%now determine the order of these points:
    stationaryPointsOrder = zeros(size(stationaryPoints));
    for j=1:length(stationaryPoints)
        ChebDng=ChebDg;
        %keep trying higher derivatives until non-zero at stationary point:
        while abs(ChebDng(stationaryPoints(j))) < ErrTol
            %consider a higher derivative:
            ChebDng=diff(ChebDng);
            stationaryPointsOrder(j) = stationaryPointsOrder(j) + 1;
        end
    end
    
    realPoints=[a repmat(stationaryPoints+1,1,2) b];
    
    pathPowers=[1 repmat(stationaryPointsOrder+1,1,2) 1];
    
    numPaths=2*length(stationaryPoints)+2;
    
    for SDpath=1:numPaths
        
        [x{SDpath}, w{SDpath}]=pathQuad( 0, inf, pathPowers(SDpath), [], N );
        Pmax=max(x{SDpath}./(freq^(1/pathPowers(SDpath))));
        
        %now approximate SD path, by solving DE with operator:
%         if mod(SDpath,2)==1 %move along SD path to infinity
            L{SDpath}=chebop(0, Pmax);
%         else
%             L{SDpath}=chebop(-Pmax,0);
%         end
        L{SDpath}.op = @(p,h) diff(h)-pathPowers(SDpath)*1i*p^(pathPowers(SDpath)-1)/Dg(h);
        L{SDpath}.lbc=realPoints(SDpath);
        h{SDpath}=L{SDpath}\0;
        dh_dp{SDpath}=diff(h{SDpath});
        
        W_{SDpath}=((-1)^(SDpath+1)/(freq^(1/pathPowers(SDpath))))*exp(1i*freq*g(realPoints(SDpath))).*dh_dp{SDpath}(x{SDpath}./(freq^(1/pathPowers(SDpath)))).*w{SDpath};
        X_{SDpath}=h{SDpath}(x{SDpath}./freq^(1/pathPowers(SDpath))); 
                
    end
    
    X=[];   W=[];
    for SDpath=1:numPaths
        W=[W; W_{SDpath}];
        X=[X; X_{SDpath};];          
    end

end