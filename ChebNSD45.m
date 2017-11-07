function [ X, W ] = ChebNSD45( a,b,freq,N,g)
%Where important information about the phase g(x) is unknown, this code
%will approximate the information using Chebfun.


%ASSUMES STATIONARY POINTS LIE IN INTERIOR, AND NO SINGULARITIES
    ErrTol=1E-3;
    RelTol=1E-3;

%first determine stationary points
    G{1}=chebfun(g,[a b], 'splitting', 'on');
    G{2}=diff(G{1});
    stationaryPoints=unique(roots(G{2}));
%now determine the order of these points:
    stationaryPointsOrder = zeros(size(stationaryPoints));
    for j=1:length(stationaryPoints)
        n=2;
        while abs(G{n}(stationaryPoints(j))) < ErrTol
            %consider a higher derivative:
            n=n+1;
            G{n}=diff(G{n-1});
            stationaryPointsOrder(j) = stationaryPointsOrder(j) + 1;
        end
    end
    
    realPoints=[a repmat(stationaryPoints,1,2) b];
    
    pathPowers=[1 repmat(stationaryPointsOrder+1,1,2) 1];
    
    numPaths=2*length(stationaryPoints)+2;
    
    

    %now choose the right one:
    
    
    for SDpath=1:numPaths
        
        
        [x{SDpath}, w{SDpath}]=pathQuad( 0, inf, pathPowers(SDpath), [], N );
        P0=[0; (x{SDpath}./(freq^(1/pathPowers(SDpath))))];
        

        [~,H] = ode45(@(t,y) NSDpathODE(t,y,pathPowers(SDpath)-1,G), P0, [NSDpathIC( pathPowers(SDpath),SDpath); realPoints(SDpath) ], odeset('RelTol',RelTol) );
        %H(:,1) contains h'(p), H(:,2) contains h(p), will throw away initial
        %point
        
      %  h{SDpath}=H(2:end,1); dhdp{SDpath}=H(2:end,2);
        
        h{SDpath}=H(2:end,1);
        if pathPowers(SDpath)>1
            dhdp{SDpath}=H(2:end,2);
        else
            %dhdpCheb=diff(chebfun(h{SDpath}));
            %error comes from here, as Chebfun doesn't know x coords of
            %points
            %*a better option is to re-insert h(p) into DE for h'(p)*
            dhdp{SDpath}=1i./G{2}(H(2:end));
        end
                
        W_{SDpath}=((-1)^(SDpath+1)/(freq^(1/pathPowers(SDpath))))*exp(1i*freq*g(realPoints(SDpath))).*dhdp{SDpath}.*w{SDpath};
        %ERROR HERE. Somehow for g(x)~x case, the i/omega is already
        %absorbed into h'(p) (causing an error), but for g(x)~x^2 it isn't 
        X_{SDpath}=h{SDpath}; 
                
    end
    
    X=[];   W=[];
    for SDpath=1:numPaths
        W=[W; W_{SDpath}];
        X=[X; X_{SDpath};];          
    end

end