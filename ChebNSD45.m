function [ X, W ] = ChebNSD45( a,b,freq,N,g)
%Where important information about the phase g(x) is unknown, this code
%will approximate the information using Chebfun and ODE45

    %tolerance for finding stationary points. Higher values make it more
    %likely to correctly determine their order, but less acurate locations
    ZeroThresh=1E-14;
    %tolerance for ODE45 solver:
    RelTol=1E-13;
    
    %get stationary points, their order, and higher derivatives of g if
    %required
    [ stationaryPoints, stationaryPointsOrder, G ] = getStationaryPoints( a,b,g,ZeroThresh );
    
    %each path will start/end at one of the following points:
    realPoints=[a repmat(stationaryPoints,1,2) b];
    
    %the polynomial degree at the start of the path:
    pathPowers=[1 repmat(stationaryPointsOrder+1,1,2) 1];
    
    %total number of paths:
    numPaths=2*length(stationaryPoints)+2;
    
    %now loop over all paths:
    for SDpath=1:numPaths
        
        %get weights and nodes
        [x{SDpath}, w{SDpath}]=pathQuad( 0, inf, pathPowers(SDpath), [], N );
        P0=[0; (x{SDpath}./(freq^(1/pathPowers(SDpath))))];

        [p,H] = ode45(@(t,y) NSDpathODE(t,y,pathPowers(SDpath)-1,G), P0, [NSDpathIC( pathPowers(SDpath),SDpath); realPoints(SDpath) ], odeset('RelTol',RelTol) );
        %H(:,1) contains h'(p), H(:,2) contains h(p), will throw away initial
        %points
        
        if pathPowers(SDpath)>1
            % h'(p) given as output of ODE45, so use it
            h{SDpath}=H(2:end,2);
            dhdp{SDpath}=H(2:end,1);
        else
            %re-insert h(p) into DE for h'(p)*
            h{SDpath}=H(2:end);
            dhdp{SDpath}=1i./G{2}(H(2:end));
        end
                
        %absorb h'(p) and other constants into weights, also negate paths
        %incoming from infinity
        W_{SDpath}=((-1)^(SDpath+1)/(freq^(1/pathPowers(SDpath))))*exp(1i*freq*g(realPoints(SDpath))).*dhdp{SDpath}.*w{SDpath};
        
        X_{SDpath}=h{SDpath}; 
                
    end
    
    X=[];   W=[];
    for SDpath=1:numPaths
        W=[W; W_{SDpath}];
        X=[X; X_{SDpath};];          
    end

end