function [ X, W ] = ChebNSD( a,b,freq,N,varargin )
%Where important information about the phase g(x) is unknown, this code
%will approximate the information using Chebfun.
%here g and it's derivative Dg are passed as anonymous functions
    g=varargin;

%ASSUMES STATIONARY POINTS LIE IN INTERIOR, AND NO SINGULARITIES
    ErrTol=1E-3;
    
%WARNING TO SELF:
%Chebfun is designed for functions on [a,b]=[-1,1], when this changes, will
%need to make other changes too
if a~=-1 || b~=1
    error('Need to account for Chebfun outside of [-1,1]');
end

%first determine stationary points
    ChebDg=chebfun(g{2},[a b], 'splitting', 'on');
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
    
    realPoints=[a repmat(stationaryPoints,1,2) b];
    
    pathPowers=[1 repmat(stationaryPointsOrder+1,1,2) 1];
    
    numPaths=2*length(stationaryPoints)+2;
    
    

    %now choose the right one:
    
    
    for SDpath=1:numPaths
        
        if mod(SDpath,2)==0
            PathOdd=false;
        else
            PathOdd=true;
        end
        
        %determine locally the branches, one of which will be the start of
        %the path:
        LocPaths=exp(1i*pi/(2*pathPowers(SDpath))*(4*(0:pathPowers(SDpath))+1));
    
        switch pathPowers(SDpath)
            case 1
                LocDh(SDpath)=1i;
            case 2
                if PathOdd
                    LocDh(SDpath)=LocPaths(1);
                else
                    LocDh(SDpath)=LocPaths(2);
                end
            case 3
                if PathOdd
                    LocDh(SDpath)=LocPaths(1);
                else
                    LocDh(SDpath)=LocPaths(2);
                end
            case 4
                error('Not sure about this 4th order stationary points yet...');
        end
        
        [x{SDpath}, w{SDpath}]=pathQuad( 0, inf, pathPowers(SDpath), [], N );
        Pmax=max(x{SDpath}./(freq^(1/pathPowers(SDpath))));
        
        %now approximate SD path, by solving DE with operator:
%         if mod(SDpath,2)==1 %move along SD path to infinity
              L{SDpath}=chebop(0, Pmax);
              %chebfun can only solve on bounded domains.
              %also complains when coefficients vanish at infinity
%         else
%             L{SDpath}=chebop(-Pmax,0);
%         end
    
        %L{SDpath}.op = @(p,h) diff(h)-pathPowers(SDpath)*1i*p^(pathPowers(SDpath)-1)/Dg(h);
        
        switch pathPowers(SDpath)
            case 1
                L{SDpath}.op = @(p,h) diff(h).*g{2}(h)-1i;
                L{SDpath}.lbc=[realPoints(SDpath)];
            case 2
                L{SDpath}.op = @(p,h) diff(h,2)*g{2}(h)+(diff(h))^2*g{3}(h)-2i;
                L{SDpath}.lbc=[realPoints(SDpath), LocDh(SDpath)];
        end
        
        %default Chebfun IVP solver fails due to 'vanishing coefficients',
        % use this one instead:
        cheboppref.setDefaults('ivpSolver', @chebcolloc2)
        %(there are probably others...)
        
        h{SDpath}=L{SDpath}\0;
        dh_dp{SDpath}=diff(h{SDpath});
        
        W_{SDpath}=((-1)^(SDpath+1)/(freq^(1/pathPowers(SDpath))))*exp(1i*freq*g{1}(realPoints(SDpath))).*dh_dp{SDpath}(x{SDpath}./(freq^(1/pathPowers(SDpath)))).*w{SDpath};
        X_{SDpath}=h{SDpath}(x{SDpath}./freq^(1/pathPowers(SDpath))); 
                
    end
    
    X=[];   W=[];
    for SDpath=1:numPaths
        W=[W; W_{SDpath}];
        X=[X; X_{SDpath};];          
    end

end