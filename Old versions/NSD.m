function [ X, W ] = NSD( a,b,N,g,singularities )
    %generic Numerical steepest descent routine
    %singularities are all along the real line, and specified by the user
    
    %error in comparing limits of arguemnt of paths at infinity
    errTol=1E-14;

    %currently defunct constant hard-coded for cases where Cauchy integral
    %is needed:
    nCauchy=15;
        
    if nargin<5
        error('Cant automatically determine inverse / derivatives / singularities of phase function yet')    
    end
    
    %check for singularities in (real) interior of integration range:
    interiorSingularities=[];
    for s=singularities
       if  a<s.position && s.position<b
           interiorSingularities=[interiorSingularities s.position];
       end
    end
    
    if ~isequal([],interiorSingularities) 
        %singularities in interior, so let's make like a banana and split
        intervalSplit=[a interiorSingularities b];
        X=[]; W=[];
        for j=1:(length(intervalSplit)-1)
            %call self recursively, without interior singularities:
            [ X_, W_ ] = NSD( intervalSplit(j),intervalSplit(j+1),N,g,singularities );
            X=[X; X_;]; W=[W; W_;];
        end
    else

        %get frequency
        freq=g.freq;
        
        %seperate the interior stationary points
        interiorStationaryPoints=getInteriorStationaryPoints( g.stationaryPoints, a, b );
        
        %determine the number of steepest descent paths
        numPaths=2*length(interiorStationaryPoints)+2;
        
        realPoints=[a repmat(interiorStationaryPoints,1,2) b];

        %initialise local polynomial degree at each path start point
        pathPowers=ones(size(realPoints));
        
        %now check for stationary points, and increase polynomial degree as
        %required:
        for j=1:length(realPoints)
            for k=1:length(g.stationaryPoints)
                if abs(realPoints(j)-g.stationaryPoints(k))<errTol
                    %this is a stationary point, needed higher polynomial
                    %degree
                    pathPowers(j)=g.order(k)+1;
                end
            end
        end

        X=[];   W=[];
            for SDpath=[1 numPaths 2:(numPaths-1)]
                if iscell(g.inverse) %multi-valued inverse, stored in cells
                       %choose the correct branch here, the one that makes sense.
                       %start with the first branch:
                       branch=1;
                       h{SDpath}=@(p) g.inverse{1}(g.eval(realPoints(SDpath))+1i*p.^pathPowers(SDpath));
                        if SDpath == 1 || SDpath == numPaths
                            %check the SD path passes through the endpoints of the
                            %original integral:
                           while abs(h{SDpath}(0)-realPoints(SDpath))>errTol
                               %try a different branch
                               branch=branch+1;
                               h{SDpath}=@(p) g.inverse{branch}(g.eval(realPoints(SDpath))+1i*p.^pathPowers(SDpath));
                           end                
                        else
                            %now compute the limit of the argument of the path as it tends
                            %to inifinty, so the next path can be chosen correctly
                            if SDpath == 2
                                pathEndArg=angle(h{1}(inf));
                            elseif SDpath == (numPaths-1)
                                pathEndArg=angle(h{numPaths}(inf));
                            else
                                error('Not sure how to handle this many stationary points');
                            end
                            %check that the SD paths join up at infinity:
                            while abs(pathEndArg-angle(h{SDpath}(inf)))>errTol
                               branch=branch+1;
                               h{SDpath}=@(p) g.inverse{branch}(g.eval(realPoints(SDpath))+1i*p.^pathPowers(SDpath));
                            end
                            %THIS APPROACH OF USING INFINITY SEEMS A BIT
                            %RISKY. In general, want to just choose the
                            %path which hugs the real axis closest
                        end
                        %function is invertible, so can define SD paths
                        DerivOfInverseOf_g=@(x) g.derivOfInverse{branch}(x);

                else
                    h{SDpath}=@(p) g.inverse(g.eval(realPoints(SDpath))+1i*p.^pathPowers(SDpath));
                    DerivOfInverseOf_g=@(x) 1./(g.inverse(g.deriv(x)));
                end

                %change position of singularities with change of variables
                COVsingularities=singularities;
                for s=1:length(singularities)
                    %THIS MAP ALSO ONLY MAKES SENSE FOR SIMPLE POLYNOMIAL
                    %CASES:
                    COVsingularities(s).position=(singularities(s).position-realPoints(SDpath))*(freq^(1/pathPowers(SDpath)))/1i;
                end
                
                % get weights and nodes for path.
                [x{SDpath}, w{SDpath}]=pathQuad( 0, inf, pathPowers(SDpath), COVsingularities, N );
                
                if pathPowers(SDpath)==1 %non degenerate stationary point or endpoint integral
                    dh_dp{SDpath}=@(p) 1i*DerivOfInverseOf_g(g.eval(realPoints(SDpath))+1i.*p.^pathPowers(SDpath));
                else
                %ULTIMATELY SHOULD INCLUDE SOME TEST HERE IN CASE THE
                %INVERSE DOESN'T CANCEL NICELY, and if it fails, then
                %compute derivatives via the Cauchy integral
                %this may happen for stationary points where Gauss-Hermite
                %is used, as there (may) be a quadrature point exactly at
                %the stationary point, where the derivative of the inverse is undefined

%                     CauchyCentre=mean(x{SDpath}./(freq^(1/pathPowers(SDpath))));
%                     CauchyRadius=max(abs(CauchyCentre-x{SDpath}));
%                     dh_dp{SDpath}=@(a) CauchyDiff( a, h{SDpath}, CauchyCentre, CauchyRadius, nCauchy);
                    dh_dp{SDpath}=@(p) pathPowers(SDpath)*p.^(pathPowers(SDpath)-1)*1i.*DerivOfInverseOf_g(g.eval(realPoints(SDpath))+1i.*p.^pathPowers(SDpath));
                                  %    =i*r*p^{r-1}*(g^-1)'(g(a)+ip^r),
                                  %    (chain rule)
                end
                
                %now append weights and nodes to the rest, noting that we subtract every even path, as these are coming back down from
                %infinity
                W_{SDpath}=((-1)^(SDpath+1)/(freq^(1/pathPowers(SDpath))))*exp(1i*freq*g.eval(realPoints(SDpath))).*dh_dp{SDpath}(x{SDpath}./(freq^(1/pathPowers(SDpath)))).*w{SDpath};
                X_{SDpath}=h{SDpath}(x{SDpath}./freq^(1/pathPowers(SDpath)));  
            end
        X=[];   W=[];
        for SDpath=1:numPaths
            W=[W; W_{SDpath}];
            X=[X; X_{SDpath};];          
        end
    
    end
end