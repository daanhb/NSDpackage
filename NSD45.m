function [ X, W ] = NSD45( a,b,freq,N,G,varargin)
%Where important information about the phase g(x) is unknown, this code
%will approximate the information using and ODE45.
%We do ask however, that the user provides derivatives of g

%if there is a stationary point just outside of [a,b] on the real line what do we do about
%it? It will get picked up by the root finder. If we can plot the level
%curvers from h_a(Pmax) to h_b(Pmax), can check for stationary points
%inside of this region... maybz

%if the points on the SD path end up outside of the initial rectangle, then
%there might be stationary points nearby which should be considered; in
%fact it may be sensible to start over with a larger rectangle, to be
%totally robust.

    %% ------------------------------------------------------------------%
    % -------------------------- KEY PARAMETERS ------------------------ %
    %--------------------------------------------------------------------%
    
    %loads of different values for 'almost zero'
        %tolerance for finding stationary points. Higher values make it more
        %likely to correctly determine their order, but less acurate locations
        errTol=1E-14;
        %tolerance for ODE45 solver:
        RelTol=1E-13;
        %tolerance for two statoionary points to become clumped together:
        RectTol=10^-3;
        %abs value estimate at edge of rectangle
        RectPerimThresh=1E-16;
        %distance from an integer to be considered not an integer
        intThresh=0.01;
        rectSretch=10;
    
    %and one value for 'quite big'
        %to identify exponential growth, instead of decay:
        wrongPathThresh=2;
        
    %flag for plotting stuff as we go
        visuals=false;
        analytic=false;
    
    %% ------------------------------------------------------------------%
    % -------------------- INTERPRET USER INPUT ------------------------ %
    %--------------------------------------------------------------------%

    if length(varargin)==1 && ~ischar(varargin{1})
        %glitchy Matlab varagin thing, only an issue when this function 
        %call's it's self recursively. Fix it here:
        varargin=varargin{1};
        %MESSES UP WHEN YOU SEND MATLAB ONLY ONE OPTION, DUH
    end
    %use NaN to store things which are not yet defined, as opposed to
    %empty (which is denoted [], and may be specified by the user):
    gStationaryPoints=NaN;
    fSingularities=[]; fSingularitiesObj=[]; %singularities in f, not g
    Mf=1; fTest=false;
    interiorSingularities=[]; fSingularities=[];
    %check through optional inputs
    for j=1:length(varargin)
        if ischar(varargin{j})
           lowerCaseArg=lower(varargin{j});
           switch  lowerCaseArg
               case 'stationary points'
                   gStationaryPoints=varargin{j+1};
                   if isempty(gStationaryPoints)
                       gSPorders=[];
                   end
               case 'fsingularities'
                   fSingularitiesObj=varargin{j+1};
                    %check for singularities in (real) interior of integration range:
                    for s=fSingularitiesObj
                        if imag(s.position)==0
                           if  a<s.position && s.position<b
                               interiorSingularities=[interiorSingularities s.position];
                           end
                        else %keep hold of these for later, incase deformation crosses them
                            fSingularities=[fSingularities s.position];
                        end
                    end
               case 'gSingularities'
                   gSingularities=varargin{j+1};
               case 'order'
                   gSPorders=varargin{j+1};               
               case 'mf' %upper bound of non-oscillatory f in complex plane
                   Mf=varargin{j+1};
               case 'ginv' %user has provided inverse of g(x)

               case 'gderinv' %user has provided inverse of g'(x)
                   
               case 'analytic'
                   analytic=varargin{j+1};
               case 'visuals on'
                   visuals=true;
                   hold on;
               case 'ftest' %not necessary, but user can input f, to find singularities etc
                   fTest=true;
                   F=varargin{j+1};
                   if length(F) <2
                       error('need f(x) and its derivative, in cell form');
                   end
           end
       end
    end
    
    
    %% ------------------------------------------------------------------%
    % -------------------------- SINGULARITIES ------------------------- %
    %--------------------------------------------------------------------%
    
    
    %This one should work better, although the '2' is quite arbitrary.
    rectRad=2*(b-a);
    initRect=[a-rectRad-rectRad*1i  b+rectRad-rectRad*1i  b+rectRad+rectRad*1i  a-rectRad+rectRad*1i];
        
    %scan for singularities, if requested:
    if analytic
        fStationaryPoints=[]; fSingularities=[]; fSPorders=[];
    elseif fTest
        [fStationaryPoints, fSingularities, fSPorders, ~] = findZerosSingsRect( F{1}, F{2}, initRect, RectTol, N , visuals);
    end
    
    
    %if there are any interior singularities, run iteratively with them at
    %the endpoint in each case
    if ~isequal([],interiorSingularities) 
        %singularities in interior, so let's make like a banana and split
        intervalSplit=[a interiorSingularities b];
        X=[]; W=[];
        for j=1:(length(intervalSplit)-1)
            %call self recursively, without interior singularities:
            %vararginCopy=varargin;
            [ X_, W_ ] = NSD45( intervalSplit(j),intervalSplit(j+1),freq,N,G,varargin);
            X=[X; X_;]; W=[W; W_;];
        end
        %and end routine 'early'
        return;
    end 
    
    
    %% ------------------------------------------------------------------- 
    % ------------------------ STAIONARY POINTS ------------------------ %
    %--------------------------------------------------------------------%
    
    if isnan(gStationaryPoints) %no stationary points specified by user
        if ~analytic
           error('Can only detect stationary points of analytic phase functions'); 
        end
        %assumes g is a cell array containing increasing derivatives of phase
        %function g. Requires at least up to G{3}:=g"(x)
        if length(G)<3
            if ~analytic
                error('Require up to second derivative of phase non-analytic functions, OR an anlytic phase function, to automatically detect stationary points');
            else
                G=finishDerivs( G, 3, N, RectTol );
            end
            
            %could potentially make a hack where only g' is required. This can
            %then be used to find singularities of g, and then a Cauchy
            %integral can be used to determine derivatives, avoiding
            %singularities / subtracting residues
        end
        
        %now find all stationary points inside of this rectangle
        [gStationaryPoints, gSingularities, gSPorders, ~] = findZerosSingsRect( G{2}, G{3}, initRect, RectTol, N , visuals);
    end
    
    %now determine non-singular branch points
    branchPoints=[];
    allStationaryPoints=[gStationaryPoints fStationaryPoints];
    allSPorders=[gSPorders fSPorders];
    for j=1:length(allSPorders)
        if ~isNearlyInt( allSPorders(j), intThresh )
            branchPoints=[branchPoints allSPorders(j)];
        end
    end
%     for j=1:length(gSingularities)
%         if ~isNearlyInt( gSingOrder(j), intThresh )
%             branchPoints=[branchPoints gSingularities(j)];
%         end
%     end
    
    %THROW AWAY STATIONARY POINTS WHICH ARE ALONG A PATH OF ASCENT
    %now this is fairly uneccessary, since the path is chosen wisely later,
    %however this will speed things up a bit
%     copyStationaryPoints=gStationaryPoints;
%     gStationaryPoints=[];
%     copyOrder=gSPorders;
%     gSPorders=[];
%     spCount=0;
%     for j=1:length(copyStationaryPoints)
%        if exp(1i*G{1}(copyStationaryPoints(j))*freq )<wrongPathThresh %SD path decreases value of integrand
%            spCount=spCount+1;
%            gStationaryPoints(spCount)=copyStationaryPoints(j);
%            gSPorders(spCount)=copyOrder(j);
%        end
%         %otherwise stationary point will not be on SD path, so can be
%         %ignored
%     end
    
    %sort stationary points in order of real components:
    %stationaryPoints=sort(stationaryPoints,'ComparisonMethod','real');
    
    numPathsSD=2*length(gStationaryPoints);
    startPath=[gStationaryPoints gStationaryPoints];
    startPath=sort(startPath,'ComparisonMethod','real');
    pathPowers=round([gSPorders gSPorders]) + 1;
    
    %now determine if a or b are at a stationary point
    if isempty(gStationaryPoints)
        startPath=[a b];
        pathPowers=[1 1];
        numPathsSD=2;
    else
        if min(abs(gStationaryPoints-a))>RectTol
            startPath=[a startPath];
            pathPowers=[1 pathPowers];
            numPathsSD=numPathsSD+1;
        else
            %delete repeated path, and swap remaining one for endpoint a
            startPath=[a startPath(3:end)];
            pathPowers=pathPowers(2:end);
            numPathsSD=numPathsSD-1;
        end
        if min(abs(gStationaryPoints-b))>RectTol
            startPath=[startPath b];
            pathPowers=[pathPowers 1];
            numPathsSD=numPathsSD+1;
        else
            %delete repeated path, and swap remaining one for endpoint a
            startPath=[startPath(1:(end-2)) b];
            pathPowers=pathPowers(1:(end-1));
            numPathsSD=numPathsSD-1;
        end
    end
    
    
    %CHANGE
    %Can we use Cauchy differentiation formula to estimate order of stationary
    %points a little better?
    %will need to locate singularities, and branch cuts... and then move
    %branch cuts to compute g^(n)(s)?=0. Pretty complicated.
    
    %seperate the interior stationary points
%     interiorStationaryPoints=getInteriorStationaryPoints( stationaryPoints, a, b );
%     %THIS FUNCTION SHOULD BE CHANGED TO JUST EXCLUDE STATIONARY POINTS AT
%     %{a,b}
%     
% %     %determine the number of steepest descent paths
% %     numPaths=2*length(interiorStationaryPoints)+2;
% % 
% %     realPoints=[a repmat(interiorStationaryPoints,1,2) b];
% 
%     %initialise local polynomial degree at each path start point
%     pathPowers=ones(size(realPoints));
%     
%     % THE REST of this section needs re-writing:
%     
%     if isempty(branchPoints)
%         %now check for stationary points, and increase polynomial degree as
%         %required:
%         for j=1:length(realPoints)
%             for k=1:length(stationaryPoints)
%                 if abs(realPoints(j)-stationaryPoints(k))<RectTol
%                     %this is a stationary point, needed higher polynomial
%                     %degree
%                     pathPowers(j)=round(order(k))+1;
%                 end
%             end
%         end
%     else
%         error('havent accounted for branch cuts/points yet');
%         %need to have a think about this. Can move branch cuts, but this is
%         %inherently messy.
%     end
%     
%     %each path will start/end at one of the following points:
%     realPoints=a;
%     for j=1:length(stationaryPoints)
%         realPoints=[ realPoints stationaryPoints(j) stationaryPoints(j)];
%     end
%     realPoints=[realPoints b];
 
    
    %% ------------------------------------------------------------------%
    % -------------------------- COMPUTE PATHS ------------------------- %
    %--------------------------------------------------------------------%
    
    % if not enough derivatives provided, approximate them by Cauchy
    %differentiation formula
    if analytic && length(G)<max(pathPowers)+1
        G = finishDerivs( G, max(pathPowers)+1, N, RectTol );
    end
    
    [prePathLengthsVec, pathEndPoint, finitePathPowers, ~] = finitePathTest( startPath, G , pathPowers);

    %initlaise matrix for path finding
    P=NaN(numPathsSD,2);
    
    %------now do a load of bodging for the finite path stuff------%
    % CREATE #SA here, equal to number of finite paths
    numPathsSA=sum(prePathLengthsVec<inf);
    numPaths=numPathsSD+numPathsSA;
    
    %now bolt the finite path powers on the end
    endPowers=finitePathPowers(:,2).';
    endPowers=endPowers(~isnan(endPowers));
    pathPowers=[pathPowers endPowers];
    
    %now split the finite path lengths in half
    for j=1:length(prePathLengthsVec)
        if prePathLengthsVec(j)<inf
            prePathLengthsVec(end+1)=prePathLengthsVec(j)/2;
            prePathLengthsVec(j)=prePathLengthsVec(j)/2;
        end
    end
    %PERHAPS DOESN'T NEED TO BE SPLIT IN HALF, could be something else
    
    %and the start points
    startPath=[startPath pathEndPoint(~isnan(pathEndPoint)).'];
    
    %--- end of finite path bodging ---------------------------------%
    
    %now loop over all paths:
    for SDpath=1:numPaths
        
          %do steepest descent paths first, then bolt the ascent paths (if
          %there are any) on at the end
          if SDpath>numPathsSD
              ascFlag=true;
          else
              ascFlag=false;
          end
          
%         try %ODE45 may fail to find certain paths, although these typically seem to be paths which we do not take
%             %change position of singularities with change of variables
            COVsingularities=fSingularitiesObj;
            for s=1:length(fSingularitiesObj)
                COVsingularities(s).position=(fSingularitiesObj(s).position-startPath(SDpath))*(freq^(1/pathPowers(SDpath)))/1i;
                %THIS MAP only works for real singularities
            end

            %get weights and nodes
            if isinf(prePathLengthsVec(SDpath))
                %[x{SDpath}, w{SDpath}]=pathQuad( 0, (prePathLengthsVec(SDpath)/2)^(1/pathPowers(SDpath)), pathPowers(SDpath), COVsingularities, N );
                [x{SDpath}, w{SDpath}]=pathQuad( 0, inf, pathPowers(SDpath), COVsingularities, N );
                P0=[0; (x{SDpath}./(freq^(1/pathPowers(SDpath))))];
            else %finite path
                [x{SDpath}, w{SDpath}] = pathQuadFinite( prePathLengthsVec(SDpath), pathPowers(SDpath), COVsingularities, freq, N );
                P0=[0; x{SDpath}];
            end
            %THE FACTOR OF A HALF FOR THE FINITE PATHS IS QUITE ARBITRARY

            %set ICs for SD path differential equation, which will determine SD
            %path:
            switch pathPowers(SDpath)
                case  1
                    ICs=startPath(SDpath) ;
                case 2
                    ICs=[ startPath(SDpath); NSDpathICv2( pathPowers(SDpath), (-1)^(SDpath+1), G,  startPath(SDpath), ascFlag ) ];
                otherwise
                    ICs=[ startPath(SDpath); NSDpathICv2( pathPowers(SDpath), (-1)^(SDpath+1), G,  startPath(SDpath), ascFlag ); zeros(pathPowers(SDpath)-2,1) ];
            end

            %now solve IVP:
            [~,H] = ode45(@(t,y) NSDpathODE(t,y,pathPowers(SDpath)-1,G, ICs, ascFlag), P0, ICs, odeset('RelTol',RelTol) );
    %         warnODE45 = warning('query','last');
    %         if 

            %NEED TO INPUT HIGHER ORDER DE's INTO NSDpathODE.m
            %H(:,1) contains h'(p), H(:,2) contains h(p), will throw away initial
            %points

            if pathPowers(SDpath)>1
                % h'(p) given as output of ODE45, so use it
                h{SDpath}=H(2:end,1);
                dhdp{SDpath}=H(2:end,2);
            else
                %re-insert h(p) into DE for h'(p)*
                h{SDpath}=H(2:end);
                dhdp{SDpath}=1i./G{2}(H(2:end)); %this is the DE for 
            end

            %absorb h'(p) and other constants into weights, also negate paths
            %incoming from infinity
            if isinf(prePathLengthsVec(SDpath))
                W_{SDpath}=(1/(freq^(1/pathPowers(SDpath))))*exp(1i*freq*G{1}(startPath(SDpath))).*dhdp{SDpath}.*w{SDpath};
            else
                W_{SDpath}=dhdp{SDpath}.*w{SDpath};
            end

            X_{SDpath}=h{SDpath}; 
            
            P(SDpath,1)=startPath(SDpath);
            P(SDpath,2)=X_{SDpath}(end);   
%         catch  %shall not be used, ODE45 failed
%             W_{SDpath}=[]; 
%             X_{SDpath}=NaN;
%             
%             P(SDpath,1)=NaN;
%             P(SDpath,2)=NaN;   
%         end
                
    end
    
    %% there is a difference between knowing the path, and walking the path%
    try
        pathOrder=findFullPath( P, G{1}, freq, 1E-10, N, numPathsSD );
    catch
        error('Could not find suitable SD path from a to b :-(');
    end
    
    X=[];   W=[]; xv=[]; yv=[];
    inOut=1;
    for SDpath=pathOrder
        W=[W; inOut*W_{SDpath}];
        X=[X; X_{SDpath};];       
        if inOut==1
            xv=[ xv; real(X_{SDpath})];
            yv=[yv; imag(X_{SDpath})];
        else
            xv=[ xv; real(flipud(X_{SDpath}))];
            yv=[ yv; imag(flipud(X_{SDpath}))];            
        end
        inOut=inOut*-1;
    end
    if visuals
        for SDpath=pathOrder
           plot(X_{SDpath},'x'); 
        end
        %plot(X,'kx','LineWidth',1.5);
        %plot(xv,yv, 'r', 'LineWidth',2);
        plot([a b],[0 0], 'b', 'LineWidth',2);
        plot(fSingularities,'ro','LineWidth',2);
        hold off;
    end
    
    %for s=imagSingularities
    if ~isempty(fSingularities)
            [fin,fon] = inpolygon(real(fSingularities),imag(fSingularities),xv,yv);
        %end
        if max(fon)
            error('f has singularity on SD path, risky business');
        end
        if max(fin)
            %got some residues to compute. Choose the radius to be small enough
            %that the integral will be non-oscillatory:
            singRad=1/freq;
            clear X_ W_;
            for s=fSingularities(fin)
                [X_, W_] = residueQuad( s,singRad, N );
                X=[X; X_]; W=[W; W_.*exp(1i*freq*G{1}(X_))];
            end
        end
    end
    
    %for s=imagSingularities
    badPoints=[gSingularities branchPoints];
    if ~isempty(badPoints)
        [gin,gon] = inpolygon(real(badPoints),imag(badPoints),xv,yv);
        if max(gin) || max(gon)
            error('Singularity in g(z) or branch point in region of deformation');
        end
    end
    
end