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
        %tolerance for ODE45 solver:
        RelTol=1E-12;
        %tolerance for two statoionary points to become clumped together:
        RectTol=10^-3;
        %distance from an integer to be considered not an integer
        intThresh=0.01;
        
    %flag for plotting stuff as we go
        visuals=false;
        gAnalytic=true;
        
    %default rectangle radius - dependence on frequency not totally clear
    %yet
        rectRad=1;%2/freq;%.5;
        
    %default settle radius
        settleRad=[];
        
        %default inf flags
        ainf=false;
        binf=false;
    
    %% ------------------------------------------------------------------%
    % -------------------- INTERPRET USER INPUT ------------------------ %
    %--------------------------------------------------------------------%

    %check that G is full of function handles - a misplaced space can mess
    %this up
    for j=1:length(G)
        if ~isa(G{j}, 'function_handle') 
            error(sprintf('%dth entry of 5th argument not function handle',j));
        end
    end
    
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
                   
               case 'ganalytic'
                   gAnalytic=varargin{j+1};
               case 'visuals on'
                   visuals=true;
                   hold on;
               case 'ftest' %not necessary, but user can input f, to find singularities etc
                   fTest=true;
                   F=varargin{j+1};
                   if length(F) <2
                       error('need f(x) and its derivative, in cell form');
                   end
               case 'rectrad'
                   rectRad=varargin{j+1};
               case 'settlerad' %radius outside of which function settles down
                   settleRad = varargin{j+1};
                   rectRad = settleRad;
               case 'ainf'
                   ainf=true;
               case 'binf'
                   binf=true;
           end
       end
    end
    
    
    %% ------------------------------------------------------------------%
    % -------------------------- SINGULARITIES ------------------------- %
    %--------------------------------------------------------------------%
         
    %scan for singularities, if requested:
    if gAnalytic
        fSingularities=[]; fSPorders=[];
    elseif fTest
        [~, fSingularities, fSPorders, ~] = findZerosSingsRect( F{1}, F{2}, initRect, RectTol, N , visuals);
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
    % ---------------- STAIONARY and CCRITICAL POINTS ------------------ %
    %--------------------------------------------------------------------%
    
    if isnan(gStationaryPoints) %no stationary points specified by user
        [gStationaryPoints, gSingularities, gSPorders,~] = getStationaryPoints(a,b,rectRad,...
                                                            gAnalytic, G, RectTol, N , visuals,...
                                                            settleRad);
    end
    
    %locate branch points
    branchPoints = getIntegrandBranchPoints(gSPorders, fSPorders, intThresh);
    
    %combine stationary points with [a,b] to make 'critical points'
    [criticalPoints, pathPowers] = makeCriticalPoints(a,b,gStationaryPoints,gSPorders+1, RectTol, gAnalytic, ainf, binf);
    
    
    %% ------------------------------------------------------------------%
    % -------------------------- COMPUTE PATHS ------------------------- %
    %--------------------------------------------------------------------%
    
    % if not enough derivatives provided:
    if length(G)< (round(max(pathPowers))+1)
        error('Need at least %d derivatives of g(x) to compute paths',round(max(pathPowers))+1);
    end
    
    %------now do a load of bodging for the finite path stuff------%
%     
    % **** currently assume that there are no singularities on finite paths
    
    % ***** current code does not account for finite paths stemming from
    % different branches of same stationary point. Will need to replace the
    % vector of distances by a matrix, add an extra loop in
    % getBranchesOfFinitePaths, and also account for the cases where one
    % path goes through three stationary points... aaghsifdhdheh@##$%
    
    % ***** not sure if current finite path code can handle finite paths
    % with multiplicity one at one (or both) ends
    
     [prePathLengthsVec, pathEndIndex, ~, ~] = finitePathTest( criticalPoints, G , pathPowers);
     
     [hFinite, dhdpFinite, Wfinite, FPindices] = getBranchesOfFinitePaths(prePathLengthsVec, pathEndIndex, criticalPoints, pathPowers, G, freq, N, [], RelTol);
    
    %--- end of finite path bodging ---------------------------------%
    
    %now loop over all paths:
    fullIndex=0;    %this will be a unique combination of critPointIndex and branchIndex
    for critPointIndex=1:length(criticalPoints)
        
        %set ICs for SD path differential equation, which will determine SD path:
        ICs = NSDpathICv3( pathPowers(critPointIndex), G, criticalPoints(critPointIndex));
        %initialise this vector map, n'th entry is all the paths indices
        %which start at the n'th critical point
        CritPointToPathInds{critPointIndex}=[];
        for branchIndex=1:pathPowers(critPointIndex)
            fullIndex=fullIndex+1;
        
            if ~FPindices{critPointIndex}(branchIndex)

                COVsingularities=fSingularitiesObj;
                for s=1:length(fSingularitiesObj)
                    COVsingularities(s).position=(fSingularitiesObj(s).position-criticalPoints(critPointIndex))*(freq^(1/pathPowers(critPointIndex)))/1i;
                    %THIS MAP only works for real singularities
                end

                %get weights and nodes
                    [x{critPointIndex, branchIndex}, w{critPointIndex, branchIndex}]=pathQuad( 0, inf, pathPowers(critPointIndex), COVsingularities, N );
                    P0=[0; (x{critPointIndex, branchIndex}./(freq^(1/pathPowers(critPointIndex))))];
                    %**** in above case, weights and nodes are independent of
                    %branchIndex

                %now solve IVP:
                [~,H] = ode45(@(t,y) NSDpathODE(t,y,pathPowers(critPointIndex)-1,G, ICs{branchIndex}, false), P0, ICs{branchIndex}, odeset('RelTol',RelTol) );

                %NEED TO INPUT HIGHER ORDER DE's INTO NSDpathODE.m
                %H(:,1) contains h'(p), H(:,2) contains h(p), will throw away initial
                %points

                if pathPowers(critPointIndex)>1
                    % h'(p) given as output of ODE45, so use it
                    h{critPointIndex,branchIndex}=H(2:end,1);
                    dhdp{critPointIndex,branchIndex}=H(2:end,2);
                else
                    %re-insert h(p) into DE for h'(p)*
                    h{critPointIndex,branchIndex}=H(2:end);
                    dhdp{critPointIndex,branchIndex}=1i./G{2}(H(2:end)); %back into the ODE 
                end

                %check that ODE45 gave a path that actually decays
                %exponentially:
                [divTF,nFinal] = divergenceTest(h{critPointIndex,branchIndex},G{1});
                if divTF
                    h{critPointIndex,branchIndex}=h{critPointIndex,branchIndex}(1:nFinal);
                    dhdp{critPointIndex,branchIndex}(1:nFinal);
                end
                
                %now check if the path we've just created is 'nearly
                %finite'
                [nearlyFinite, nfParamPoint, nfCritPointIndex] = nearlyFiniteCheck(P0, h{critPointIndex,branchIndex}, dhdp{critPointIndex,branchIndex}, criticalPoints(critPointIndex), criticalPoints);
                if nearlyFinite
                    %mash this nearly finite path into a finite path.
                    % update the FPindices vector
                    FPindices{critPointIndex}(branchIndex)=nfCritPointIndex;
                    [hFinite{critPointIndex,branchIndex}, dhdpFinite{critPointIndex,branchIndex}, Wfinite{critPointIndex,branchIndex}]...
                            = nearlyFinitePathFix(nfParamPoint, criticalPoints(nfCritPointIndex), criticalPoints(critPointIndex), pathPowers(critPointIndex),...
                                ICs{branchIndex}, G, freq, N, COVsingularities, RelTol);
                    X_{fullIndex} = hFinite{critPointIndex,branchIndex};
                    W_{fullIndex} = dhdpFinite{critPointIndex,branchIndex}.*Wfinite{critPointIndex,branchIndex};
                    %
                else
                    %absorb h'(p) and other constants into weights, also negate paths
                    %incoming from infinity
                    W_{fullIndex}=...
                        (1/(freq^(1/pathPowers(critPointIndex))))*exp(1i*freq*G{1}(criticalPoints(critPointIndex)))...
                        .*dhdp{critPointIndex,branchIndex}.*w{critPointIndex,branchIndex};

                    X_{fullIndex}=h{critPointIndex,branchIndex}; 
                end
            else %finite paths have already been computed, necessarily to detect if they were finite
                %finite paths shouldn't diverge (in IVP sense) - else they become infinite
                %paths anyway
                X_{fullIndex}=hFinite{critPointIndex,branchIndex};
                W_{fullIndex}=dhdpFinite{critPointIndex,branchIndex}.*Wfinite{critPointIndex,branchIndex}*exp(1i*freq*G{1}(criticalPoints(critPointIndex)));
                %
            end
            
            P(fullIndex,2)=X_{fullIndex}(end);   
            P(fullIndex,1)=criticalPoints(critPointIndex);
            CritPointToPathInds{critPointIndex}=[CritPointToPathInds{critPointIndex} fullIndex];
            PathIndsToCritPoint(fullIndex)=[critPointIndex];
            CritPointEndPoint{critPointIndex}(branchIndex)=X_{fullIndex}(end);
        end
                
    end
    
    numPathsSD=fullIndex; %record total number of paths
    
    %% there is a difference between knowing the path, and walking the path%
    
    %now make a useful array which is indexed by the fullIndex, and returns
    %the fullIndex of every connected finite path
    fullIndex=0;
    for critPointIndex=1:length(criticalPoints)
        for branchIndex=1:pathPowers(critPointIndex)
            fullIndex=fullIndex+1;
            if FPindices{critPointIndex}(branchIndex)>0
                FPfullIndices{fullIndex}=CritPointToPathInds{FPindices{critPointIndex}(branchIndex)};
            else
                FPfullIndices{fullIndex}=[];
            end
        end
    end
    
    [X, W] = choosePathV2(a,b,criticalPoints, CritPointEndPoint, FPindices, P, G, freq, N, numPathsSD, visuals, X_, W_, hFinite, settleRad, ainf, binf);
    
    if ~gAnalytic || ~isempty(fSingularities)
        %add tiny circles around singular points if they lie inside the region
        %of deformation:
        [XsmallDisk, WsmallDisk] = handleSingularities(fSingularities, gSingularities, branchPoints, freq, G, N, visuals);
        X=[X; XsmallDisk];  W=[W; WsmallDisk];
    end
    
    %stop adding stuff to the lovely diagram
    if visuals
        hold off;
    end
end