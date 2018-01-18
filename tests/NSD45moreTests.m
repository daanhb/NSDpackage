clear; clc;

% decide what to do about complex stationary points close to endpoints

%finally add check for number of oscillations of function

%add visuals to the uneccessary Cauchy integral bit... and try to fix it

%add the old inverse functionality back in at some point

%singularities close to the SD path will still cause problems

%pimp out the generalised gauss bit

Npts=25;
Nbruce=150;

freqs=[100 10000 100000];

F={@(z) ones( size(z)), @(z) 1./(z-.2i),  @(z) z.^2,  @(z) log(z+1)};
Fsingularities={[],singularity(1,0.2i,'?','point'),[], singularity(1, -1, 'log', 'point'), [], []};
FbruceSingularities={[],[],[],-1};
smallBit=1;%-(.5*1i + 1)^2;


gCount=1;
% %start with analytic phases
%      G{gCount}={@(x) x, @(x) 1, @(x) 0}; gCount=gCount+1;
%      G{gCount}={@(x) .5*x.^2 + x*(-1+smallBit*1i), @(x) x-1-smallBit*1i, @(x) ones(size(x)), @(x) zeros(size(x))}; gCount=gCount+1;
%      G{gCount}={ @(x)x.^3  ,  @(x)3*x.^2 ,  @(x)6*x  ,  @(x)6};gCount=gCount+1;
%       
      G{gCount}={@(x) (1/3)*x.^3 + smallBit*x ,   @(x) x.^2 + smallBit  ,  @(x) 2*x  ,  @(x) 2};  gCount=gCount+1;
%     G{gCount}={@(x) x.^2}; gCount=gCount+1; %use Cauchy derivs to approximate higher derivs here

    %G{gCount}={@(z) 1./(100i-z), @(z)  1./(-z + 100i).^2, @(z)  -2./(z -
    %100i).^3}; gCount=gCount+1; %fails, because you can't remove the
    %residue of exp(1/x), it's infinite, I think...

%     G{gCount}={@(x) x.^2, @(x) 2*x, @(x) 2*ones(size(x))}; Gsingularities{gCount}=[]; gCount=gCount+1;
%     G{gCount} = SqrtTypePhase( .5,1 ); %non-cheating version of what comes later
%     G{gCount}={@(x) x.^.5, @(x) .5*x.^-.5, @(x) -.25*x.^-1.5}; a{gCount}=-0.1; b{gCount}=1; Gsingularities{gCount}=0; gCount=gCount+1;
    
     G{gCount}={@(x) x.^2.*sqrt(x-1), @(x)  (x.* (5*x - 4))./(2 .*sqrt(x - 1)), @(x)  (15*x.^2 - 24*x + 8)./(4* (x - 1).^(3/2))}; gCount=gCount+1;

%the following phase is bad because it has a large negative imaginary
%component at the endpoints, which transforms to a large exponential
%argument on the NSD path:
%    G{gCount}={@(x) x.^3./3-(3i*x.^2)./8-x./8, @(x) (x-1i/2).*(x-1i/4), @(x)  2*x+3i/4, @(x) 2, @(x)1}; gCount=gCount+1;


analyticCount=gCount-1;
%now do the special cases, non-analytic functions, which require more info currently
    %Simon's nasty one:
    G{gCount}={@(x) (x+1).^-1, @(x) -(x+1).^-2, @(x) 2*(x+1).^-3, @(x) -6*(x+1).^-4};
    a{gCount}=0; b{gCount}=1;
    SPs{gCount}=[]; 
    gPoles{gCount}=[];
    gCount=gCount+1;

    %now do an HNA-type integral. This function will provide all of the
    %necessary inputs, although it is cheating a little.
    [G{gCount}, SPs{gCount}, SPOs{gCount} ] = SqrtTypePhase( .5,1 );Gsingularities{gCount}=[];
    a{gCount}=-1; b{gCount}=1;
    gPoles{gCount}=[];
    gCount=gCount+1;
    

%create blank variables to store timings:
TocBruce=inf(length(F),length(G));
TocNSD45=TocBruce;

for j=1:length(F)
    for k=1:length(G)
        %also compute it with extra info, and compare timings:
        if length(G{k})==1
            higherDerivs='(without higher derivatives specified)';
        else
            higherDerivs='';
        end
        if ~isempty(Fsingularities{j})
            singWarning='(f is singular, so brute force will be less acurate)';
        else
            singWarning='';
        end
        
        if isempty(a{k})
            a{k}=-1; b{k}=1;
        end
            
        fprintf('f=%s, g=%s %s%s\n',func2str(F{j}),func2str(G{k}{1}),higherDerivs,singWarning);
        for freq=freqs
            %compute exact value by brute force:
            tic;
            [x, wBruce] = Bruce( a{k},b{k},Nbruce,freq );
            exactVal=(F{j}(x).*exp(1i*freq*G{k}{1}(x))).'*wBruce;
            TocBruce(j,k)=toc;
            
            %now compute NSD value in 'black box' mode:
                tic;
            if k<=analyticCount
                [ z, wNSD ] = NSD45( 0,b{k},freq,Npts,G{k}, 'analytic',true, 'fSingularities',Fsingularities{j}, 'visuals on', true); 
            else
                [ z, wNSD ] = NSD45( a{k},b{k},freq,Npts,G{k}, 'fSingularities',Fsingularities{j},'stationary points',SPs{k}, 'order',SPOs{k}, 'poles',gPoles{k});                
            end
            NSDval=F{j}(z).'*wNSD;
            TocNSD45(j,k)=toc;
            
            fprintf('\tfreq:%.0f\terr:%e\tBruce time:%.3f\tNSD time:%.3f\n', freq, abs(exactVal-NSDval)./abs(exactVal),TocBruce(j,k),TocNSD45(j,k));
            hold off;
        end
    end
end