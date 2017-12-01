clear; clc;

% adapt brute force to handle singularities better

% distinguish better between analytic f and g.

% decide what to do about complex stationary points close to endpoints

%finally add check for number of oscillations of function

%add visuals to the uneccessary Cauchy integral bit

Npts=15;
Nbruce=30;

freqs=[100 1000 10000 100000];

F={@(z) z.^2, @(z) ones(size(z)), @(z) log(z+1)};
Fsingularities={[],[],singularity(1, -1, 'log', 'point'), [], [], [], [], [], [], [], [], [], []};

% %start with analytic phases
%     G{1}={@(x) x, @(x) 1, @(x) 0};
%     G{2}={ @(x)x.^3  ,  @(x)3*x.^2 ,  @(x)6*x  ,  @(x)6};
      smallBit=0.001;
      G{1}={@(x)x.^3+smallBit*x ,   @(x)3*x.^2+smallBit  ,  @(x)6*x  ,  @(x)6};
%     G{4}={@(x) x.^2}; %use Cauchy derivs to approximate higher derivs here
%     G{5}={@(x) x.^2, @(x) 2*x, @(x) 2*ones(size(x))};
%     G{6} = SqrtTypePhase( .5,1 ); %non-cheating version of what comes later
    G{2}={@(x) x.^.5,@(x) .5*x.^-.5, @(x) -.25*x.^-1.5};
    G{3}={@(x) x.^2.*sqrt(x-1), @(x)  (x.* (5*x - 4))./(2 .*sqrt(x - 1)), @(x)  (15*x.^2 - 24*x + 8)./(4* (x - 1).^(3/2))};
    %above function has complex valued for real argument, in particular
    %at it's stationary point at x=0.8. This ballses up the NSD code
    %massively, blowing up the value of the integrand. Currently this gets
    %binned off by the NSD code, as it assumes that it's not on the SD
    %path.

analyticCount=length(G);

%now do the special cases, non-analytic functions, which require more info currently
    %Simon's nasty one:
    G{analyticCount+1}={@(x) (x+1).^-1, @(x) -(x+1).^-2, @(x) 2*(x+1).^-3, @(x) -6*(x+1).^-4};
    a{analyticCount+1}=0; b{analyticCount+1}=1;
    SPs{analyticCount+1}=[]; 
    gPoles{analyticCount+1}=[];

    %now do an HNA-type integral. This function will provide all of the
    %necessary inputs, although it is cheating a little.
    [G{analyticCount+2}, SPs{analyticCount+2}, SPOs{analyticCount+2} ] = SqrtTypePhase( .5,1 );
    a{analyticCount+2}=-1; b{analyticCount+2}=1;
    gPoles{analyticCount+2}=[];
    

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
                [ z, wNSD ] = NSD45( a{k},b{k},freq,Npts,G{k}, 'analytic',true, 'singularities',Fsingularities{j}, 'visuals on', true); 
            else
                [ z, wNSD ] = NSD45( a{k},b{k},freq,Npts,G{k}, 'singularities',Fsingularities{j},'stationary points',SPs{k}, 'order',SPOs{k}, 'poles',gPoles{k});                
            end
            NSDval=F{j}(z).'*wNSD;
            TocNSD45(j,k)=toc;
            
            fprintf('\tfreq:%.0f\terr:%e\tBruce time:%.3f\tNSD time:%.3f\n', freq, abs(exactVal-NSDval)./abs(exactVal),TocBruce(j,k),TocNSD45(j,k));
        end
    end
end