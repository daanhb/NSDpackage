clear;

Npts=10;
Nbruce=25;

freqs=[100, 1000, 10000, 100000];

F={@(z) z.^2, @(z) 1, @(z) log(abs(z))};
singularities={[],[],singularity(1, -2, 'log', 'point')};


G{1}={ @(x)x.^3-.001*x  ,  @(x)3*x.^2-.001  ,  @(x)6*x  ,  @(x)6};
G{2}={@(x)x.^3+.1*x ,   @(x)3*x.^2+.1  ,  @(x)6*x  ,  @(x)6};

TBruce=inf(length(f),length(G));

for j=1:length(f)
    for k=1:length(G)
        for freq=freqs
            %compute exact value by brute force:
            tic;
            [x, wBruce] = Bruce( -1,1,Nbruce,freq );
            exactVal=(F{j}(x).*exp(1i*freq*G{k}{1}(x))).'*w;
            TBruce(j,k)=toc;
            
            %now compute NSD value in 'black box' mode:
            [ z, wNSD ] = NSD45( -1,1,freq(j),Npts,G{k}); 
            NSDval=F{j}(z).'*w;
            
            %also compute it with extra info, and compare timings:
            
        end
    end
end