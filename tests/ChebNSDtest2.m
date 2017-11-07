clear classes;
%chebfun gives some warnings that I don't understand
warning('off','MATLAB:catenate:DimensionMismatch');
%add the quadrature paths
addpath ../../hoipack/quadrules;
addpath ../../hoipack/OPQ;
addpath ..;
Npts=5; %should be enough
a=-1/2; b=1;
%define some functions for integrands and their corresponding singularities
f{1}=@(x) x;
f{2}=f{1};

%now define some oscillators:
%oscillator(dim,freq,  eval,   deriv, inverse, derivOfInverse, stationaryPoints, orders)
g{1}=@(x) x;
g{2}=@(x) x.^2;

Dg{1}=@(x) 1;
Dg{2}=@(x) 2*x;
DDg{1}=@(x) 0;
DDg{2}=@(x) 2;
exactValue{1}=0.001062211671542962804369482732100805419035433845875674184144423022375204668126136414511545122173897912 - 0.0001200953318397544835610113383386917186534531648898521634009286027931729869868520391739820773603641203i;
exactValue{2}=0.0008987037800369039741257327735620026308282383202138016148530725247030183051850584058902696796105665873 - 0.0001606953855027221739260645887400173433411324630531716205270150067835736162010385266953309413837534089i;

numTests=length(f);
errs=NaN(numTests,1);

for j=1:length(f)
    tic;
    [ x, w ] = ChebNSD( a,b,1000,Npts,g{j},Dg{j},DDg{j} );
    relErr(j)=abs(f{j}(x).'*w-exactValue{j});
    fprintf('\nAbs err ChebNSD %d: %e\n\t, took %.2f secs', j, relErr(j),toc);
    
    tic;
    [ x45, w45 ] = ChebNSD45( a,b,1000,Npts,g{j} );
    relErr45(j)=abs(f{j}(x45).'*w45-exactValue{j});
    fprintf('\nAbs err ChebNSD45 %d: %e\n\t, took %.2f secs', j, relErr45(j),toc);
end