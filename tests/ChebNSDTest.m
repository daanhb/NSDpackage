clear classes;
%chebfun gives some warnings that I don't understand
warning('off','MATLAB:catenate:DimensionMismatch');
%add the quadrature paths
addpath ../../hoipack/quadrules;
addpath ../../hoipack/OPQ;
addpath ..;
Npts=5; %should be enough
a=-1; b=1;
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
exactValue{1}=-0.0011231043935003419770359866783i;
exactValue{2}=0;

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