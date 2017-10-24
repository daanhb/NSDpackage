clear classes;
%add the quadrature paths
addpath ../../hoipack/quadrules;
addpath ../../hoipack/OPQ;
addpath ..;
Npts=5; %should be enough
a=-1; b=1;
%define some functions for integrands and their corresponding singularities
f{1}=@(x) x;
fSingBit{1}=[];
% 
% f{2}=@(x) log(2+x);
% fSingBit{2}=singularity(1, -2, 'log', 'point');
f{2}=f{1};
fSingBit{2}=[];

f{3}=@(x) ones(size(x));
fSingBit{3}=[];

%now define some oscillators:
%oscillator(dim,freq,  eval,   deriv, inverse, derivOfInverse, stationaryPoints, orders)
g{1}=oscillator(1, 1000, @(x) x,  @(x) 1, @(x) x,     @(x) 1,      [],[]);
g{2}=oscillator(1, 100, @(x) x.^2, @(x) 2.*x, {@(x) sqrt(x), @(x) -sqrt(x)}, {@(x) .5*x.^(-1/2), @(x) -.5*x.^(-1/2)},  0, 1);
g{3}=oscillator(1, 1000, @(x) x.^3, @(x) 3*x.^2, {@(x) x.^(1/3), @(x) exp(2i*pi/3).*x.^(1/3), @(x) exp(4i*pi/3).*x.^(1/3)}, {@(x) x.^(-2/3)/3, @(x) exp(2i*pi/3).*x.^(-2/3)/3, @(x) exp(4i*pi/3).*x.^(-2/3)/3},  0, 2);

exactValue(1)=-0.0011231043935003419770359866783i;
exactValue(2)=0;
exactValue(3)=0.15521959088497665452635471192667;

numTests=length(f);
errs=zeros(numTests,1);

for t=1:numTests
    errs(t) = GenericTestCode( f{t},g{t}, exactValue(t), Npts, fSingBit{t}, a, b );
end

%one day it'd be good to test against a brute force qudarature, and loop
%over everything, to be really thorough, like so:
% for f_index=1:length(f)
%     for g_index=1:length(g)
%         errs(f_index,g_index) = GenericTestCode( f,g, trueValue, Npts, fSingBit, a, b );
%     end
% end