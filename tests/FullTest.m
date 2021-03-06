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

f{4}=@(x) log(2+x);
fSingBit{4}=singularity(1, -2, 'log', 'point');

f{5}=f{4};
fSingBit{5}=fSingBit{4};

f{6}=@(x) log(1+x);
fSingBit{6}=singularity(1, -1, 'log', 'point');

f{7}=f{6};
fSingBit{7}=fSingBit{6};

%singularity in the interior
f{8}=@(x) log(x);
fSingBit{8}=singularity(1, 0, 'log', 'point');

f{9}=f{1};
fSingBit{9}=[];

f{10}=f{1};
fSingBit{10}=[];

%now define some oscillators:
%oscillator(dim,freq,  eval,   deriv, inverse, derivOfInverse, stationaryPoints, orders)
g{1}=oscillator(1, 1000, @(x) x,  @(x) 1, @(x) x,     @(x) 1,      [],[]);
g{2}=oscillator(1, 100, @(x) x.^2, @(x) 2.*x, {@(x) sqrt(x), @(x) -sqrt(x)}, {@(x) .5*x.^(-1/2), @(x) -.5*x.^(-1/2)},  0, 1);
g{3}=oscillator(1, 1000, @(x) x.^3, @(x) 3*x.^2, {@(x) x.^(1/3), @(x) exp(2i*pi/3).*x.^(1/3), @(x) exp(4i*pi/3).*x.^(1/3)}, {@(x) x.^(-2/3)/3, @(x) exp(2i*pi/3).*x.^(-2/3)/3, @(x) exp(4i*pi/3).*x.^(-2/3)/3},  0, 2);
g{4}=g{1};%just copy this for the singularity test
g{5}=g{2};
g{6}=g{2};
g{7}=g{3};
g{8}=g{1};
%try one with a slightly more interesting phase:
g{9}=oscillator(1, 1000, @(x) x.^2+.5*x, @(x) 2*x+.5, {@(x)-0.25+sqrt(4*x+1/4)/2, @(x)-0.25-sqrt(4*x+1/4)/2}, {@(x) 1./(4*x + 1/4).^(1/2), @(x) -1./(4*x + 1/4).^(1/2)}, -.25, 1);
%now try one with a stationary point at an endpoint:
g{10}=oscillator(1, 1000, @(x) x.^2-1, @(x) 2.*x, {@(x) sqrt(x+1), @(x) -sqrt(x+1)}, {@(x) .5*(x+1).^(-1/2), @(x) -.5*(x+1).^(-1/2)},  [-1 1], [1 1]);

exactValue(1)=-0.0011231043935003419770359866783i;
exactValue(2)=0;
exactValue(3)=0.15521959088497665452635471192667;
exactValue(4)=0.0009080460249256141 - 0.0006167335598745538i;
exactValue(5)=0.08413231737930602 + 0.08208095157921345i;
exactValue(6)=0.00677373541265399 + 0.02598982464981512i;
exactValue(7)=- 0.002467416491118143 + 0.006527549390019811i;
exactValue(8)=- 0.001765639564955236 + 0.002597718689939043i;
exactValue(9)=- 0.006225484373974832 - 0.0131410355084497i;
exactValue(10)=0;

numTests=length(f);
errs=NaN(numTests,1);

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