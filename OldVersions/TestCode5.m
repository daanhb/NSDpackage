addpath ../hoipack/quadrules;
addpath ../hoipack/OPQ;

%https://www.wolframalpha.com/input/?i=int_-1%5E1exp(i*x%5E3*1000)dx
% compute HOI, f(x)=1, g(x)=x^3, \omega=1000
trueValue =  0.15521959088497665452635471192667;
%true value taken from:
f=@(x) ones(size(x));

%number of points to take in SD code:
Npts=3;

fSingBit=[];

%inputs are as follows:
%oscillator(dim,freq,  eval,   deriv, inverse, derivOfInverse, stationaryPoints, orders)
gPhase=oscillator(1, 1000, @(x) x.^3, @(x) 3*x.^2, {@(x) x.^(1/3), @(x) exp(2i*pi/3).*x.^(1/3), @(x) exp(4i*pi/3).*x.^(1/3)}, {@(x) x.^(-2/3)/3, @(x) exp(2i*pi/3).*x.^(-2/3)/3, @(x) exp(4i*pi/3).*x.^(-2/3)/3},  0, 2);

[ x, w ] = NSDv3( 1000,-1,1,Npts, gPhase, fSingBit);

%true value is zero
err= abs(f(x).'*w-trueValue)/abs(trueValue);
fprintf('Relative error in approximation is %e', err);