addpath ../hoipack/quadrules;
addpath ../hoipack/OPQ;

% compute HOI, f(x)=x, g(x)=x^2, \omega=1000
trueValue = 0.120225 + 0.116734i;
%true value taken from:
% https://www.wolframalpha.com/input/?i=int_-1%5E1+exp(100*i*x%5E2)dx
f=@(x) ones(size(x));


%number of points to take in SD code:
Npts=5;

%no singularities
fSingBit=[];

%inputs are as follows:
%oscillator(dim,freq,  eval,   deriv, inverse, stationaryPoints)
gPhase=oscillator(1, 100, @(x) x.^2, @(x) 2.*x, {@(x) sqrt(x), @(x) -sqrt(x)}, {@(x) .5*x.^(-1/2), @(x) -.5*x.^(-1/2)},  0, 1);

[ x, w ] = NSDv3( 100,-1,1,Npts, gPhase, fSingBit);

err= abs(f(x).'*w -trueValue)./abs(trueValue);
fprintf('Relative rror in approximation is %e', err);