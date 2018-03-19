% compute HOI, f(x)=log(x+2), g(x)=x, \omega=1000
addpath ../hoipack/quadrules;
addpath ../hoipack/OPQ;
trueValue =  -0.017347650573975*i;
%true value taken from:
% https://www.wolframalpha.com/input/?i=int_-1%5E1x*exp(100*i*x)dx
f=@(x) x;


%number of points to take in SD code:
Npts=5;

%define the logarithmic singularity at x=-2
fSingBit=[];%singularity(1, 0, 'log', 'point');

%inputs are as follows:
%oscillator(dim,freq,       eval,    deriv,  inverse, stationaryPoints,order)
gPhase=oscillator(1, 100, @(x) x,  @(x) 1, @(x) x,     [],        []);

[ x, w ] = NSDv2( 100,-1,1,Npts, gPhase, fSingBit);

err= abs(f(x).'*w   -  trueValue)/abs(trueValue);
fprintf('Error in approximation is %e', err);