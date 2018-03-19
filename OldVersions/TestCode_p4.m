addpath ../hoipack/quadrules;
addpath ../hoipack/OPQ;

% compute HOI, f(x)=x, g(x)=x^3, \omega=1000
trueValue = 0.917958 + 0.433517i;%0;% 0.0155972 + 0.143156i;
%true value taken from:
% https://www.wolframalpha.com/input/?i=int_0%5E1+log(x)exp(1000*i*x)dx
f=@(x) x.^3+1;


%number of points to take in SD code:
Npts=15;

fSingBit=[];

%inputs are as follows:
%oscillator(dim,freq,  eval,   deriv, inverse, stationaryPoints)
gPhase=oscillator(1, 10, @(x) x.^4, @(x) 4*x.^3, {@(x) x.^(1/4), @(x) 1i.*x.^(1/4), @(x) -1.*x.^(1/4),@(x) -1i.*x.^(1/4)},   0, 3);

[ x, w ] = NSDv2( 10,-1,1,Npts, gPhase, fSingBit);

err= abs(f(x).'*w -trueValue)/abs(trueValue);
fprintf('Relative rror in approximation is %e', err);

%currently gives correct answer for odd polynomials, and is very wrong for
%even polynomials, of neither polynomials.