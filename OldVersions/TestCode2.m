% compute HOI, f(x)=log(x+2), g(x)=x, \omega=1000
trueValue =  -0.00157023 - 0.00748414*i;
%true value taken from:
% https://www.wolframalpha.com/input/?i=int_0%5E1+log(x)exp(1000*i*x)dx
f=@(x) log(x);


%number of points to take in SD code:
Npts=5;

%define the logarithmic singularity at x=-2
fSingBit=singularity(1, 0, 'log', 'point');

%inputs are as follows:
%oscillator(dim,freq,  eval,   deriv, inverse, stationaryPoints, orders)
gPhase=oscillator(1, 1000, @(x) x,  @(x) 1, @(x) x,     [],         []);

[ x, w ] = NSDv2( 1000,0,1,Npts, gPhase, fSingBit);

err= abs(f(x).'*w   -  trueValue)/abs(trueValue);
fprintf('Error in approximation is %e', err);