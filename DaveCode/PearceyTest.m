%test a general polynomial case, for Dave

%changes to make:
    % get Freud quads from Daan
    % actually test some of these
     
%% play-with-able parameters:
Npts=15; %number of quadrature points per SD path
x=6; y=2; % x in {-8,-6,-4,-2,0,2,4,6,8} & y in {0,2,4,6,8}
polyCoeffs = [1 0 x y 0]; %a_1*x^N + ... a_N*x + a_{N+1}, if a:=polyCoeffs
a=-1; b=1; %directions of valleys

%% ----------------------------------------------------------------------------------------- %%
addpath ../..;
StandardPaths;
%make the polynomial
order = length(polyCoeffs)-1;
fprintf('\nPhase is an order %d polynomial',order);
%first derivative
D1polyCoeffs = polyCoeffs(1:order).*fliplr(1:order);
%second derivative
D2polyCoeffs = D1polyCoeffs(1:(order-1)).*fliplr(1:(order-1));

G = {@(x) polyval(polyCoeffs,x),@(x) polyval(D1polyCoeffs,x),@(x) polyval(D2polyCoeffs,x)};

%% ----------------------------------------------------------------------------------------- %%
%some quick error tests:
if aValleyIndex==bValleyIndex
    error('Start and end valleys are the same, choose something more interesting');
end
if max(aValleyIndex,bValleyIndex)>order
    error(sprintf('There are only %d valleys', order));
end
%% ----------------------------------------------------------------------------------------- %%

%this radius could be considered artificial infinity, outside of which all SD
%paths won't deviate much from a straight line
R = monomialSettleRadius(polyCoeffs)+rand;
fprintf('\nSettled radius = %.1f\n',R);

%close; %close figure, if there is one
freq=1; %stick to low frequency regime for now
[ X, W ] = NSD45( a, b, freq, Npts, G, 'analytic', true, 'visuals on', ...
                    'settleRad', R, 'ainf', 'binf');
         %second row are all options for infinite endpoints
fprintf('\nPink quad points are on chosen path, blue line is discarded path(s)');
         %now time it (without making visuals)
tic;
NSD45( a, b, freq, Npts, G, 'analytic', true, ...
                    'settleRad', R, 'ainf', 'binf');
T=toc;
fprintf('\nTook %f seconds',T);

I_GHH=sum(W); %Gibbs-Hewett-Huybrechs estimate
I_CHK=KirkPearceyData(x,y); %Conor-Hobbs-Kirk estimate
fprintf('\nRelative error (against 6dp of Pearcey data) %e',abs(I_GH-I_K)/abs(I_K));