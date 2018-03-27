function [G,stationaryPoints,stationaryPointsOrders] = NSDeetsFromPoly(polyCoeffs, stationaryPointMinDist)
%returns all the relevant info needed by NSD given only the coefficients of
%the polynomial
    
    if nargin == 1
        %if two stationary points are this close, make them one
        stationaryPointMinDist=0.001;
    end

    order = length(polyCoeffs)-1;
    
    %use Matlab's roots function to get the data about the stationary
    %points.
    %first factor out the t^2 and compute the roots of what's left:
    
    Dpolycoeffs=polyCoeffs(1:(order)).*fliplr(1:(order));
    stationaryPointsInit = sort(roots(Dpolycoeffs));
    distinctCount=0;
    %determine the orders of the cubic points, by checking for repetitions:
    for n = 1:length(stationaryPointsInit)
        if n>1 & abs(stationaryPointsInit(n)-stationaryPointsInit(n-1)) < stationaryPointMinDist
            stationaryPointsOrders(distinctCount)=stationaryPointsOrders(distinctCount)+1;
        else
            distinctCount = distinctCount+1;
            stationaryPoints(distinctCount) = stationaryPointsInit(n);
            stationaryPointsOrders(distinctCount) = 1;
        end
    end
    
    
    G{1}=@(x) polyval(polyCoeffs,x);
    for n = 2:2*(order+1)
        G{n} = @(x) polyval(Dpolycoeffs,x);
        Dpolycoeffs=Dpolycoeffs(1:(order-n+1)).*fliplr(1:(order-n+1));
    end
    
end

