function [gStationaryPoints, gSingularities, gSPorders,gPoleOrders] = getStationaryPoints(a,b,rectRad,...
                                                                     analytic, G, RectTol, N , visuals)
    %This one should work better, although the '2' is quite arbitrary.
    if isempty(rectRad)
        rectRad=.75*(b-a);
    end
    initRect=[a-rectRad-rectRad*1i  b+rectRad-rectRad*1i  b+rectRad+rectRad*1i  a-rectRad+rectRad*1i];
   

    if ~analytic
       error('Can only detect stationary points of analytic phase functions'); 
    end
    %assumes g is a cell array containing increasing derivatives of phase
    %function g. Requires at least up to G{3}:=g"(x)
    if length(G)<3
        if ~analytic
            error('Require up to second derivative of phase non-analytic functions, OR an anlytic phase function, to automatically detect stationary points');
        else
            G=finishDerivs( G, 3, N, RectTol );
        end
    end

    %now find all stationary points inside of this rectangle
    [gStationaryPoints, gSingularities, gSPorders, gPoleOrders] = findZerosSingsRect( G{2}, G{3}, initRect, RectTol, N , visuals);
end

