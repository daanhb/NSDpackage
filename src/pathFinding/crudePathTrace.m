function c = crudePathTrace(g, startPoint, r1, r2, r4, visuals)
%crudley traces SD path
    if nargin <=5
        visuals = false;
    elseif visuals
        hold on;
    end
    N=max(8,ceil(2*pi*r1/r2));
    initCircle = smallCircle(startPoint, r1, N);
    c=startPoint;
    
    loopCount=0;
    while r4 > abs(c-startPoint)
        if visuals
            plot(initCircle);
        end
        [~,minPoint] = max(imag(g(initCircle)));
        c = initCircle(minPoint);
        initCircle = smallCircle(c, r1, N);
        loopCount=loopCount+1;
        if loopCount>10000
            error('Cannot trace SD path for some reason');
        end
    end
end