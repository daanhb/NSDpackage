function c = crudePathTrace(g, startPoint, avoidPoints, r2, r4, visuals)
%crudley traces SD path
    if nargin <=5
        visuals = false;
    elseif visuals
        hold on;
    end
    
    c=startPoint;
    
    finalN = 20;
    
    %initialise the while loop stuff
    r1Max=1; %really not sure how to justify this choice?
    r1=0;
    loopCount=0;
    
    while r4 > abs(c-startPoint)
        if r1<r1Max
             r1=min(min(.5*min(abs(avoidPoints-c)),r1Max));
        end
        N=max(8,ceil(2*pi*r1/r2));
        initCircle = smallCircle(c, r1, N);
        if visuals
            plot(initCircle);
        end
        [~,minPoint] = max(imag(g(initCircle)));
        c = initCircle(minPoint);
        
        loopCount=loopCount+1;
        if loopCount>10000
            error('Cannot trace SD path for some reason');
        end
    end
    finalCircle = smallCircle(c, r2, finalN);
    [~,minPoint] = max(imag(g(finalCircle)));
    c = finalCircle(minPoint);
    if visuals
        plot(finalCircle);
    end
end