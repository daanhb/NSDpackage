function [X, W] = choosePath(a,b,P, G, freq, N, numPathsSD, pathPowers, visuals, X_, W_)
       [pathOrder, pathCost]=findCheapestPath( P, G{1}, freq, N, numPathsSD, pathPowers );%findFullPath( P, G{1}, freq, 1E-10, N, numPathsSD );
        
        if pathCost>.1
            warning('Cost of truncating SD path seems quite high, may not be optimal path');
        end
    
    X=[];   W=[]; xv=[]; yv=[];
    inOut=1;
    for fullIndex=pathOrder
        W=[W; inOut*W_{fullIndex}];
        X=[X; X_{fullIndex};];       
        if inOut==1
            xv=[ xv; real(X_{fullIndex})];
            yv=[yv; imag(X_{fullIndex})];
        else
            xv=[ xv; real(flipud(X_{fullIndex}))];
            yv=[ yv; imag(flipud(X_{fullIndex}))];            
        end
        inOut=inOut*-1;
    end
    if visuals
        for fullIndex=1:numPathsSD
            if ismember(fullIndex,pathOrder)
                plot(X_{fullIndex},'xm'); 
            else
                plot(X_{fullIndex},'oc'); 
            end
        end
        plot([a b],[0 0], 'b', 'LineWidth',2);
    end
    
end

