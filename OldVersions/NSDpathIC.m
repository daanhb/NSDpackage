function LocDh = NSDpathIC( r,SDpath )
        
    
    if mod(SDpath,2)==0
        PathOdd=false;
    else
        PathOdd=true;
    end
        
    %determine locally the branches, one of which will be the start of
    %the path:
    LocPaths=exp(1i*pi/(2*r)*(4*(0:r)+1));

    switch r
        case 1
            %LocDh=1i;
            LocDh=[];
        case 2
            if PathOdd
                LocDh=LocPaths(1);
            else
                LocDh=LocPaths(2);
            end
        case 3
            %this will cause errors, at least at first...
            if PathOdd
                LocDh=LocPaths(1);
            else
                LocDh=LocPaths(2);
            end
        case 4
            error('Not sure about this 4th order stationary points yet...');
    end

end

