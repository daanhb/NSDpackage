function [ stationaryPoints, stationaryPointsOrder, G ] = getStationaryPoints( a,b,g,ZeroThresh )
%if stationary points are not provided, this function estimates them, along
%with derivatives, using Chebfun
    stationaryPoints=[];
    %threshold for how close to zero points need to be to be considered
    %stationary:
    if nargin <=3
        ZeroThresh=1E-15;
    end
    
    %account for multiple types of input:
    if ~iscell(g)
        g_=g;
        clear g;
        g{1}=g_;
        clear g_;
    end

    %if derivative provided, use it, else approximate with Chebfun
    if length(g)>1
        G{2}=chebfun(g{2}, [a b], 'splitting', 'on');
    else
        G{1}=chebfun(g{1}, [a b], 'splitting', 'on');
        G{2}=diff(G{1});
    end
    
    %now look for zeros in a little interval around the true
    %zeros, hence allowing for some machine precision error:
    errRootsAbove=roots(G{2}+ZeroThresh).';
    errRootsBelow=roots(G{2}-ZeroThresh).';
    psudoRoots=[errRootsBelow errRootsAbove ];
    
    if ~isempty(psudoRoots) %not empty
        distinctPsudoRoots(1)=psudoRoots(1);
        k=2;
        for j=2:length(psudoRoots)
            if abs(psudoRoots(j)-psudoRoots(j-1))>2*ZeroThresh
                distinctPsudoRoots(k)=psudoRoots(j);
                k=k+1;
            end
        end
    else
        distinctPsudoRoots=[];
    end
    
    
    for root=1:length(distinctPsudoRoots)
       stationaryPoints(root)=fminbnd(@(x) abs(G{2}), distinctPsudoRoots(root)-2*ZeroThresh, distinctPsudoRoots(root)+2*ZeroThresh);
    end
    
    %now determine the order of these points:
    stationaryPointsOrder = ones(size(stationaryPoints));
    for j=1:length(stationaryPoints)
        n=2;
        while abs(G{n}(stationaryPoints(j))) < ZeroThresh
            %add one to stationary point order
            stationaryPointsOrder(j) = stationaryPointsOrder(j) + 1;
            
            %consider a higher derivative:
            n=n+1;
            %approximate it if not given:
            if length(G)<n
                G{n}=diff(chebfun(G{n-1}, [a b], 'splitting', 'on'));
            end
        end
    end

end