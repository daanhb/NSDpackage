function [ distinctRoots, stationaryPointsOrder, G ] = getStationaryPoints( a,b,g,ZeroThresh )
%if stationary points are not provided, this function estimates them, along
%with derivatives, using Chebfun

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
    
    ChebRoots=roots(G{2},'all');
    
    %this also finds complex roots... which is partly in error when they're
    %small, but in some cases may be significant as they correspond to a
    %stationary point close to the real line
    ChebRootsAreReal=imag(G{2}(ChebRoots))<ZeroThresh;
    
    %now remove imaginary components of roots with sufficiently small
    %imaginary component:
    %1,...,numRoots
    indices=(1:length(ChebRootsAreReal)).';
    %now keep indices corresponding to roots with small imaginary part
    realIndices=indices.*ChebRootsAreReal;
    %now act only on these indices of 
    psudoRoots(realIndices)=real(ChebRoots(realIndices));
    
    %we still need to check if the roots are distinct.
    if ~isempty(psudoRoots) %not empty
        distinctRoots(1)=psudoRoots(1);
        k=2;
        for j=2:length(psudoRoots)
            if abs(G{2}(psudoRoots(j))-G{2}(psudoRoots(j-1)))>ZeroThresh
                distinctRoots(k)=psudoRoots(j);
                k=k+1;
            end
        end
    else
        distinctRoots=[];
    end
    
    %now determine the order of these points:
    stationaryPointsOrder = zeros(size(distinctRoots));
    for j=1:length(distinctRoots)
        n=2;
        while abs(G{n}(distinctRoots(j))) < ZeroThresh
            %add one to stationary point order
            stationaryPointsOrder(j) = stationaryPointsOrder(j) + 1;
            
            %consider a higher derivative:
            n=n+1;
            %approximate it if not given:
            if length(G)<n
                G{n}=diff(chebfun(g, [a b], 'splitting', 'on'), n-1);
            end
        end
    end

end