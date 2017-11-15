function [ distinctRootsAB, stationaryPointsOrder, G ] = getStationaryPoints( a,b,g,ZeroThresh )
%if stationary points are not provided, this function estimates them, along
%with derivatives, using Chebfun

    %threshold for how close to zero points need to be to be considered
    %stationary:
    if nargin <=3
        ZeroThresh=1E-12;
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
    ChebRootsAreReal=abs(imag((ChebRoots)))<ZeroThresh;
    
    %now remove imaginary components of roots with sufficiently small
    %imaginary component:
    %1,...,numRoots
    indices=(1:length(ChebRootsAreReal)).';
    %now keep indices corresponding to roots with small imaginary part
    realIndices=indices(ChebRootsAreReal); 
    %now act only on these indices of 
    realRoots=real(ChebRoots(realIndices));
    
    %we still need to check if the roots are distinct.
    if ~isempty(realRoots) %not empty
        distinctRoots(1)=realRoots(1);
        k=2;
        for j=2:length(realRoots)
            if abs((realRoots(j))-(realRoots(j-1)))>ZeroThresh%abs(G{2}(realRoots(j))-G{2}(realRoots(j-1)))>ZeroThresh
                distinctRoots(k)=realRoots(j);
                k=k+1;
            end
        end
    else
        distinctRoots=[];
    end
    
    %if root detected is nearly a or b, just switch it to a or b.
    %also, delete it if it's outside of the range [a.b] on real line
    for root=1:length(distinctRoots)
        if abs(a-distinctRoots(root))<ZeroThresh
              distinctRoots(root)=a;
        elseif abs(b-distinctRoots(root))<ZeroThresh
              distinctRoots(root)=b;              
        end
    end
    distinctRootsAB=[];
    for root=1:length(distinctRoots)
        if ~( distinctRoots(root)<a  ||  distinctRoots(root)>b )
              %if inside of interval, keep entry
              distinctRootsAB = [distinctRootsAB distinctRoots(root)] ;              
        end
    end
    
    %now determine the order of these points:
    stationaryPointsOrder = zeros(size(distinctRootsAB));
    for j=1:length(distinctRootsAB)
        n=2;
        while abs(G{n}(distinctRootsAB(j))) < ZeroThresh
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