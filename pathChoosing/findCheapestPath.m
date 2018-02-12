function [path, cost] = findCheapestPath( P, g,  freq, Npts, m )
%looks at start and 'end' points of each SD path computed, and finds the
%shortest path from a to b

    A=inf(m);

    for j=1:m
       for ell=(j+1):m
           if P(j,1)==P(ell,1) %start points are the same
              A(ell,j)=0;
           else  %endpoints are essentially the same
              A(ell,j)=pathCost(P(j,2),P(ell,2), g,  freq, Npts);               
              A(j,ell)=A(ell,j); %symmetric obvz
           end
       end
    end
    
    %now compute shortest (cheapest) path from a to b
    [cost,path]=dijkstra(A,m,1);
    
end