function [path, cost] = findCheapestPath( P, g,  freq, Npts, m, pathPowers)
%looks at start and 'end' points of each SD path computed, and finds the
%shortest path from a to b

    A=inf(m);

    for j=1:m
       for ell=(j+1):m
           if P(j,1)==P(ell,1) %start points are the same
              A(ell,j)=1E-255; %need this bodge - matlab only connects non-zero entries
           else  %endpoints are possibly 'connected'
              A(ell,j)=pathCost(P(j,2),P(ell,2), g,  freq, Npts);  
           end             
           A(j,ell)=A(ell,j); %symmetric obvz
       end
    end
    
    %now compute shortest (cheapest) path from a to b
    G = graph(A.','upper');
    
    cost=inf;
    for startPath=1:pathPowers(1)
        for endPath=(m-pathPowers(end)+1):m
            [path_, cost_] =shortestpath(G,1,m);
            if cost_<cost
                path=path_;
                cost=cost_;
            end
        end
    end
    
end