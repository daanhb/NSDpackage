function path = findFullPath( P, g,  freq, thresh, Npts )
%looks at start and 'end' points of each SD path computed, and finds the
%shortest path from a to b

    % init adjacency matrix:
    [m, n]=size(P);
    if n~=2
        error('Matrix of start/end points must by Nx2');
    end

    A=eye(m);

    for j=1:m
       for ell=(j+1):m
           if P(j,2)==P(ell,2) %start points are the same
              A(ell,j)=1;
           elseif negligablePath(P(j,1),P(ell,1), g,  freq, thresh, Npts) %endpoints are essentially the same
              A(ell,j)=1;               
           end
       end
    end
    
    %now compute shortest path from a to b
    [~,path]=dijkstra(A,6,1);
    
end