function costBd = pathCost( p1,p2, g,freq, Npts )
%compute line integral between two points, to decide if they can be
%considered connected in the context of the SD adjacency matrix
    %basicPath=linspace(p1,p2,Npts);
    [z,w] = quad_gauss(ceil(Npts*abs(p2-p1)), p1, p2);

    costBd = abs(w.'*exp(-freq*imag(g(z))));
    
end