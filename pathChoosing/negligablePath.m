function tf = negligablePath( p1,p2, g, freq, thresh, Npts )

    basicPath=linspace(p1,p2,Npts);

    if max(abs(exp(1i*freq*g(basicPath))))<thresh
        tf=true;
    else
        tf=false;
    end
    
end