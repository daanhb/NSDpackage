function [X_,W_] = abortPathSearch(a,b,ainf,binf,freq,Npts)
%if no path can be found, this file is called
    if ainf
        a=a*R;
    end
    if binf
        b=b*R;
    end
    [X_{1},W_{1}] = oscQuadExpensive(a,b,freq, Npts);
end

