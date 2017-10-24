function err = GenericTestCode( f,g, trueValue, Npts, fSingBit, a, b )

    if nargin<=4
        fSingBit=[];
    end
    if nargin<=3
        Npts=5;
    end
    if nargin<=5
        a=-1;   b=1;
    end
    %add the quadrature paths
    addpath ../../hoipack/quadrules;
    addpath ../../hoipack/OPQ;

    %oscillator(dim,freq,  eval,   deriv, inverse, derivOfInverse, stationaryPoints, orders)
    [ x, w ] = NSD( a, b, Npts, g, fSingBit);

    fprintf('Estimating integral over [%.1f,%.1f], f=%s, oscillator=exp(i*%.0f%s),\n',a,b,func2str(f),g.freq,func2str(g.eval));
    if trueValue==0
        err= abs(f(x).'*w-trueValue);
        fprintf('\tAbsolute error in approximation is %e\n', err);
    else
        err= abs(f(x).'*w-trueValue)/abs(trueValue);
        fprintf('\tRelative error in approximation is %e\n', err);
    end
end

