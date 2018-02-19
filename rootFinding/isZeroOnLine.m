function isZero = isZeroOnLine( a,b,F, lineLength )

    if nargin<=3
        lineLength=abs(b-a);
    end

    thresh=1E-1;

    lineMap=@(s) a+exp(angle(b-a)*1i)*lineLength.*s;
    
    f=@(s) F(lineMap(s));

    [realZeros, id0Real] = AllZeros(@(x) real(f(x)),0,1,1000);
    
    [imagZeros, id0Complex] = AllZeros(@(x) imag(f(x)),0,1,1000);
    
    if ~id0Real && ~id0Complex
        [X,Y] = meshgrid(realZeros,imagZeros);
        zeroDists = min(min(abs(X-Y)));
        if min(min(zeroDists)) < thresh
            isZero = true;
        else
            isZero = false;
        end
    elseif ~id0Real && ~isempty(realZeros)
        isZero = true;
    elseif ~id0Complex && ~isempty(imagZeros)
        isZero = true;
    elseif id0Complex && id0Real
        error('Function seems to be identically zero');
    else
        isZero = false;
    end
    
end

