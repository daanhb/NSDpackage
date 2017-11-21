function z0 = isZeroOnLine( a,b,F, lineLength )
%uses Chebfun to see if complex valued f is zero on line in complex plane from a to b

if nargin<=3
    lineLength=abs(b-a);
end

    lineMap=@(s) a+lineLength.*(s+1)./2;
    
    f=@(s) F(lineMap(s));

    ChebF=chebfun(f);
    
    s0=roots(ChebF);
    
    z0=lineMap(s0);
end

