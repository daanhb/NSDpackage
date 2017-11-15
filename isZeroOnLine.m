function z0 = isZeroOnLine( a,b,F )
%uses Chebfun to see if complex valued f is zero on line in complex plane from a to b

    lineMap=@(s) a+(b-a).*(s+1)./2;
    
    f=@(s) F(lineMap(s));

    ChebF=chebfun(f);
    
    s0=roots(ChebF);
    
    z0=lineMap(s0);
end

