function I = testLineInt(p1,p2,G,freq)
%quick test function to compute an straight line integral by brute forece,
%from p1 to p2
    L=abs(p2-p1);
    dhdp=(p2-p1)/L;
    [x, w] = quad_gauss(15, 0, L);
    z = p1 + x*dhdp;
    I = w.'*exp(1i*freq*G{1}(z))*dhdp;
end

