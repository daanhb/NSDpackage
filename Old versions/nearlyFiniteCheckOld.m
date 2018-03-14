function [yn,p] = nearlyFiniteCheck(G, minDist, c1, c2, r1, r2)
%a necessary (and cheap) condition to determine it a path is nearly finite,
%i.e. almost hits another stationary point
    %local approximation of g close to critical point c2
    %g2loc = @(z) G{r2+1}(c2)/factorial(r2) * (z-c2).^r2;
    g=G{1};
    %theta = acos(real(g(c1)));
    NptsPerArc=100;
    thresh=1E-12;
    
    N =ceil(NptsPerArc*2*pi*minDist);
    theta=linspace(0,2*pi,N+1);
    theta=theta(1:N);
    z = c2 + minDist*exp(1i*theta);
    P = (g(c2+z)-g(c1))/1i;
    %now only keep stuff which is (at least almost) real
    P=P(abs(imag(P))<thresh);
    if isempty(P)
        yn=false;
    else
        yn=true;
        p=min(real(P));
    end
    
    C2=G{r2+1}(c2)/factorial(r2);
    psi2=angle(C2); R2=abs(C2);
    
    theta = (acos(real(g(c1))/(R2*minDist))-psi2+2*(1:r2)*pi)/r2;
    theta = [theta -theta];
    
    z = c2 + minDist*exp(1i*theta);
    P = (g(c2+z)-g(c1))/1i;
    %now only keep stuff which is (at least almost) real
    P=P(abs(imag(P))<thresh);
    if isempty(P)
        yn=false;
    else
        yn=true;
        p=min(real(P));
    end
    
end