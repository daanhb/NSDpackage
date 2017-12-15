clear classes; clc;

Npts=15;
Nbruce=150;

freqs=[100 1000 10000 100000];

%have decided there will be general type of HNA phase:
%  g(z) = (z^2+az+b)^.5 + (z^2+cz+d)^.5
%where:
%  a^2-4b >=0,   c^2-4d>=0.

% will feed code vectors of form J=[a b c d]

%can do everything with anonymous functions, like so:
CoDist=@(t,a,b) sqrt(t.^2+a*t+b);
CoDistD1=@(t,a,b) (a + 2*t)./(2*sqrt(t.*(a + t) + b));
CoDistD2=@(t,a,b) (4*b - a.^2)./(4*(t.*(a + t) + b).^(3/2));
%define phase function:
g=@(z,a,b,c,d) dist(z,a,b)+dist(z,c,d);
%and its derivatives:
gD1=@(z,a,b,c,d) CoDistD1(z,a,b) + CoDistD1(z,c,d);
gD2=@(z,a,b,c,d) CoDistD2(z,a,b) + CoDistD2(z,c,d);
%and non-oscillatory bit of integrand:
f=@(z,k,a,b,c,d) -1/8*besselh(0,1,k*dist(z,a,b)).*besselh(0,1,k*dist(z,c,d))./exp(1i*k*g(z,a,b,c,d));

JALL{1}=[2 1/2 3 2];
Jcount=0;

for J=JALL
    Jcount=Jcount+1;
    G={@(z) g(z,J{1}(1),J{1}(2),J{1}(3),J{1}(4)), @(z) gD1(z,J{1}(1),J{1}(2),J{1}(3),J{1}(4)), @(z) gD2(z,J{1}(1),J{1}(2),J{1}(3),J{1}(4))};
    
    for k=freqs
        F=@(z) f(z,k,J{1}(1),J{1}(2),J{1}(3),J{1}(4)) ;
        %do the bruce thing:
        tic;
        [x, wBruce] = Bruce( -1,1,Nbruce,k );
        exactVal=(F(x).*exp(1i*k*G{1}(x))).'*wBruce;
        TocBruce(Jcount)=toc;
            
        tic;
        [ z, wNSD ] = NSD45( -1,1,k,Npts,G, 'analytic',true,'visuals on', true); 
        NSDval=F(z).'*wNSD;
        TocNSD45(Jcount)=toc;
        clear F;
    end    
    clear G;
end
