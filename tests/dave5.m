clear; clc;

x=10; y= -12; z=-3;

%g=@(t) t.^5 + z*t.^3 + y*t.^2 +x.*t;

G{1} = @(t) t.^5 + z*t.^3 + y*t.^2 +x.*t ;
G{2} = @(t) 5*t.^4 + 3*z*t.^2 + 2*y*t +x ;
G{3} = @(t) 20*t.^3 + 6*z*t + 2*y ;

a=-1; b=1; freq=1;

Npts=50;
Nbruce=5000;

[ ~, wNSD ] = NSD45( a, b, freq, Npts, G, 'visuals on','analytic', true);


[x, wBruce] = Bruce( a,b,Nbruce,1 );
exactVal=(exp(1i*freq*G{1}(x))).'*wBruce;


NSDval=sum(wNSD);

relErr=abs(NSDval-exactVal)./abs(exactVal)