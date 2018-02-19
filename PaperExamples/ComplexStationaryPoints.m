clear classes;
syms x f g;

N_Q_set=[4];
omega_set=10.^[2];

%symbolic functions:
%g = symfun(.5*x^2-.5i*x, x); %has stationary points at +i,-i
g = symfun(x^3/3 +x/4,x);
f = symfun(1,x);

%function handles:
%G={@(x) .5.*x.^2-.5i*x, @(x) (x-.5i), @(x) 1};%has stationary points at +i,-i
G={@(z) z.^3./3+z./4, @(z) (z-1i/2).*(z+1i/2), @(z) 2*z, @(z) 2};%has stationary points at +i,-i

Qcount=0;
for N_Q=N_Q_set
    Qcount=Qcount+1;
    omegaCount=0;
    for omega=omega_set
        omegaCount=omegaCount+1;
        I=int(f(x)*exp(1i*omega*g(x)),-1,1);
        [z,w] = NSD45( -1, 1, omega, N_Q, G, 'analytic', true, 'visuals on');
        Q=sum(w);
        R(omegaCount,Qcount)=abs(Q-vpa(I));
    end
end
log10(R)