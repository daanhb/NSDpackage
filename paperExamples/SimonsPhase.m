clear classes;
syms x f g;

N_Q_set=[ 4];
omega_set=10.^[2 4 6];

%symbolic functions:
g = symfun((x+1).^-1, x);
f = symfun(1,x);

%function handles:
G={@(x) (x+1).^-1, @(x) -(x+1).^-2, @(x) 2*(x+1).^-3, @(x) -6*(x+1).^-4};

Qcount=0;
for N_Q=N_Q_set
    Qcount=Qcount+1;
    omegaCount=0;
    for freq=omega_set
        omegaCount=omegaCount+1;
        I=int(f(x)*exp(1i*freq*g(x)),0,1);
        [z,w] = NSD45( 0, 1, freq, N_Q, G, 'analytic', true, 'visuals on', 'rectRad',.5/freq);
        Q=sum(w);
        R(omegaCount,Qcount)=abs(Q-vpa(I));
    end
end
log10(R)