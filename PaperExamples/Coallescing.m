clear classes;
syms x f g c;

freq = 100;
N_Q = 5;

%symbolic function:
g = symfun(1/3*x^3-c*x, [x c]);
f = symfun(1,x);

counter=0;
for smallBit = [ 5E-4 .01 5E-10 0]
    counter=counter+1;
    G{counter}={@(x) (1/3)*x.^3 - smallBit*x ,   @(x) x.^2 - smallBit  ,  @(x) 2*x  ,  @(x) 2, @(x) 0};
    [z,w] = NSD45( -1, 1, freq, N_Q, G{counter}, 'ganalytic', true, 'visuals on', 'RectRad', 1);
    I=int(f(x)*exp(1i*freq*g(x, smallBit)),-1,1);
    Q=sum(w);
    R(counter)=abs(Q-vpa(I))/abs(vpa(I));
end

log10(R)