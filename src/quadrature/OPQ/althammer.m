% ALTHAMMER Althammer polynomials generated by the Stieltjes algorithm.
%
n=21; g=100; s=1; nd=[n n]; a0=0; same=1;
ab=r_jacobi(n); 
zw=gauss(n,ab);
xw=[zw(:,1) zw(:,1) zw(:,2) g*zw(:,2)];
B=stieltjes_sob(n,s,nd,xw,a0,same);
k=6:5:n;
B(2,k)'