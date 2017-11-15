clear classes;
%chebfun gives some warnings that I don't understand
warning('off','MATLAB:catenate:DimensionMismatch');
%add the quadrature paths
addpath ../../hoipack/quadrules;
addpath ../../hoipack/OPQ;
addpath ..;
Npts=50; %should be enough
a=[-1 -1 0 -1]; b=[1 1 1 1];
small=0.1;
%define some functions for integrands and their corresponding singularities
f=@(x) 1;

%now define some phase functions:
% HNA-type phase
g{1}=@(x) sqrt(x.^2+.5)+x;          
%practical, for me...

% Similar, but seems problematic, unsure why, phase
g{2}=@(x) sqrt(x.^2+.5)+small*x;  
%this has stationary points close to, but not on, the real line. This causes
%Chebfun and my code some issues. Also the function is only very slowly
%oscillating in this region.
%derivative is x/sqrt(x^2 + 1) + 0.1

% Simon's evil phase
g{3}=@(x) 1/(1+x);                  
%causes SD path to curve around and reach its limit at -1
%derivative:  -1/(x + 2)^2

%two-become-one phase
g{4}=@(x) x.^3/3 -small*x;

exactValue{1}= 0.00567695 + 0.000130121i;
exactValue{2}=0.0533585 + 0.0376739i;
exactValue{3}=0.001856901409128378686718514335539748454092417207053726714198799162335708236657735982032609549963554128 - 0.003542796240661998124139831392186183786458610519172310965028181506558171523581193217212479122316590340i;
exactValue{3}=0.002682651243826470967382891147144876496100393707644691966017005260122044730449068044524117089883785086 - 0.004106825681982276378253378158197093219057350340027575172231809930945066211595450567539187834994658126i;
exactValue{4}=0.0269582;

numTests=length(f);
errs=NaN(numTests,1);

for j=4%1:length(g)    
    tic;
    [ x, w ] = ChebNSD45( a(j),b(j),1000,Npts,g{j} );
    relErr45(j)=abs(sum(w)-exactValue{j})/exactValue{j};
    fprintf('\nRel err ChebNSD45 %d: %e\n\t, took %.2f secs', j, relErr45(j),toc);
end