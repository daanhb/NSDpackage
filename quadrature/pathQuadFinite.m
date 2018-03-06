function [ x,w ] = pathQuadFinite( b,S,freq,N )
%special quadrature rule for finite paths, basically just a Gaussian
%quadrature rule, but this might change
m=1; %this used to be an input before I rewrote the code
if ~isempty(S)
    error('Cant handle singularities on finite paths yet');
end
    
        [ x, w_preExp ] =  quad_gauss(N, 0, b^(1/m));
        %no rescaling of this integral:
        w=w_preExp.*exp(-freq*x.^m);

end

