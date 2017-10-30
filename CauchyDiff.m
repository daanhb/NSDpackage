function fDiff = CauchyDiff( a,f,centre,radius, N )
%computes the Cauchy Derivative of an analytic function
    if nargin==4
        N=15;
    end
    
    %evaluation points need to be flat vector for this to work
    a_upright=a(:).';

    if min(abs(a_upright-centre))>=radius
       error('Points X must lie inside disk'); 
    end
    theta=2*pi*(1:N)/(N+1);
    [theta, w]=quad_gauss(N, 0, 2*pi);
    z=centre+radius*exp(1i*theta);
    %w=2*pi/N;
    
    %I thought the trapezium rule would work well, but it didn't. So using
    %Gaussian quadrature insstead.
    fDiff=w.'*1/(2*pi*1i)*(1i*radius*exp(1i*theta).*(f(z)./((z-a_upright).^2)));
    
    %now transform output into same shape as input
    fDiff=reshape(fDiff,size(a));
end