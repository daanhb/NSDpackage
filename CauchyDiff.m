function fDiff = CauchyDiff( a,f,centre,radius, diffs, N )
%computes the Cauchy Derivative of an analytic function
    if nargin<=4
        diffs=1;
    end
    if nargin<=5
        N=15;
    end
    
    %evaluation points need to be flat vector for this to work
    a_upright=a(:).';

    if min(abs(a_upright-centre))>=radius
       error('Points X must lie inside disk'); 
    end
    %theta=(2*pi*(1:N)/(N+1)).'; w=2*pi/N*ones(N,1);
    [theta, w]=quad_gauss(N, 0, 2*pi);
    z=centre+radius*exp(1i*theta);
    
    %I thought the trapezium rule would work well, but it didn't. So using
    %Gaussian quadrature insstead.
    fDiff=w.'*(1/(2*pi*1i)*(1i*radius*exp(1i*theta).*(f(z)./((z-a_upright).^(diffs+1)))));
    
    %create 'one' to divide by at the end
    one=w.'*1/(2*pi*1i)*(1i*radius*exp(1i*theta).*(1./((z-a_upright))));
    
    %now transform output into same shape as input
    fDiff=reshape(fDiff,size(a))/one;
    
end