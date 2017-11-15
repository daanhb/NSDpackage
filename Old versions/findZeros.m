function allZeros = findZeros( f,df, centre,radius, N )
%computes the Cauchy Derivative of an analytic function
    if nargin<=4
        N=15;
    end

    %theta=(2*pi*(1:N)/(N+1)).'; w=2*pi/N*ones(N,1);
    [theta, w]=quad_gauss(N, 0, 2*pi);
    z=centre+radius*exp(1i*theta);
    
    %I thought the trapezium rule would work well, but it didn't. So using
    %Gaussian quadrature insstead.
    allZeros=w.'*(1/(2*pi*1i)*(1i*radius*exp(1i*theta).*(df(z)./f(z))));
    
    %create 'one' to divide by at the end
    one=w.'*1/(2*pi*1i)*(1i*radius*exp(1i*theta).*(1./((z-centre))));
    
    %now transform output into same shape as input
    allZeros=allZeros/one;
    
end