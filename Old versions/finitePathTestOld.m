function [P, A] = finitePathTest( SPs, g, thresh )
%function takes in stationary points as 'SPs' and returns matrix of path
%lengths, if finite

    if nargin==2
        thresh=1E-12;
    end
    N=length(SPs);
    P=inf(N,1);
%     g=G{1};
%     Dg=G{2};

    for n=1:N
       for m=1:N
           test=g(SPs(m)-g(SPs(n)))/1i;
          if m==n
              A(m,n)=0;
          elseif abs(real(g(SPs(m)))-real(g(SPs(n))))<thresh && real(test)>0 && abs(imag(test))<thresh
              A(m,n)=real(test);
              P(n)=A(m,n);
          else
              A(m,n)=inf;
          end
       end
    end


end

