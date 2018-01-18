function [P, endPoint, finitePathPowers, A] = finitePathTest( SPs, G, pathPowers, thresh )
%function takes in stationary points as 'SPs' and returns matrix of path
%lengths, if finite

    if nargin==3
        thresh=1E-12;
    end
    N=length(SPs);
    P=inf(N,1);
    finitePathStart=false(N,1);
    finitePathEnd=false(N,1);
    endPoint=NaN(N,1);
    finitePathPowers=NaN(N,2);
    g=G{1};
%     Dg=G{2};
    %approximate h_\sigma'(0+small)
    Dh=@(n) SPs(n)  +  thresh*NSDpathICv2( pathPowers(n), (-1)^(n+1), G,  SPs(n) );

    for n=1:N
       for m=1:N
              test=(g(SPs(m))-g(SPs(n)))/1i;
          if SPs(m)==SPs(n)
              A(m,n)=0;
          elseif abs(real(g(SPs(m)))-real(g(SPs(n))))<thresh && real(test)>0 && abs(imag(test))<thresh && imag(G{2}(Dh(n)))>=0
              A(m,n)=real(test);
              %P(n)=A(m,n);
          else
              A(m,n)=inf;
          end
       end
    end
    
    %in case it's possible that one finite path crosses two stationary
    %points, the following will truncate the path length at the first:
    for n=1:N
        for m=1:N
           if A(m,n)<P(n) && 0<A(m,n)
                P(n)=A(m,n);
                finitePathStart(n)=true;
                finitePathEnd(m)=true;
                endPoint(n)=SPs(m);
                finitePathPowers(n,1)=pathPowers(n);
                finitePathPowers(n,2)=pathPowers(m);
           end
        end
    end

end