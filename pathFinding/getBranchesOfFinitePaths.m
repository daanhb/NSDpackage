function [h,dhdp, W, FPindices] = getBranchesOfFinitePaths(prePathLengthsVec, pathEndIndex, criticalPoints, pathPowers, G, freq, N, COVsingularities, RelTol)

 
FPindices=false(length(prePathLengthsVec),2);

    %initialise output
    for j=1:length(prePathLengthsVec)
        for m=1:pathPowers(j)
            h{j,m}=[];
            dhdp{j,m}=[];
            W{j,m}=[];
        end
    end
    for j=1:length(prePathLengthsVec)
       if ~isinf(prePathLengthsVec(j))
          %there's a finite path here
          P1=(prePathLengthsVec(j)/2)^(1/pathPowers(j));
          ICsSD = NSDpathICv3( pathPowers(j), G, criticalPoints(j), false);
          [Psd{j}, Wsd{j}] = pathQuadFinite( P1, COVsingularities, freq, N );
          for m=1:pathPowers(j)
               [~,Hsd{m}] = ode45(@(t,y) NSDpathODE(t,y,pathPowers(j)-1,G, ICsSD{m}, false), Psd{j}, ICsSD{m}, odeset('RelTol',RelTol) );      
          end
          
          P2=(prePathLengthsVec(j)-P1^(pathPowers(j)))^(1/pathPowers(pathEndIndex(j)));
          ICsSA = NSDpathICv3( pathPowers(pathEndIndex(j)), G, criticalPoints(pathEndIndex(j)), true);
          [Psa{j}, Wsa{j}] = pathQuadFinite( P2, COVsingularities, freq, N );
          for n=1:pathPowers(pathEndIndex(j))
              [~,Hsa{n}] = ode45(@(t,y) NSDpathODE(t,y,pathPowers(pathEndIndex(j))-1,G, ICsSA{n}, true), Psa{j}, ICsSA{n}, odeset('RelTol',RelTol) );
          end
       
       %now check which was the best match. See which midpoints of which
       %branches end up the closest
       dist=inf;
           for  m=1:pathPowers(j)
               for n=1:pathPowers(pathEndIndex(j))
                    dist_=abs(Hsd{m}(end,1)-Hsa{n}(end,1));
                    if dist_<dist
                        dist=dist_;
                        mClosest=m;
                        nClosest=n;
                    end
               end
           end
           
            h{j,mClosest}=[Hsd{mClosest}(:,1); Hsa{nClosest}(:,1);];
            dhdp{j,mClosest}=[Hsd{mClosest}(:,2); Hsa{nClosest}(:,2);];
            W{j,mClosest}=[Wsd{j}; Wsa{j};];
            FPindices(j,mClosest)=true;
           clear ICsSA Hsd Hsa Wsa Wsd dist;
       end
    end
end