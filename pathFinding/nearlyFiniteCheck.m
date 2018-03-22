function [yn,P, nearCP] = nearlyFiniteCheck(p0, h0, h1, startCP, CPs, freq)
%a condition to determine it a path is nearly finite

    jumpThresh=1;
    distThresh=1/freq;
    
    p=p0(2:end);
   
    N=length(h1);
    otherCPs = CPs(startCP~=CPs);
    diff_h1=(abs(abs(h1(1:(N-1)))-abs(h1(2:N)))) ./ (p(2:N)-p(1:(N-1)));
    dist_h1=abs(otherCPs-h0);
    
    if max(diff_h1)>jumpThresh && min(min(dist_h1))<distThresh
        yn=true;
        [~,pIndex] = max(diff_h1);
        P=p(pIndex);
        dist=inf;
        for CPindex=1:length(CPs)
            %only check other critical points
            if ~ismember(startCP,CPs(CPindex))
                dist_=abs(h0(pIndex)-CPs(CPindex));
                if dist_<dist
                    nearCP=CPindex;
                    dist=dist_;
                end
            end
        end
    else
        yn=false;
        P=[];
        nearCP=[];
    end
    
end