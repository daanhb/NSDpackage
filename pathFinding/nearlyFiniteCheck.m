function [yn,P, nearCP] = nearlyFiniteCheck(p0, h0, h1, startCP, CPs)
%a condition to determine it a path is nearly finite

    jumpThresh=1;
    
    p=p0(2:end);
   
    N=length(h1);
    diff_h1=(abs(abs(h1(1:(N-1)))-abs(h1(2:N)))) ./ (p(2:N)-p(1:(N-1)));
    if max(diff_h1)>jumpThresh
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