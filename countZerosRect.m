function allZeros = countZerosRect( f,df, rect, N )
%detects roots of analytic function inside of a given rectangle in complex
%domain
%N is number of points per unit length of rectangle
    if nargin<=4
        N=15;
    end
    
    if length(rect)~=4
        error('rectangle must consist of four complex values, ordered anti-clockwise');
    end
    
    rect(5)=rect(1); %add 5th value as first to complete loop
    
    %two things to integrate:
    allZeros=0;
    
    for j=1:4
        width=abs(rect(j+1)-rect(j));
        quadPts=max(ceil(N*width),5);
        [z_, w] = quad_gauss(quadPts, 0, width);
        %want unit direction of contour along edge
        dir=round((rect(j+1)-rect(j))/width);
        z = rect(j) + dir*z_;
        
        allZeros=allZeros+w.'*dir*(df(z)./f(z))/(2*pi*1i);
        
    end
    
    %will always be an integer, so...
    allZeros=round(allZeros);
    
end