function allZeros = findZerosRect( f, df, rectIn, thresh, N )
%detects roots of analytic function inside of a given rectangle in complex
%domain
%N is number of points per unit length of rectangle
%thresh is the width of the rectangle which we pinpoint the zeros to

%need to account for multiple zeros, by calling self recursively

%need to check for zeros on edge of domain

Nedge=10000;

    if nargin<=4
        N=15;
    end
    
    if length(rectIn)~=4
        error('rectangle must consist of four complex values, ordered anti-clockwise');
    end
    
    %first check if zeros lie on edge of rectangle
    %...
    
    height=max(imag(rectIn))-min(imag(rectIn));
    width=max(real(rectIn))-min(real(rectIn));
    
    fullRectIn=[rectIn rectIn(1)];
    for j=1:4
%         Redge=linspace(fullRectIn(j),fullRectIn(j+1),Nedge);
%         [minVal, minInd]=min(abs(f(Redge)));
        zeroOnLine=~isempty(isZeroOnLine( fullRectIn(j), fullRectIn(j+1), f ));
        if zeroOnLine
            %retry with a larger rectangle
            newRect=mean(rectIn) + (1+rand/10)*[-width-1i*height  width-1i*height  width+1i*height  -width+1i*height]/2;
            allZeros = findZerosRect( f, df, newRect, thresh, N );
            return;
        end
    end
    
    allZeros=rand+rand*i;
    rect=rectIn;
    
    while max(abs(f(allZeros)))>thresh %max(height,width)>thresh
        
         BLcorner=min(real(rect))+1i*min(imag(rect));
%         if width>height
%             %take bototm left corner, and add half (width) rectangle to it
%             newRectA=BLcorner + [0  .5*width  (.5*width+1i*height)  1i*height];
%             newRectB=newRectA + .5*width;
%             bisectingLine=BLcorner+.5*width+[0  1i*height];
%         else
%             %take bototm left corner, and add half (height) rectangle to it
%             newRectA=BLcorner + [0  width  (width+.5i*height)  .5i*height];
%             newRectB=newRectA + .5i*height;
%             bisectingLine=BLcorner+.5i*height+[0  width];
%         end
% 
%         bisectingLinePoints=linspace(bisectingLine(1),bisectingLine(2),Nedge);
%         [minVal, minInd]=min(abs(f(edge)));


        bisectRatio=.5; %start by bisecting in the middle, but possibly change this if we land on a zero
        minVal=0;
        
        zeroOnLine=true;
        while zeroOnLine
            if width>height
                %take bototm left corner, and add half (width) rectangle to it
                newRectA=BLcorner + [0  bisectRatio*width  (bisectRatio*width+1i*height)  1i*height];
                newRectB=newRectA + bisectRatio*width;
                bisectingLine=BLcorner+bisectRatio*width+[0  1i*height];
            else
                %take bototm left corner, and add half (height) rectangle to it
                newRectA=BLcorner + [0  width  (width+bisectRatio*1i*height)  bisectRatio*1i*height];
                newRectB=newRectA + bisectRatio*1i*height;
                bisectingLine=BLcorner+bisectRatio*1i*height+[0  width];
            end
            
%             bisectingLinePoints=linspace(bisectingLine(1),bisectingLine(2),Nedge);
%             [minVal, minInd]=min(abs(f(Redge)));
            zeroOnLine=~isempty(isZeroOnLine( fullRectIn(j), fullRectIn(j+1), f ));
            bisectRatio=.5+.5*(rand-.5);% in (.25,.75)
        end
        
        
        Acount=countZerosRect( f, df, newRectA, N );
        Bcount=countZerosRect( f, df, newRectB, N );
        
        if Acount>0 && Bcount>0
            allZeros=[findZerosRect( f, df, newRectA, thresh, N )  findZerosRect( f, df, newRectB, thresh, N )];
            return;
        elseif Acount>0
            rect=newRectA;
        elseif Bcount>0
            rect=newRectB;
        elseif Acount==0 && Bcount==0
            allZeros=[];
            break
        end
        
        height=max(imag(rect))-min(imag(rect));
        width=max(real(rect))-min(real(rect));
        
        plotRects( newRectA, newRectB );
        allZeros=mean(rect);
        
    end
    
    
    
end