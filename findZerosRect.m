function allZeros = findZerosRect( f, df, rectIn, thresh, N )
%detects roots of analytic function inside of a given rectangle in complex
%domain
%N is number of points per unit length of rectangle
%thresh is the width of the rectangle which we pinpoint the zeros to

    if nargin<=4
        N=15;
    end
    
    if length(rectIn)~=4
        error('rectangle must consist of four complex values, ordered anti-clockwise');
    end
    
    
    height=max(imag(rectIn))-min(imag(rectIn));
    width=max(real(rectIn))-min(real(rectIn));
    
    
    %first check if zeros lie on edge of rectangle
    %...
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
    
    allZeros=rand+rand*1i;
    rect=rectIn;
    
    while max(height,width)>thresh % max(abs(f(allZeros)))>thresh %was an alternative
        
        BLcorner=min(real(rect))+1i*min(imag(rect));

        bisectRatio=.5; %start by bisecting in the middle
        
        zeroOnLine=true; %this value keeps choosing new bisection point until bisection line doesn't land on a zero
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
            
            zeroOnLine=~isempty(isZeroOnLine( fullRectIn(j), fullRectIn(j+1), f ));
            bisectRatio=.5+.5*(rand-.5);% in (.25,.75)
        end
        
        %count the number of zeros in each subdivision:
        Acount=countZerosRect( f, df, newRectA, N );
        Bcount=countZerosRect( f, df, newRectB, N );
        
        if Acount>0 && Bcount>0
            %call self recursively over each sub-tangle, until zeros are
            %located
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