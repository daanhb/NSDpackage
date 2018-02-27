function [allZeros, orders] = findZerosRect( f, df, rectIn, thresh, N, visuals)
%detects roots of analytic function inside of a given rectangle in complex
%domain
%N is number of points per unit length of rectangle
%thresh is the width of the rectangle which we pinpoint the zeros to

    randScale=0.5; %needs to be sufficiently bigger than thresh2 of isZeroOnLine
    if nargin<=4
        N=15;
    end
    
    if nargin <=5
        visuals=false;
    end
    
    if length(rectIn)~=4
        error('rectangle must consist of four complex values, ordered anti-clockwise');
    end
    
    intThresh=0.01; %copied from NSD45 codez
    
    
    height=max(imag(rectIn))-min(imag(rectIn));
    width=max(real(rectIn))-min(real(rectIn));
    
    
    %first check if zeros lie on edge of rectangle
    %...
    fullRectIn=[rectIn rectIn(1)];
    zeroOnLine=false;
    for j=1:4
        zeroOnLine=max(zeroOnLine,isZeroOnLine( fullRectIn(j), fullRectIn(j+1), f ));
    end
    failCount=0;
    while zeroOnLine
        if failCount>10
            error('Keep failing to find a line without any zeros');
        end
        %retry with a larger rectangle
        newRect=mean(rectIn) + (1+rand*randScale)*[-width-1i*height  width-1i*height  width+1i*height  -width+1i*height]/2;
        %[allZeros, orders] = findZerosRect( f, df, newRect, thresh, N );
        fullRectIn=[rectIn rectIn(1)];  
        zeroOnLine=false;
        for j=1:4
            zeroOnLine=max(zeroOnLine,isZeroOnLine( fullRectIn(j), fullRectIn(j+1), f ));
        end
        failCount=failCount+1;
    end
%     
    fprintf('\t%d\n',failCount);
        
    allZeros=rand+rand*1i;
    rect=rectIn;
    
    while max(height,width)>thresh % max(abs(f(allZeros)))>thresh %was an alternative
        
        BLcorner=min(real(rect))+1i*min(imag(rect));

        bisectRatio=.5; %start by bisecting in the middle
        
        zeroOnLine = true; %this value keeps choosing new bisection point until bisection line doesn't land on a zero
        failCount=0;
        while zeroOnLine
            if failCount>10
                error('Keep failing to find a line without any zeros');
            end
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
%             
            zeroOnLine=isZeroOnLine( bisectingLine(1), bisectingLine(2), f );
            bisectRatio=.5+randScale*.5*(rand-.5);% in randScale*(.25,.75)
            failCount=failCount+1;
        end
        fprintf('\t%d\n',failCount);
        %plot the new rectangles
        if visuals
            plotRects( newRectA, newRectB );
        end
                
        %count the number of zeros in each subdivision:
        countOutA = countZerosRect( f, df, newRectA, N, 2 );
        Acount = round(countOutA(1)); Apos=countOutA(2);
        countOutB = countZerosRect( f, df, newRectB, N, 2 );
        Bcount = round(countOutB(1)); Bpos=countOutB(2);
        
        if Acount>0 && Bcount>0
            %call self recursively over each sub-tangle, until zeros are
            %located
            [allZerosA, ordersA] = findZerosRect( f, df, newRectA, thresh, N, visuals );
            [allZerosB, ordersB] = findZerosRect( f, df, newRectB, thresh, N, visuals );
            allZeros = [allZerosA allZerosB];
            orders=round([ordersA ordersB]);
            
            return;
        elseif Acount>0
            %first check if this is the only zero
            if Acount==1
                allZeros=Apos;
                orders=1;
                return
            else
                %guess at the weighted average, if there's only one
                %stationary point in this rect, this is where it'll be
                maybeCentreA = Apos/Acount;
                focusRectA = [maybeCentreA-thresh-thresh*1i  maybeCentreA+thresh-thresh*1i  maybeCentreA+thresh+thresh*1i  maybeCentreA-thresh+thresh*1i];
                focusAcount = round(countZerosRect( f, df, focusRectA, N ));
                if focusAcount==Acount
                    allZeros=mean(focusRectA);
                    orders=Acount;
                    return;
                end
            end
            rect=newRectA;
        elseif Bcount>0
             if Bcount==1
                allZeros=Bpos;
                orders=1;
                return
             else
                %guess at the weighted average, if there's only one
                %stationary point in this rect, this is where it'll be
                maybeCentreB = Bpos/Bcount;
                focusRectB = [maybeCentreB-thresh-thresh*1i  maybeCentreB+thresh-thresh*1i  maybeCentreB+thresh+thresh*1i  maybeCentreB-thresh+thresh*1i];
                focusBcount = round(countZerosRect( f, df, focusRectB, N ));
                if focusBcount==Bcount
                    allZeros=mean( focusRectB );
                    orders=Bcount;
                    return;
                end
            end
            rect=newRectB;
        elseif Acount==0 && Bcount==0
            allZeros=[];
            orders=[];
            break
        end
        
        height=max(imag(rect))-min(imag(rect));
        width=max(real(rect))-min(real(rect));
         
        allZeros=mean(rect);
        
    end
    
end