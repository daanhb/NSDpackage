function isZero = isZeroOnLine( a,b,F, thresh2, lineLength )

    if nargin<=4
        lineLength=abs(b-a);
    end

    thresh1=.01; %threshold for how close the complex zero has to be to the real one
    %thresh2=0.01; %threshold for how how close a zero can be to the line

    lineMap=@(s) a+exp(angle(b-a)*1i)*lineLength.*s;
    
    f=@(s) F(lineMap(s));
    
    t=linspace(0,1,1000*ceil(max(1,lineLength)));
    
    if min(abs(f(t)))<thresh2
        isZero=true;
%         if isZero
%             fprintf('zeroOnLine\n');
%         end
    else
        isZero=false;
    end
    return;
    
    
    % --- replaced the below with the above
    
    % --- it's faster and more effective, despite seeming less
    % sopphisticated
    
    %check for zeros around shifted versions of f, incase it just almost
    %x-axis
    shifts = [-thresh2 0 thresh2];
    realZeros=[];
    imagZeros=[];
    id0Real=false;
    id0Imag=false;
    j=0;
    for zeroShiftReal = shifts
        for zeroShiftImag = shifts
            j=j+1;
            [realZeros_{j}, id0Real_{j}] = AllZeros(@(x) real(f(x)+zeroShiftReal),0,1,1000);
            [imagZeros_{j}, id0Imag_{j}] = AllZeros(@(x) imag(f(x)+zeroShiftImag),0,1,1000);
            realZeros=[realZeros realZeros_{j}];
            imagZeros=[imagZeros imagZeros_{j}];
            id0Real=max(id0Real, id0Real_{j});
            id0Imag=max(id0Imag, id0Imag_{j});
        end
    end
    
    if ~id0Real && ~id0Imag
        [X,Y] = meshgrid(realZeros,imagZeros);
        zeroDists = min(min(abs(X-Y)));
        if min(min(zeroDists)) < thresh1
            isZero = true;
        else
            isZero = false;
        end
    elseif ~id0Real && ~isempty(realZeros)
        isZero = true;
    elseif ~id0Imag && ~isempty(imagZeros)
        isZero = true;
    elseif id0Imag && id0Real
        error('Function seems to be identically zero');
    else
        isZero = false;
    end
%     
%     if isZero
%         fprintf('zeroOnLine\n');
%     end
%     
end