function dHdp = NSDpathODE(p,h,n,g,ICs,ascFlag)
    %computes NSD path h and h', at point order n

    AscDescDir=1; %
    if nargin == 6
        if ascFlag
            AscDescDir=-1;
        end
    end
    
    if length(g)<(n+1)
        error('Higher derivatives of g required as input');
    end
    pSmallThresh=1E-12;
    if n==0
         dHdp =  AscDescDir*1i./g{2}(h); 
    elseif p>pSmallThresh && abs(g{2}(h(1)))>pSmallThresh
        switch n
            case 1
                dHdp = [h(2); (AscDescDir*2i-h(2).^2.*(g{3}(h(1))))./g{2}(h(1))]; 
                % H=[h,h']

            case 2
                dHdp = [h(2); h(3); (AscDescDir*6i-h(2).^3.*(g{4}(h(1)))-3*g{3}(h(1)).*h(2).*h(3))./g{2}(h(1))];
                % H=[h,h',h'']
            otherwise
                error('Cant go above 2nd order stationary points yet');
        end
    else
        %very small p, can get unstable... so approximate Taylor style:
        dHdp = [ICs(2); zeros(n,1)];
%         switch n 
%             case 1
%                 [AscDescDir*2i/(g{3}*ICs(1)*ICs(2))];
%             case 2
%                 [h(2); h(3); (AscDescDir*6i-h(2).^3.*(g{4}(h(1)))-3*g{3}(h(1)).*h(2).*h(3))./(g{3}())];
%         end
    end
end