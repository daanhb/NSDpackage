function dHdp = NSDpathODE(p,h,n,g,ICs,ascFlag)
    %computes NSD path h and h', at point order n

    AscDescDir=1; %
    if nargin == 6
        if ascFlag
            AscDescDir=-1;
        end
    end
    
    if length(g)<(n+2)
        error('Higher derivatives of g required as input');
    end
    pSmallThresh=1E-12;
    if n==0
         dHdp =  AscDescDir*1i./g{2}(h); 
    elseif p>=pSmallThresh && abs(g{2}(h(1)))>=pSmallThresh
        switch n
            case 1
                dHdp = [h(2); (AscDescDir*2i-h(2).^2.*(g{3}(h(1))))./g{2}(h(1))]; 
                % H=[h,h']

            case 2
                dHdp = [h(2); h(3); (AscDescDir*6i-h(2).^3.*(g{4}(h(1)))-3*g{3}(h(1)).*h(2).*h(3))./g{2}(h(1))];
                % H=[h,h',h'']
            case 3
                dHdp = [h(2); h(3); h(4); (AscDescDir*24i - g{5}(h(1)).*h(2).^4 - 6*g{4}(h(1)).*h(2).^2.*h(3) - g{3}(h(1)).*(3*h(3).^2 + 4*h(4).*h(2)))./g{2}(h(1))];
                % H=[h,h',h'', h''']
            case 4
                dHdp = [h(2); h(3); h(4); h(5); ...
                        (AscDescDir*120i - g{6}(h(1)).*h(2).^5 - 10*g{5}(h(1)).*h(2).^3.*h(3) - 10*h(4).*g{4}(h(1)).*h(2).^2 - h(3).*g{3}(h(1)) - 5*h(2).*(3*g{4}(h(1))).*h(3).^2 + h(5).*g{3}(h(1)))./g{2}(h(1));
                        ];
                
            otherwise
                error('Cant go above 3rd order stationary points yet');
        end
    else
        %very small p, can get unstable... so approximate Taylor style:
        %dHdp = [ICs(2); zeros(n,1)];
        dHdp = [ICs(2:(n+1)).'; 0];
        % ***** should go to higher order Taylor approximation at some point,
        % there's not really anything stopping me with NSD_ICsV3
%         dHdp(1) = ICs(2);
%         for n_ = 1:n
%             dHdp(n_+1)=ICs(2)
%         switch n 
%             case 1
%                 [AscDescDir*2i/(g{3}*ICs(1)*ICs(2))];
%             case 2
%                 [h(2); h(3); (AscDescDir*6i-h(2).^3.*(g{4}(h(1)))-3*g{3}(h(1)).*h(2).*h(3))./(g{3}())];
%         end
    end
end