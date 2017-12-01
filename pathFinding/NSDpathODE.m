function dHdp = NSDpathODE(p,h,n,g,ICs)
    %computes NSD path h and h', at point order n

    if length(g)<(n+1)
        error('Higher derivatives of g required as input');
    end
    pSmallThresh=1E-12;
    if n==0
         dHdp =  1i./g{2}(h); 
    elseif p>pSmallThresh
        switch n
            case 1
                dHdp = [h(2); (2i-h(2).^2.*(g{3}(h(1))))./g{2}(h(1))]; 
                % H=[h,h']

                %not sure if this longer equation is necessary... could instead
                %just try writing:
                %dHdp=1i*r*p^(r-1)./g{2}(h) for a general ODE
            case 2
                dHdp = [h(2); h(3); (6i-h(2).^3.*(g{4}(h(1)))-3*g{3}(h(1)).*h(2).*h(3))./g{2}(h(1))];
                % H=[h,h',h'']
            otherwise
                error('Cant go above 2nd order stationary points yet');
        end
    else
        %very small p, can get unstable... so approximate Taylor style:
        dHdp = [ICs(2); zeros(n,1)];
    end
end