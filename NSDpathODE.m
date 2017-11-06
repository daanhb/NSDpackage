function dhdp = NSDpathODE(p,h,n,g)
    %computes NSD path h and h', at point order n

    if length(g)<(n+2)
        error('Higher derivatives of g required as input');
    end
    switch n
        case 0
            dhdp =  1i./g{2}(h);
        case 1
            dhdp = [h(2); (2i-h(2).^2.*(g{3}(h(1))))./g{2}(h(1))];
            %not sure if this longer equation is necessary... could instead
            %just write:
            %dhdp=1i*r*p^(r-1)./g{2}(h) for a general ODE
        otherwise
            error('Cant go above 1st order stationary points yet');
                
    end
    
end