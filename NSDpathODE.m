function dHdp = NSDpathODE(p,h,n,g)
    %computes NSD path h and h', at point order n

    if length(g)<(n+2)
        error('Higher derivatives of g required as input');
    end
    switch n
        case 0
            dHdp =  1i./g{2}(h); 
            %H=h
            
        case 1
%             dHdp = [h(2); (2i-h(2).^2.*(g{3}(h(1))))./g{2}(h(1))]; 
%             %H=[h',h], aparantly... changed to:
            dHdp = [(2i-h(1).^2.*(g{3}(h(2))))./g{2}(h(2)); h(1)]; 
            %H=[h',h], hopefully...

            
            %not sure if this longer equation is necessary... could instead
            %just try writing:
            %dHdp=1i*r*p^(r-1)./g{2}(h) for a general ODE
        case 2
            dHdp = [h(2); h(3); (6i-h(2).^3.*(g{4}(h(1)))-3*g{3}(h(1)).*h(2).*h(3))./g{2}(h(1))];
            %H=[h,h',h'']
        otherwise
            error('Cant go above 2nd order stationary points yet');
    end
    
end