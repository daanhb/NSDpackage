function [hLevel]=NSDlevelCurves(ha0,hb0, G, freq)
% 
% %ODE for the level curves of SD paths
%     dHdp=@(sigma, h_) G{2}(sigma)./G{2}(h_);
% %     res=@(ha, hb) [ha(1)-haP; hb(1)-hbP];
% %     solinit = bvpinit(linspace(-1,1,100),[haP hbP]);
%  %   res=@(ha) ha-haP;
%    % solinit = bvpinit(linspace(-1,1,100),haP);
%     
%     [~,hLev] = ode45(dHdp,[-1,1], haP);

%^^ doesn't work so well ^^


%     hold on;
%     for zeta=linspace(-1,1,10)
%         [~,h_] = ode45(@(t,y) NSDpathODE(t,y,0,G, zeta), P0, zeta, odeset('RelTol',RelTol) ); plot(h_,'x');
%     end
% ^^ also doesn't work so well ^^

%try a Newton iteration:
    Npts=100;

    gThresh=log(1E-15)/freq;
    %NewtonThresh=0;%1E-12;
    
    F    = @(x,y) imag(G{1}(x+1i*y)) - gThresh;
    dFdy = @(x,y) real(G{2}(x+1i*y));
    
    x=linspace(real(ha0),real(hb0),Npts);
    y=[imag(ha0) zeros(1,Npts-2) imag(hb0)];
    
    for j=2:length(y)
       y(j)=y(j-1); %first guess at previous value
       
       for n=1:20 %F(x(j),y(j)) > 1 
           y(j) = y(j) - F(x(j),y(j))/dFdy(x(j),y(j));
       end       
        
    end
    
    hLevel=x+1i*y;
    
end

