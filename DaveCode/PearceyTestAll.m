%test a Pearchey integral with polynomial phase, for Dave
     
clc; %clear screen
clear classes; %clear all variables
close; %close any images currently open, so we don't draw all over them

Npts=30; %number of quadrature points per SD path
X=[-8 -6 -4 -2 2 4 6 8]; Y=[2 4 6 8]; %}
a=-1; b=1; %directions of valleys
order=4;
freq=1;
mkdir('PearceyImages')

xCount=0;
for x=X
    xCount=xCount+1;
    fprintf('x=%d\n',x);
    yCount=0;   
    for y=Y
        yCount=yCount+1;
        fprintf('\ty=%d ',y);
        polyCoeffs = [1 0 x y 0]; %a_1*x^N + ... a_N*x + a_{N+1}, if a:=polyCoeffs
        %% ----------------------------------------------------------------------------------------- %%
        %first derivative
        D1polyCoeffs = polyCoeffs(1:order).*fliplr(1:order);
        %second derivative
        D2polyCoeffs = D1polyCoeffs(1:(order-1)).*fliplr(1:(order-1));
        %third derivative
        D3polyCoeffs = D2polyCoeffs(1:(order-1)).*fliplr(1:(order-1));

        %define g and it's first three derivatives in format required:
        G = {@(x) polyval(polyCoeffs,x),@(x) polyval(D1polyCoeffs,x),@(x) polyval(D2polyCoeffs,x), @(x) polyval(D3polyCoeffs,x)};

        %% ----------------------------------------------------------------------------------------- %%

        %this radius could be considered artificial infinity, outside of which all SD
        %paths won't deviate much from a straight line
        R = monomialSettleRadius(polyCoeffs)+rand;

        [ Z, W ] = NSD45( a, b, freq, Npts, G, 'analytic', true, 'visuals on', ...
                            'settleRad', R, 'ainf', 'binf');
        title(sprintf('x=%d, y=%d', x, y));
        %imPathName=sprintf('PearceyImages\\x%dy%d',x,y);
        imPathNamJHPG=sprintf('PearceyImages\\x%dy%d.jpg',x,y);
%        savefig(imPathName);
        saveas(gcf,imPathNamJHPG)
        close;
        I_GHH=sum(W); %Gibbs-Hewett-Huybrechs estimate
        I_CHK=KirkPearceyData(x,y); %Conor-Hobbs-Kirk estimate
        err(xCount,yCount)=abs(I_CHK-I_GHH)/abs(I_CHK);
        fprintf('\tabs err=%e\n',err(xCount,yCount));
        %fprintf('\nRelative error (against 6dp of Pearcey data) %e',abs(I_CHK-I_GHH)/abs(I_CHK));
    end
end
%%plot some final error estimates
imagesc(1:8,1:4,log10(err).');
xticks=1:8;
xticklabels({'-8','-6','-4','-2','2','4','6','8'});
yticks=1:4;
yticklabels({'','2','','4','','6','','8'});
%view(0, 90);
xlabel('x'); ylabel('y'); title('(log10) error compared with CHK 6.d.p results');
colorbar;
saveas(gcf,'PearceyImages\\ErrorChecks.jpg');