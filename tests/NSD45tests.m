clear;

Npts=15;

%add functions with branch cuts

%do timing comparisons against brute force approx, with/without specifying
%stationary points

%need to test singular f(x) too

%need to include higher derivatives atm
G{1}={ @(x)x.^3-.001*x  ,  @(x)3*x.^2-.001  ,  @(x)6*x  ,  @(x)6};
G{2}={ @(x)x.^3-.0000001*x  ,  @(x)3*x.^2-.0000001   , @(x)6*x ,   @(x)6};
G{3}={@(x)x.^3+.1*x ,   @(x)3*x.^2+.1  ,  @(x)6*x  ,  @(x)6};
G{4}=G{1};

f=@(x) 1; %it doesn't really matter what f(x) is yet

%exact to 5 digits, that's all wolfram alpha gives you for free :/
exactVal=[ 0.16246 
           0.15522 
           0.0676662 
           0.0437498 ];
       
freq=[1000,1000,100,100000,1000];
       
 for j=1:length(exactVal)
    tic;
    [ X, w ] = NSD45( -1,1,freq(j),Npts,G{j}); 
    T=toc;
    if exactVal(j)~=0
        err(j)=abs(sum(w)-exactVal(j))/exactVal(j);
    else %do absolute error in this case
        err(j)=abs(sum(w)-exactVal(j));
    end
    fprintf('Case %d error: %e\t took %.2f s\n,',j, err(j), T)
    clear X w;
    jCount=j;
 end
 
 
 %return;
 %-----------------------------------------------------------------------%
 %first one doesn't work so well, instead try specific values:
 jCount=jCount+1;
 SPs = [-sqrt(.001/3)  sqrt(.001/3)];
 Os  = [1 1];
 
 [ X, w ] = NSD45( -1,1,freq(1),Npts,G{1},'stationary points',SPs,'order',Os, 'poles',[]);
 
  err(jCount)=abs(sum(w)-exactVal(1))/exactVal(1);
  fprintf('Case %d error: %e\n',jCount, err(jCount));
  
  %even including the correct stationary points fails, and makes very little
  %differece to the error.
  
  %throwing away the middle two vertical integrals increases the error to
  %10%, they clearly count for a lot
  
  %It seems that this fails because there are two oscillations
  %at once here, and they are at very different rates. Increasing the
  %frequency (see test case four) seems to 'fix' things. Ulitmately, need a
  %test and an alternative method for such cases.
  
 jCount=jCount+1;
 %try absorbing slow oscillations into f, see if this helps
 G{jCount}={@(x) x.^3, @(x) 3*x.^2, @(x) 6*x, @(x) 6};
 [ X, w ] = NSD45( -1,1,freq(1),Npts,G{jCount});
 F=@(x) f(x).*exp(1i*freq(1)*-.001*x);
 
 err(jCount)=abs(w.'*F(X)-exactVal(1))/exactVal(1);
 fprintf('Case %d error: %e\n',jCount, err(jCount));
 % Yayyyy this fixes things. However, implementing this in a black box type
 % scheme... that's gonna be tricky, or just impossible.
 
 %---------------------------------------------------------------------------
 jCount=jCount+1;
 
 %Simon says this one is tricky:
 
 G{jCount}={@(x) (x+1).^-1, @(x) -(x+1).^-2, @(x) 2*(x+1).^-3, @(x) -6*(x+1).^-4};
 %use 100 points just to check how well it hugs the path
 [ X, w ] = NSD45( 0,1,1000,25,G{jCount},'stationary points',[], 'poles',[], 'visuals on');
 SWCexact=0.002682651243826470 - 0.004106825681982276i;
 
 err(jCount)=abs(sum(w)-SWCexact)/abs(SWCexact);
 fprintf('Case %d, (Simon special) error: %e\n', jCount, err(jCount));
 
 
 %---------------------------------------------------------------------------
jCount=jCount+1;
%now try some square root tests:
a=.5; b=1;

%G{jCount}={@(x) (x.^2-.5).^.5+x, @(x) x./sqrt(x.^2 - 0.5) + 1 , @(x) -0.5./(x.^2 - 0.5).^(3/2) };
[G{jCount}, SPs, SPOs ] = SqrtTypePhase( a,b );

% %get exact value by brute force:
% % [xBruce, wBruce] = Bruce( -1,1,20,freq(5));
% % exactValue=exp(1i*freq(5)*G{jCount}{1}(xBruce)).'*wBruce;

% ???for some reason brute force computations would not converge???
exactValue=0.00567695 + 0.000130121i;

[ X, w ] = NSD45( -1,1,1000,Npts,G{jCount},'stationary points',SPs,'order',SPOs, 'poles',[]);
err(jCount)=abs(sum(w)-exactValue)/abs(exactValue);
fprintf('Case %d, (HNA-type) error: %e\n', jCount, err(jCount));