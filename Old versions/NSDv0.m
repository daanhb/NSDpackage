function [ X, W ] = NSD( freq,a,b,N,g,singularities )
    %generic Numerical steepest descent routine
    %singularities are all along the real line, and specified by the user,
    %input in vector form
    if nargin<6
        error('Cant automatically determine inverse / derivatives / singularities of phase function yet')    
    end

    %parameter for 'being near to a singularity':
    Rs=.1;
    
    %parameters for the graded quadrature:
    p_max=10;   GradDelta=.15;
    
    if length(g.deriv)==1
        %function is invertible, so can define SD paths
        h=@(p,q) g.inverse(g.eval(q)+1i*p);
        dh_dp=@(p,q) 1i*g.derivOfInverse(g.eval(q)+1i*p);
        
        %depending on the locations of the sinularities, choose suitable
        %weights and nodes along steepest descent paths:
        
        endPoint=[a b];
        for path=[1 2]
            if isequal(singularities,[])
                [x{path}, w{path}] = GaussLaguerre(N, 0);
            else
                for s=singularities
                    if endPoint(path)==s.position
                        %endpoint is on singularity
                        if strcmp(s.blowUpType,'log')
                            [x{path}, w{path}]=quad_gengauss_loglaguerre(N);
                        else
                            error('have only coded for logarithimc singularities so far');
                        end
                    elseif abs(endPoint(path)-s.position)<=Rs %nearly singular
                        pathSplit=sqrt(Rs^2-(s.position-endPoint(path))^2);
                        [x1, w1]= GradedQuad( N, p_max, GradDelta );
                        %now scale the graded weights and nodes from [0 1] to [0 pathSplit]
                        x1=x1*pathSplit; w1=w1*pathSplit.*exp(-x1);
                        %now do Gauss Laguerre on [pathSplit \infty]
                        [x2, w2] = GaussLaguerre(N, 0);
                        x2=x2+pathSplit; w2=w2*exp(-pathSplit);
                        %combine all weights and nodes for first SD path
                        x{path}=[x1; x2;]; w{path}=[w1; w2];
                    else    %not even nearly singular
                        [x{path}, w{path}] = GaussLaguerre(N, 0);
                    end
                    clear x1 x2 w1 w2;
                end
            end
        end
        
        %weights
        W=(1/freq)*[exp(1i*freq*g.eval(a)).*dh_dp(x{1}./freq,a)*w{1};  -exp(1i*freq*g.eval(b)).*dh_dp(x{2}./freq,b)*w{2};];

        %nodes
        X=[h(x{1}./freq,a); h(x{2}./freq,b)];
    else
        error('Cant do anything harder than linear phase yet');
    end

end

