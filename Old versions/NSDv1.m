function [ X, W ] = NSDv2( freq,a,b,N,g,singularities )
    %generic Numerical steepest descent routine
    %singularities are all along the real line, and specified by the user,
    %input in vector form
    if nargin<6
        error('Cant automatically determine inverse / derivatives / singularities of phase function yet')    
    end
    
    %error in comparing limits of arguemnt of paths at infinity
    errTol=1E-14;
    
    %determine the number of steepest descent paths
    numPaths=2*length(g.stationaryPoints)+2;
    
    realPoints=[a repmat([g.stationaryPoints],1,2) b];
    
    pathPowers=[1 repmat(g.order+1,1,2) 1];
        
    X=[];   W=[];
        for path=[1 numPaths 2:(numPaths-1)]
            if iscell(g.inverse) %multi-valued inverse, stored in cells
                   %choose the correct branch here, the one that makes sense.
                   %start with the first branch:
                   branch=1;
                   h{path}=@(p) g.inverse{1}(g.eval(realPoints(path))+1i*p.^pathPowers(path));
                    if path == 1 || path == numPaths
                        %check the SD path passes through the endpoints of the
                        %original integral:
                       while abs(h{path}(0)-realPoints(path))>errTol
                           %try a different branch
                           branch=branch+1;
                           h{path}=@(p) g.inverse{branch}(g.eval(realPoints(path))+1i*p.^pathPowers(path));
                       end                
                    else
                        %now compute the limit of the argument of the path as it tends
                        %to inifinty, so the next path can be chosen correctly
                        if path == 2
                            pathEndArg=angle(h{1}(inf));
                        elseif path == (numPaths-1)
                            pathEndArg=angle(h{numPaths}(inf));
                        else
                            error('Not sure how to handle this many stationary points');
                        end
                        %check that the SD paths join up at infinity:
                        while abs(pathEndArg-angle(h{path}(inf)))>errTol
                           branch=branch+1;
                           h{path}=@(p) g.inverse{branch}(g.eval(realPoints(path))+1i*p.^pathPowers(path));
                        end
                    end
                    %function is invertible, so can define SD paths
                    %can write the derivative of the inverse of g using that a-level trick
                    %that I forgot about
                    DerivOfInverseOf_g=@(x) 1./(g.inverse{branch}(g.deriv(x)));
                    dh_dp{path}=@(p) 1i*DerivOfInverseOf_g(g.eval(realPoints(path))+1i.*p.^pathPowers(path));

            else
                h{path}=@(p) g.inverse(g.eval(realPoints(path))+1i*p.^pathPowers(path));
                DerivOfInverseOf_g=@(x) 1./(g.inverse(g.deriv(x)));
                dh_dp{path}=@(p) 1i*DerivOfInverseOf_g(g.eval(realPoints(path))+1i.*p.^pathPowers(path));
            end
            
            %change position of singularities with change of variables
            COVsingularities=singularities;
            for s=1:length(singularities)
                COVsingularities(s).position=(singularities(s).position-realPoints(path))*freq/1i;
            end
            % get weights and nodes for path.
            [x{path}, w{path}]=pathQuad( 0, inf, pathPowers(path), COVsingularities, N );
            %now append weights and nodes to the rest, noting that we subtract every even path, as these are coming back down from
            %infinity
            W=[W; ((-1)^(path+1)/(freq^(1/pathPowers(path))))*exp(1i*freq*g.eval(realPoints(path))).*dh_dp{path}(x{path}./freq).*w{path};];
            %nodes
            X=[X; h{path}(x{path}./freq^(1/pathPowers(path)));];

        end

end

