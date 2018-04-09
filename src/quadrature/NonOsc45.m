function [x, w] = NonOsc45(a,b,freq,N,g,singularity, width)
%alternative to NSD45 when the integrand is sufficiently non-oscillatory
    if nargin ==6
        width= b-a;
    end
    maxSingDist=.2;
    p_max=12;
    delta=.15;
    subFlag = false;
    
    Npts=max(ceil(abs(b-a)*N),N);
    if isempty(singularity)
        [x, w0] =gauss_quad(a,b,Npts);
    elseif length(singularity)==1
        if strcmp(singularity.blowUpType,'log')
            if singularity.position == a
                [x, w0] = genGaussLog( Npts, a, b, width, 'L');
            elseif singularity.position == b
                [x, w0] = genGaussLog( Npts, a, b, width, 'R');
            elseif a < singularity.position && singularity.position < b
                %split the integral and call recursively
                [x1, w1] = NonOsc45(a,singularity.position,freq,Npts,g,singularity);
                [x2, w2] = NonOsc45(singularity.position,b,freq,Npts,g,singularity);
                x = [x1; x2];
                w0 = [w1; w2];
                subFlag = true;
            else
                %could be a near singularity, check for this:
                [singDist, nearestEndpoint] = min([abs(a-singularity.position)  abs(b-singularity.position)]);
                if singDist<maxSingDist
                    [x1, w1]=GradedQuad( Npts,p_max,delta );
                    w0 = width*w1;
                    if nearestEndpoint==1
                        %close to a
                        x = a + width*x1;
                    else
                        %close to b
                        x = b - width*x1;
                    end
                else
                    [x, w0] =gauss_quad(a,b,Npts);
                end
            end
        else
            error('Cannot handle any singularities other than a log yet');
        end
    else
        error('Cannot do standard quadrature for multiple singularities yet');
    end
    
    %only scale by oscillator on the top level of the stack
    if ~subFlag
        w = w0.*exp(1i*freq*g(x));
    else
        w = w0;
    end
end

