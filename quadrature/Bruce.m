function [x, w] = Bruce( a,b,Npts,freq, abWidth )
%brute force oscillatory integrator, used to validate NSD routines
    if nargin<5
        abWidth=b-a;
    end
    
    %first determine how many wavelengths are in interval
        wavelengths=ceil(abWidth/(2*pi/freq));
    %create some blank vectors
%         x=zeros(wavelengths*Npts,1);
%         w=x;
    %now fill them with composite Gaussian quadrature
    x=[];
    w=[];
    
    ChunkSize=abWidth/wavelengths;
        for j=1:wavelengths
            [x_, w_] = quad_gauss(Npts, a+(j-1)*ChunkSize, a+j*ChunkSize);
            x = [x; x_];
            w = [w; w_];
            
            %[x(((j-1)*Npts+1):(j*Npts)),w(((j-1)*Npts+1):(j*Npts))] = quad_gauss(Npts, a+(j-1)*abWidth/wavelengths, a+j*abWidth/wavelengths);
        end
        %now top it up to the end
%         [x(((wavelengths-1)*Npts+1):end),w(((wavelengths-1)*Npts+1):end)] = quad_gauss(Npts, a+(wavelengths-1)*abWidth/wavelengths, b);
end

