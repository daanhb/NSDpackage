classdef Integral
    %all the information you could ever need about an integral
    
    properties
        dim
        integrand %anonymous function which takes 'dim' intputs and returns value of integrand
        integrandNonOsc
        weightType =[]
        domain
        width      %width of each interation domain, useful when withs are really tiny but numbers are large
        %--- now define more detailed properties for sophisticated
        %integration routines ----------------------------------%
        singularity % a vector, each containing info about a different singularity
        oscillator % info about the oscillatory nature of the integrand
        analytic
    end
    
    methods
        function self = Integral(dim, integrand, domain, width, analytic, singularityIn, oscillatorIn, integrandNonOsc)
            if nargin==5
                self.integrandNonOsc=[];
                oscillatorIn=[];
            end
            if ~isequal(oscillatorIn,[]) && ~isa(oscillatorIn,'oscillator')
                error('6th input to integral.m constructor must be of oscillator.m type');
            end
            if isequal(integrandNonOsc,[])
                self.integrandNonOsc=@(varargin) self.integrand.*exp(-1i*oscillatorIn(varargin));
            else
                self.integrandNonOsc=integrandNonOsc;
            end
            if ~isequal(singularityIn,[]) && ~isa(singularityIn,'singularity')
               error('4th input to integral.m constructor must be of singularity.m type') 
            end
            
            %set all of these values up
            self.dim=dim;
            self.domain=domain;
            if isequal(width,[])
                width=zeros(dim,1);
                for n=1:dim
                    width(n)=domain{n}(2)-domain{n}(1);
                end
                        
            else
                self.width=width;
            end
            
            %now make integrand variable number of inputs:
           integrandText='self.integrand = @(varargin) integrand(varargin{1}';
           for n=2:dim
               integrandText=strcat(integrandText, sprintf(',varargin{%d}',n));
           end
           integrandText=strcat(integrandText, ');');
           eval(integrandText);
            
            self.singularity=singularityIn;
            self.oscillator=oscillatorIn;            
        end
    end
    
end

