classdef oscillator
    %information about oscillatory components of integral
    
    properties
        dim
        freq
        phase
        phaseType
        phaseTypeData
        phaseDeriv
        phaseInverse
        stationaryPoints
    end
    
    methods
        function self=oscillator(dim,freq,phase,phaseType,phaseTypeData)
           self.dim=dim;
           self.freq=freq;
           self.phase=phase;
           if nargin==5
               self.phaseType=phaseType;
               self.phaseTypeData=phaseTypeData;
           else
               self.phaseType=[];
               self.phaseTypeData=[];
           end
           
           %store phase with variable number of input args, via this stupid
           %bodge:
           phaseText='self.phase = @(varargin) phase(varargin{1}';
           for n=2:dim
               phaseText=strcat(phaseText, sprintf(',varargin{%d}',n));
           end
           phaseText=strcat(phaseText, ');');
           eval(phaseText);
        end
        
        function oscOut=add(osc1,osc2)
            if osc1.freq ~= osc2.freq
               error('Cannot add two oscillatory of differnt frequencies, yet') ;
            else
                oscOut.freq=osc1.freq;
            end
            
            if osc1.dim ~= osc2.dim
                error('Oscillators added need to be of same dimension');
            end
            
            if isequal(osc1.phaseType,osc2.phaseType)
                oscOut.phaseType=osc1.phaseType;
                if isequal(oscOut.phaseType,'linear')
                    oscOut.phaseTypeData=osc1.phaseTypeData+osc2.phaseTypeData;
                end
            end
            
            oscOut.phase=@(varargin) osc1.phase(varargin)+osc2.phase(varargin);
            
            
        end
    end
    
end

