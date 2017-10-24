classdef singularity
    %information about singularity
    
    properties
        dim     %dimension of singularity, 0 for a single point, etc
        position   %can be anonymous function, or just a vector,
                    %depending on the dimension of the singularity
        blowUpType    %can be 'log', or 'a' where x^a is the blowup type
        pointType     %can be 'manifold' or 'point'
    end
    
    methods
        function self=singularity(dim, position, blowUpType, pointType)
            self.dim=dim;
            if strcmp(pointType,'manifold')
               if nargin(position)>dim
                  error('dimensions of manifold too high'); 
               end
            end
            self.position=position;
            self.blowUpType=blowUpType;
            self.pointType=pointType;
        end
    end
    
end

