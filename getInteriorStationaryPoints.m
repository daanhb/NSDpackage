function [ interiorStationaryPoints ] = getInteriorStationaryPoints( stationaryPoints, a, b )
    %determine which stationary points are in the interior, as these
    %must be handled differently to the special case of endpoint
    interiorStationaryPoints=[];
    for j=1:length(stationaryPoints)
        if ~ismember(stationaryPoints(j),[a b])
            interiorStationaryPoints=[interiorStationaryPoints stationaryPoints(j)];
        end
    end


end

