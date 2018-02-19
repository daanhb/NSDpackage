function [criticalPoints,pathPowers] = makeCriticalPoints(a,b,gStationaryPoints,spPowers, RectTol,gAnalytic)
    %now determine if a or b are at a stationary point
    if isempty(gStationaryPoints)
        criticalPoints=[a b];
        pathPowers=[1 1];
    else
        criticalPoints=gStationaryPoints;
        pathPowers=spPowers;
        if min(abs(gStationaryPoints-a))>RectTol
            criticalPoints=[a criticalPoints];
            pathPowers=[1 pathPowers];
        end
        if min(abs(gStationaryPoints-b))>RectTol
            criticalPoints=[criticalPoints b];
            pathPowers=[pathPowers 1];
        end
    end
    if gAnalytic
        pathPowers=round(pathPowers);
    end
end