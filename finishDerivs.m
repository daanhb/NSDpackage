function G = finishDerivs( G, newLength, Npts, zeroThresh )
%newLength=(max(pathPowers)+1)
%zeroThresh=rectTol
    origGlength=length(G);
    for extraG=(origGlength+1) : newLength
        G{extraG}=@(z) arrayfun(@(Z) CauchyDiff( Z,G{origGlength}, Z, zeroThresh, extraG-origGlength, Npts ), z);
    end

end