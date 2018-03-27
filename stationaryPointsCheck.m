function stationaryPointsCheck(G,gStationaryPoints)
    dgdx = G{2};
    thresh = 1E-10;
    for s = gStationaryPoints
        if abs(dgdx(s))>thresh
            error(sprintf('%f+%fi is not a stationary point of phase provided',real(s),imag(s)));
        end
    end
end

