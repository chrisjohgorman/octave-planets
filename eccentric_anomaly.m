function [e1] = eccentric_anomaly(M, e, tol)
    e0 = M + (180/pi) * e * sind(M) * (1 + e + cosd(M));
        e1 = e0 - (e0 - (180/pi) * e * sind(e0) - M) / (1 - e * cosd(e0));
        while abs(e0-e1) > tol
        e0 = e1;
        e1 = e0 - (e0 - (180/pi) * e * sind(e0) -M) / (1 -e * cosd(e0));
    end
end
