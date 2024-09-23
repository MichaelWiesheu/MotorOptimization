function ControlPoints = getC(geometryFunction, params, spaceGeo)
    params.draw_geometry = false;
    srf = geometryFunction(params);

    for iPatch  = 1:spaceGeo.npatch
        ind_loc = spaceGeo.gnum{iPatch};
        ControlPoints(ind_loc, 1) = reshape(srf(iPatch).coefs(1, :, :, :)./srf(iPatch).coefs(4, :, :, :), [], 1);
        ControlPoints(ind_loc, 2) = reshape(srf(iPatch).coefs(2, :, :, :)./srf(iPatch).coefs(4, :, :, :), [], 1);
        ControlPoints(ind_loc, 3) = reshape(srf(iPatch).coefs(3, :, :, :)./srf(iPatch).coefs(4, :, :, :), [], 1);
        ControlPoints(ind_loc, 4) = reshape(srf(iPatch).coefs(4, :, :, :), [], 1);
    end
end
