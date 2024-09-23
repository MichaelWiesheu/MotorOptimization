function dCdP = getdCdP(geometryFunction, params, spaceGeo, Steps)
    nParameters = numel(fieldnames(params));
    if ~exist("Steps", "var")
        Steps = ones(nParameters, 1)*1e-6;
    end

    C0 = getC(geometryFunction, params, spaceGeo);
    dCdP = zeros([size(C0), nParameters]);
    
    for iParam = 1:nParameters
        Step = Steps(iParam);
        parameterNames = fieldnames(params);
        opts1 = params;
        opts1.(parameterNames{iParam}) = opts1.(parameterNames{iParam})+Step;
        opts2 = params;
        opts2.(parameterNames{iParam}) = opts2.(parameterNames{iParam})-Step;
        
        opts1.draw_geometry = false;
        opts2.draw_geometry = false;
        C1 = getC(geometryFunction, opts1, spaceGeo);
        C2 = getC(geometryFunction, opts2, spaceGeo);
    
        dCdP(:, :, iParam) = (C1 - C2)/(2*Step);
        % dCdP(:, :, iParam) = (C1 - C0)/(Step);
    end
end
