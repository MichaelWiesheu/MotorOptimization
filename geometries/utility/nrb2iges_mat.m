function nrb2iges_mat(nrb, materials, filename)
    if nargin < 3
        filename = 'nrb.igs';
    end

    [geo, bnd, intrfc] = mp_geo_load (nrb);

    cnt = 1;

    for ibnd = 1:numel(bnd)
        nrb_export(cnt) = nrbextract(geo(bnd(ibnd).patches).nurbs, bnd(ibnd).faces);
        cnt = cnt + 1;
    end

    mats = fieldnames(materials);

    for iref = 1:numel(intrfc)
        patch1 = intrfc(iref).patch1;
        patch2 = intrfc(iref).patch2;
        matsEqual = false;
        for iMat = 1:numel(mats)
            if any(patch1 == materials.(mats{iMat})) && any(patch2 == materials.(mats{iMat}))
                matsEqual = true;
            end
        end
        if matsEqual == true
            continue
        end
        nrb_export(cnt) = nrbextract(geo(intrfc(iref).patch1).nurbs, intrfc(iref).side1);
        cnt = cnt + 1;
    end

    scaling = 1000;
    for inrb = 1:numel(nrb_export)
        nrb_export(inrb).coefs(1:2,:) = nrb_export(inrb).coefs(1:2,:)*scaling;
    end

    nrb2iges(nrb_export, filename)
end

