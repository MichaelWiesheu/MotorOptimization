function nrb2iges_bound(nrb, filename)
    if nargin < 2
        filename = 'nrb.igs';
    end

    [geo, bnd, intrfc] = mp_geo_load (nrb);

    cnt = 1;

    for ibnd = 1:numel(bnd)
        nrb_export(cnt) = nrbextract(geo(bnd(ibnd).patches).nurbs, bnd(ibnd).faces);
        cnt = cnt + 1;
    end
    for iref = 1:numel(intrfc)
        nrb_export(cnt) = nrbextract(geo(intrfc(iref).patch1).nurbs, intrfc(iref).side1);
        cnt = cnt + 1;
    end
    %     for i = 1:numel(geo)
    %         for j = 1:4
    %             nrb_export(cnt) = nrbextract(geo(i).nurbs, j);
    %             cnt = cnt +1;
    %         end
    %     end

    scaling = 1000;
    for inrb = 1:numel(nrb_export)
        nrb_export(inrb).coefs(1:2,:) = nrb_export(inrb).coefs(1:2,:)*scaling;
    end

    nrb2iges(nrb_export, filename)
end