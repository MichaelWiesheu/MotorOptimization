function values = shp_eval_mp(sp_mp, geometry, eval_points, eval_patches, options)
    nopts = numel (options);
    for iopt = 1:nopts
        switch (lower (options{iopt}))
            case 'value'
                % shp x points
                values{iopt} = sparse(size(eval_points,2), sp_mp.ndof);
            case 'gradient'
                % 2x grad x points
                values{iopt} = sptensor([size(eval_points,2), 2, sp_mp.ndof]);
        end
    end
    
    for iPoint = 1:size(eval_points,2)
        p_i = eval_patches(iPoint);
        pnt = {eval_points(1,iPoint), eval_points(2,iPoint)};

        vals_i = shp_eval_sp(sp_mp.sp_patch{p_i}, geometry(p_i), pnt, options);

        for iopt = 1:nopts
            switch (lower (options{iopt}))
              
                case 'value'
                    values{iopt}(iPoint, sp_mp.gnum{p_i}) = vals_i{iopt};
                case 'gradient'
                    values{iopt}(iPoint, :, reshape(sp_mp.gnum{p_i}, 1, [])) = sptensor(squeeze(vals_i{iopt}));
            end
        end
    end
end