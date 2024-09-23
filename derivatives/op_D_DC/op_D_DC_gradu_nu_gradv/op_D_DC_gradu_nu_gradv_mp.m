% OP_D_DC_GRADU_NU_GRADV_MP: assemble the derivative of the stiffness matrix A_dP = [a(i,j)], a(i,j) = (nu grad u_j, grad v_i), in a multipatch domain with respect to the geometry control points of g_k.
%
%   mat = op_D_DC_gradu_nu_gradv_mp (spu, spv, spg, msh, u, material, patches);
%
% INPUT:
%
%   spu:     object representing the space of trial functions (see sp_multipatch)
%   spv:     object representing the space of test functions (see sp_multipatch)
%   spg:     object representing the geometry space
%   msh:     object defining the domain partition and the quadrature rule (see msh_multipatch)
%   u:       discrete magnetic vector potential
%   mat:     material class with reluctivity functions
%   patches: list of patches where the integrals have to be computed. By default, all patches are selected.
%
% OUTPUT:
%
%   mat:    derivative of assembled stiffness matrix by control points

function A = op_D_DC_gradu_nu_gradv_mp(spu, spv, spg, msh, u, material, patch_list)
    
    if ((spu.npatch~= spv.npatch)  || (spv.npatch~= spg.npatch)  || (spu.npatch ~= msh.npatch))
        error ('op_D_DC_gradu_gradv_mp: the number of patches does not coincide')
    end
    
    rws = [];
    cls = [];
    tns = [];
    dms = [];
    vals = [];
    ncounter = 0;
    for iptc = patch_list
    
        [indices, vs] = op_D_DC_gradu_nu_gradv_tp (spu.sp_patch{iptc}, spv.sp_patch{iptc}, spg.sp_patch{iptc}, msh.msh_patch{iptc}, u(spu.gnum{iptc}), material);

        if isempty(indices)
            continue
        end
    
        rs = indices(:,1);
        cs = indices(:,2);
        ts = indices(:,3);
        dm = indices(:,4);
    
        rws(ncounter+(1:numel (rs))) = spu.gnum{iptc}(rs);
        cls(ncounter+(1:numel (cs))) = spv.gnum{iptc}(cs);
        tns(ncounter+(1:numel (ts))) = spg.gnum{iptc}(ts);
        dms(ncounter+(1:numel (dm))) = dm;
    
        if (~isempty (spu.dofs_ornt))
            vs = spu.dofs_ornt{iptc}(rs)' .* vs;
        end
        if (~isempty (spv.dofs_ornt))
            vs = vs .* spv.dofs_ornt{iptc}(cs)';
        end
        if (~isempty (spg.dofs_ornt))
            warning("op_D_DC_gradu_gradv_mp: dofs_ornt not tested yet")
            vs = vs .* spg.dofs_ornt{iptc}(ts)';
        end
    
        vals(ncounter+(1:numel (rs))) = vs;
        ncounter = ncounter + numel (rs);
    end

    rws = reshape(rws, [], 1);
    cls = reshape(cls, [], 1);
    tns = reshape(tns, [], 1);
    dms = reshape(dms, [], 1);
    vals = reshape(vals, [], 1);
    indices = [rws, cls, tns, dms];
    
    A = sptensor (indices, vals, [spv.ndof, spu.ndof, spg.ndof, msh.ndim]);

end
