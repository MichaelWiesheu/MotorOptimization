% OP_D_DC_OMEGA_MP: assemble the derivative of the patch volume w.r.t the control points of g.
%
%   mat = op_D_DC_Omega_mp (spg, msh, patches);
%
% INPUT:
%
%   spg:     object representing the geometry space
%   msh:     object defining the domain partition and the quadrature rule (see msh_multipatch)
%   patches: list of patches where the integrals have to be computed. By default, all patches are selected.
%
% OUTPUT:
%
%   mat:    derivative of assembled stiffness matrix by control points


function A = op_D_DC_Omega_mp(spg, msh, patch_list)
    
    if (spg.npatch~= msh.npatch)
        error ('op_D_DC_Omega_mp: the number of patches does not coincide')
    end
    
    rows = [];
    dms = [];
    vals = [];
    ncounter = 0;
    for iptc = patch_list
    
      [rs, dm, vs] = op_D_DC_Omega_tp (spg.sp_patch{iptc}, msh.msh_patch{iptc});

        if isempty(rs)
            continue
        end
    
        rows(ncounter+(1:numel (rs))) = spg.gnum{iptc}(rs);
        dms(ncounter+(1:numel (dm))) = dm;
    
        if (~isempty (spg.dofs_ornt))
            warning("op_D_DC_Omega_mp: dofs_ornt not tested yet")
            vs = vs .* spg.dofs_ornt{iptc}(vs)';
        end
    
        vals(ncounter+(1:numel (rs))) = vs;
        ncounter = ncounter + numel (rs);
    end

    rows = reshape(rows, [], 1);
    dms = reshape(dms, [], 1);
    vals = reshape(vals, [], 1);

    A = sparse(rows, dms, vals, spg.ndof, msh.ndim);

end
