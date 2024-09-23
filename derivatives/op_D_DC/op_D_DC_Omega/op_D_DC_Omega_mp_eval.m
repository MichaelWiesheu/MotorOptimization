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


function A = op_D_DC_Omega_mp_eval(spg, spgEval, msh, mshEval, patch_list)

    rows = cell(msh.npatch, 1);
    dims = cell(msh.npatch, 1);
    vals = cell(msh.npatch, 1);
    
    for iptc = patch_list
    
        [rs, dm, vs] = op_D_DC_Omega(spgEval{iptc}, mshEval{iptc});
    
        rows{iptc} = spg.gnum{iptc}(rs);
        dims{iptc} = dm;
    
        if (~isempty (spg.dofs_ornt))
            warning("op_D_DC_Omega_mp: dofs_ornt not tested yet")
            vs = vs .* spg.dofs_ornt{iptc}(vs)';
        end
    
        vals{iptc} = vs;
    end
    
    A = sparse(cell2mat(rows), cell2mat(dims), cell2mat(vals), spg.ndof, msh.ndim);

end
