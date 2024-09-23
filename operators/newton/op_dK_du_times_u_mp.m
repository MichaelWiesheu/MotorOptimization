% OP_DK_DU_TIMES_U_MP: assemble the derivative of the magnetic stiffness matrix A =
% [a(i,j)], a(i,j) = (nu grad u_j, grad v_i) w.r.t magnetic vector
% potential multiplier for the newton scheme in a multipatch domain.
% Directly multiplied by u for faster evaluation
%
%   mat = op_dK_du_times_u_mp(sp_mp1, sp_mp2, msh, u, material, mu_curve, [patches]);
%
% INPUT:
%
%   sp_mp1:     object representing the space of trial functions (see sp_multipatch)
%   sp_mp2:     object representing the space of test functions (see sp_multipatch)
%   msh:        object defining the domain partition and the quadrature rule (see msh_multipatch)
%   u:          solution vector for the magnetic vector potential of the domain
%   material:   struct containing a B/H curve TBD!
%   mu_curve:   "Iron" so far for TBD!
%   patches: list of patches where the integrals have to be computed. By default, all patches are selected.
%
% OUTPUT:
%
%   mat:    derivative of assembled stiffness matrix magnetic vector potential


function A = op_dK_du_times_u_mp(sp_mp1, sp_mp2, msh, u, material, patch_list)

    if (nargin < 6)
        patch_list = 1:msh.npatch;
    end
    
    if ((sp_mp1.npatch~= sp_mp2.npatch)  || (sp_mp1.npatch ~= msh.npatch))
        error ('op_dK_du_times_u_mp: the number of patches does not coincide')
    end
    
    ncounter = 0;
    rws = [];
    cls = [];
    vals = [];
    for iptc = patch_list
        [rs, cs, vs] = op_dK_du_times_u_tp (sp_mp1.sp_patch{iptc}, sp_mp2.sp_patch{iptc}, msh.msh_patch{iptc}, u(sp_mp1.gnum{iptc}), material);
        if isempty(rs)
            continue
        end
        
        rws(ncounter+(1:numel (rs))) = sp_mp1.gnum{iptc}(rs);
        cls(ncounter+(1:numel (cs))) = sp_mp2.gnum{iptc}(cs);
    
        if (~isempty (sp_mp1.dofs_ornt))
            vs = sp_mp1.dofs_ornt{iptc}(rs)' .* vs;
        end
        if (~isempty (sp_mp2.dofs_ornt))
            vs = vs .* sp_mp2.dofs_ornt{iptc}(cs)';
        end
    
        vals(ncounter+(1:numel (rs))) = vs;
        ncounter = ncounter + numel (rs);
    end
    rws = reshape(rws, [], 1);
    cls = reshape(cls, [], 1);
    vals = reshape(vals, [], 1);
    
    A = sparse(rws, cls, vals, sp_mp1.ndof, sp_mp2.ndof);
    
    clear rws cls vals rs cs vs

end