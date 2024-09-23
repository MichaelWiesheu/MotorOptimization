% OP_GRADU_NU_GRADV_MP: assemble the stiffness matrix A = [a(i,j)], a(i,j) = (nu grad u_j, grad v_i)
% where nu = nu(B) is given by b-h-curve for a multipatch domain
%
%   mat = op_gradu_nu_gradv_tp (spu, spv, msh, u, material, mu_curve, [patches]);
%
%
% INPUT:
%
%   spu:     object representing the space of trial functions (see sp_multipatch)
%   spv:     object representing the space of test functions (see sp_multipatch)
%   msh:     object defining the domain partition and the quadrature rule (see msh_multipatch)
%   u:     magnetic vector potential multiplier
%   material:   struct containing a B/H curve TBD!
%   mu_curve:   "Iron" so far for TBD!
%   patches: list of patches where the integrals have to be computed. By default, all patches are selected.

% OUTPUT:
%
%   mat:    assembled stiffness matrix

function A = op_gradu_nu_gradv_mp (spu, spv, msh, u, material, patch_list)

  if (nargin < 6)
    patch_list = 1:msh.npatch;
  end
  rows = [];
  cols = [];
  vals = [];
  
  ncounter = 0;
  for iptc = patch_list

      [rs, cs, vs] = op_gradu_nu_gradv_tp (spu.sp_patch{iptc}, spv.sp_patch{iptc}, msh.msh_patch{iptc}, u(spu.gnum{iptc}), material);

    rows(ncounter+(1:numel (rs))) = spv.gnum{iptc}(rs);
    cols(ncounter+(1:numel (rs))) = spu.gnum{iptc}(cs);

    if (~isempty (spv.dofs_ornt))
      vs = vs .* spv.dofs_ornt{iptc}(rs)';
    end
    if (~isempty (spu.dofs_ornt))
      vs = vs .* spu.dofs_ornt{iptc}(cs)';
    end
    
    vals(ncounter+(1:numel (rs))) = vs;
    ncounter = ncounter + numel (rs);
  end

  A = sparse (rows, cols, vals, spu.ndof, spu.ndof);
  clear rows cols vals rs cs vs
  
end