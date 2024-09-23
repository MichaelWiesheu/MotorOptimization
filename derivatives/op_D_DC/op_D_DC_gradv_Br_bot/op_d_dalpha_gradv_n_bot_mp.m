
function rhs = op_d_dalpha_gradv_n_bot_mp (space, msh, alpha, patch_list)

  if (nargin < 4)
    patch_list = 1:msh.npatch;
  end

  if (space.npatch ~= msh.npatch)
    error ('op_gradv_n_bot_mp: the number of patches does not coincide')
  end
  
  rhs = zeros (space.ndof, 1);
  for iptc = patch_list
    rhs_loc = op_d_dalpha_gradv_n_bot_tp (space.sp_patch{iptc}, msh.msh_patch{iptc}, alpha);
    
    if (~isempty (space.dofs_ornt))
      rhs_loc = space.dofs_ornt{iptc}(:) .* rhs_loc(:);
    end
    rhs(space.gnum{iptc}) = rhs(space.gnum{iptc}) + rhs_loc;
  end

end