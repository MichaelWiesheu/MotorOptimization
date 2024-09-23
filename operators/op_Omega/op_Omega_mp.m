% OP_OMEGA_MP: calculate the volume of the mesh
%
%   vol = op_Omega_mp (msh);
%
% INPUT:
%    msh:   structure containing the domain partition and the quadrature rule
%
% OUTPUT:
%
%   vol:    volume/(area)

function vol = op_Omega_mp(msh, patch_list)
    
    vol = 0;
    
    for iptc = patch_list
      vol = vol + op_Omega(msh_precompute(msh.msh_patch{iptc}));
    end
end
