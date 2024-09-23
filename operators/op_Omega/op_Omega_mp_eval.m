% OP_OMEGA_M_EVAL: helper function that takes the preevaluated mesh and calculates the volume of the mesh
%
%   vol = op_Omega_mp (msh);
%
% INPUT:
%    msh:   structure containing the domain partition and the quadrature rule
%
% OUTPUT:
%
%   vol:    volume/(area)

function vol = op_Omega_mp_eval(mshEval, patch_list)
    
    vol = 0;
    
    for iptc = patch_list
      vol = vol + op_Omega(mshEval{iptc});
    end
end
