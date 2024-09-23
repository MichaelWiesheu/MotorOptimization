% OP_OMEGA: calculate the volume of the mesh
%
%   vol = op_Omega (msh);
%
% INPUT:
%    msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%
% OUTPUT:
%
%   vol:    volume/(area)

function vol = op_Omega(msh)

    vol = sum(msh.jacdet.*msh.quad_weights, 'all');

end

