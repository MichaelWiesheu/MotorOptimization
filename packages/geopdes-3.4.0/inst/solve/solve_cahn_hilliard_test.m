% SOLVE_CAHN_HILLIARD: solve the Cahn-Hilliard equation, with a generalized alpha discretization in time.
%
% The functions solves the problem of finding u such that
%
%  du/dt - Delta (f(u) - lambda*Delta u) = 0
%
% with Delta the Laplacian, and f(u) = alpha u^3 - beta u, and periodic boundary conditions.
%
% The values of alpha and beta (or f itself) can be changed in op_gradmu_gradv_tp.
%
% For details on the problem and the formulation, see
%  H. Gomez, V.M. Calo, Y. Bazilevs, T.J.R. Hughes, CMAME 197 (2008), 4333-4352.
%  H. Gomez, A. Reali, G. Sangalli, J. Comput. Physics 262 (2014), 153-171.
%
% USAGE:
%
%   [geometry, msh, space, results] = solve_cahn_hilliard (problem_data, method_data, save_info)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - periodic_directions: parametric directions along which to apply periodic conditions (may be empty)
%    - lambda:       parameter representing the length scale of the problem, and the width of the interface
%    - Time_max:     final time
%    - fun_u:        initial condition. Equal to zero by default.
%    - fun_udot:     initial condition for time derivative. Equal to zero by default.
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree of the spline functions.
%    - regularity: continuity of the spline functions.
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%    - dt:         time step size for generalized-alpha method
%    - rho_inf_gen_alpha: parameter in [0,1], which governs numerical damping of the generalized alpha method
%
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      mesh object that defines the quadrature rule (see msh_cartesian)
%  space:    space object that defines the discrete space (see sp_scalar)
%  results:  a struct with the saved results, containing the following fields:
%    - time: (array of length Ntime) time at which the solution was saved
%    - u:    (size ndof x Ntime) degrees of freedom for the solution
%    - udot: (size ndof x Ntime) degrees of freedom for the time derivative
%
% Only periodic and Neumann boundary conditions are implemented. Neumann
%  conditions are considered by default.
%
% Copyright (C) 2023 Michele Torre, Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [geometry, msh, space, results] = solve_cahn_hilliard_test (problem_data, method_data, save_info)

%%-------------------------------------------------------------------------
% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

%%-------------------------------------------------------------------------
% Construct geometry structure
geometry  = geo_load (geo_name);
[knots, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);

%%-------------------------------------------------------------------------
% Check for periodic conditions, and consistency with other boundary conditions
if (exist('periodic_directions', 'var'))
  knots = kntunclamp (knots, degree, regularity, periodic_directions);
else
  periodic_directions = [];
end

if (exist ('nmnn_sides','var') && ~isempty (nmnn_sides))
  disp('User defined Neumann sides deleted')
  clear nmnn_sides
end

%%-------------------------------------------------------------------------
% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry);
  
% Construct space structure
space    = sp_bspline (knots, degree, msh, [], periodic_directions);

%%-------------------------------------------------------------------------
% Generalized alpha parameters
a_m = .5*((3-rho_inf_gen_alpha)/(1+rho_inf_gen_alpha));
a_f =  1/(1+rho_inf_gen_alpha);
gamma =  .5 + a_m - a_f;

%%-------------------------------------------------------------------------
% No flux b.c. (essential boundary condition)
% Set Neumann boundary conditions for non-periodic sides
nmnn_sides   = [];
for idir = 1:msh.ndim
  if (~ismember(idir, periodic_directions))
    nmnn_sides = [nmnn_sides, 2*(idir-1)+[1 2]];
  end
end

% Matrix and data for Neumann
[BC, indBC, dof_free] = impose_essential_neumann_weak(nmnn_sides, space, msh);

%%-------------------------------------------------------------------------
% Assemble the mass matrix, the Laplacian matrix, and the boundary term
mass_mat = op_u_v_tp (space, space, msh);
lapl_mat = op_laplaceu_laplacev_tp (space, space, msh, lambda);
term4 = int_boundary_term (space, msh,  lambda, nmnn_sides);

%%-------------------------------------------------------------------------
% Initial conditions
time = 0;
if (exist('fun_u', 'var') && ~isempty(fun_u))
  rhs = op_f_v_tp (space, msh, fun_u);
  u_n = mass_mat\rhs;
else
  u_n = zeros(space.ndof, 1);
end

if (exist('fun_udot', 'var') && ~isempty(fun_udot))
  rhs = op_f_v_tp(space, msh, fun_udot);
  udot_n = mass_mat\rhs;
else
  udot_n = zeros(space.ndof, 1);
end

norm_flux_initial = check_flux_phase_field(space, msh, u_n);
disp(strcat('initial flux =',num2str(norm_flux_initial)))

%%-------------------------------------------------------------------------
% Initialize structure to store the results
save_id = 1;
results.u = zeros(length(u_n), length(save_info)+1);
results.udot = zeros(length(u_n), length(save_info)+1);
results.time = zeros(length(save_info)+1,1);
flag_stop_save = false;

% Save initial conditions
results.u(:,1) = u_n;
results.udot(:,1) = udot_n;
results.time(1) = time;

%%-------------------------------------------------------------------------
% Loop over time steps
while time < Time_max
  disp('----------------------------------------------------------------')
  disp(strcat('time step t=',num2str(time)))

  [u_n1, udot_n1] = generalized_alpha_step(u_n, udot_n, dt, a_m, a_f, gamma, mass_mat, lapl_mat, term4, BC, indBC, dof_free, space, msh);

  % Check flux through the boundary
  norm_flux = check_flux_phase_field(space, msh, u_n1);
  disp(strcat('flux variation =',num2str(norm_flux - norm_flux_initial)))

  % Time step update
  time = time + dt;
  u_n = u_n1;
  udot_n = udot_n1;

  % Check max time
  if (time + dt > Time_max)
    dt = Time_max-time;
  end

  % Store results
  if (flag_stop_save == false)
    if (time >= save_info(save_id))
      save_id = save_id + 1;
      results.u(:,save_id) = u_n;
      results.udot(:,save_id) = udot_n;
      results.time(save_id) = time;
      if (save_id > length(save_info))
        flag_stop_save = true;
      end
    end
  end
end

disp('----------------------------------------------------------------')
disp(strcat('END ANALYSIS t=',num2str(time)))
disp('----------------------------------------------------------------')

% Crop results
results.u = results.u(:,1:save_id);
results.udot = results.udot(:,1:save_id);
results.time = results.time(1:save_id);

end

%--------------------------------------------------------------------------
% FUNCTIONS
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Neumann boundary conditions
%--------------------------------------------------------------------------
function [BC, indBC, dof_free] = impose_essential_neumann_weak(nmnn_sides, space, msh)

if (~isempty(nmnn_sides))

    BC =  spalloc (space.ndof, space.ndof, 3*space.ndof);
    indBC = [];

    for iside = nmnn_sides
        indBC = union (indBC, space.boundary(iside).dofs);
        msh_side = msh_eval_boundary_side (msh, iside);
        msh_side_int = msh_boundary_side_from_interior (msh, iside);
        sp_side = space.constructor (msh_side_int);
        sp_side = sp_precompute (sp_side, msh_side_int, 'value', true, 'gradient', true);
        coe_side = ones(msh_side.nqn, msh_side.nel);
        BC =  BC + op_gradv_n_u(sp_side ,sp_side ,msh_side, coe_side);
    end
    BC = BC';
    dof_free = setdiff(1:space.ndof, indBC);

else
    BC = [];
    indBC = [];
    dof_free = 1:space.ndof;
end


end

%--------------------------------------------------------------------------
% One step of generalized alpha-method
%--------------------------------------------------------------------------
function [u_n1, udot_n1] = generalized_alpha_step(u_n, udot_n, dt, a_m, a_f, gamma, mass_mat, lapl_mat, term4, BC, ind_BC, dof_free, space, msh)

% Convergence criteria
  n_max_iter = 20;
  tol_rel_res = 1e-10;

% Predictor step
  u_n1 = u_n;
  udot_n1 = (gamma-1)/gamma * udot_n; 

% Newton loop
  for iter = 0:n_max_iter
  % Field at alpha level
    udot_a = udot_n + a_m *(udot_n1-udot_n);
    u_a = u_n + a_f *(u_n1-u_n);

  % Compute the residual (internal)
    [Res_gl, stiff_mat] = Res_K_cahn_hilliard (space, msh, mass_mat, lapl_mat, term4, u_a, udot_a);

  % Convergence check
    if iter == 0
      norm_res_0 = norm(Res_gl(dof_free));
    end
    norm_res = norm(Res_gl(dof_free));
    

    if norm_res/norm_res_0 < tol_rel_res
      disp(strcat('iteration n° = ',num2str(iter)))
      disp(strcat('norm of absolute residual = ',num2str(norm_res)))
      break
    end
    if iter == n_max_iter
      disp(strcat('Newton reached the maximum number of iterations'))
      disp(strcat('norm of absolute residual = ',num2str(norm_res)))
    end
    
  % Compute the update, and update the solution
    A_gl = a_m * mass_mat + a_f * gamma *dt * stiff_mat;
    if (~isempty(ind_BC))
      A_gl(ind_BC,:)=0;
      A_gl = A_gl +BC;
      Res_gl(ind_BC) = 0;
    end
    d_udot = - A_gl\Res_gl;

    udot_n1 = udot_n1 + d_udot;
    u_n1 = u_n1 + gamma * dt* d_udot;
  end

end

%--------------------------------------------------------------------------
% Canh-Hilliard residual and tangent matrix
%--------------------------------------------------------------------------
function [Res_gl, stiff_mat] = Res_K_cahn_hilliard(space, msh, mass_mat, lapl_mat, term4, u_a, udot_a)

    % Double well (matrices)
    [term2, term2K] = op_gradmu_gradv_tp(space, msh, u_a);    
 
    % Residual
    Res_gl = mass_mat*udot_a + term2*u_a  + lapl_mat*u_a;

    % Tangent stiffness matrix (mass is not considered here)
    stiff_mat = term2 + term2K + lapl_mat;
    
    if (size(term4,1) > 0)
      Res_gl = Res_gl - term4 * u_a;
      stiff_mat = stiff_mat - term4;
    end

end


%--------------------------------------------------------------------------
% Integral of the double-well function
%--------------------------------------------------------------------------
function [A, B] = op_gradmu_gradv_tp (space, msh,  uhat)

% Coefficients of the double well function.
  alpha = 1;
  beta = 1;

  for idim = 1:msh.ndim
    size1 = size (space.sp_univ(idim).connectivity);
    if (size1(2) ~= msh.nel_dir(idim))
      error ('The discrete space is not associated to the mesh')
    end
  end

  A = spalloc (space.ndof, space.ndof, 3*space.ndof);
  B = spalloc (space.ndof, space.ndof, 6*space.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col, 'gradient', true);

    % Evaluate the field and its gradient at the Gaussian points
    utemp = sp_eval_msh (uhat, sp_col, msh_col, {'value', 'gradient'});
    u = utemp{1};
    gradu = utemp{2};

    % Polynomial formulation for the double-well
    coeffs_A = 3.* alpha .* u.^2 - beta;
    A = A + op_gradu_gradv (sp_col, sp_col, msh_col, coeffs_A);

    coeffs_B = 6.* alpha .* u;
    coeffs_Bv = gradu;
    for idim = 1:msh.ndim
      coeffs_Bv(idim,:,:) = coeffs_Bv(idim,:,:) .* reshape(coeffs_B, 1, size(coeffs_B,1), size(coeffs_B,2));
    end
    B = B + op_vel_dot_gradu_v (sp_col, sp_col, msh_col, coeffs_Bv)';
  end
end

%--------------------------------------------------------------------------
% Integral of the boundary term
%--------------------------------------------------------------------------
function [A] = int_boundary_term (space, msh,  lambda, nmnn_sides)

  if (~isempty (nmnn_sides))
    A =  spalloc (space.ndof, space.ndof, 3*space.ndof);
    for iside=1:length(nmnn_sides)
      msh_side = msh_eval_boundary_side (msh, nmnn_sides(iside));
      msh_side_int = msh_boundary_side_from_interior (msh, nmnn_sides(iside) );
      sp_side = space.constructor (msh_side_int);
      sp_side = sp_precompute (sp_side, msh_side_int, 'gradient', true, 'laplacian', true);

      for idim = 1:msh.rdim
        x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
      end
      coe_side = lambda (x{:});

      A =  A + op_gradv_n_laplaceu (sp_side, sp_side, msh_side, coe_side);
    end

    A = A';
  else
    A = [];
  end
  
end



%--------------------------------------------------------------------------
% Check flux through the boundaries
%--------------------------------------------------------------------------

function norm_flux = check_flux_phase_field(space, msh, uhat)

sides = [1,2,3,4];
norm_flux = 0;

for iside=1:length(sides)   

    msh_side = msh_eval_boundary_side (msh, sides(iside) ) ;
    msh_side_int = msh_boundary_side_from_interior (msh, sides(iside) ) ;
    sp_side = space.constructor ( msh_side_int) ;
    sp_side = sp_precompute ( sp_side , msh_side_int , 'gradient' , true );


    utemp = sp_eval_msh (uhat, sp_side, msh_side, {'gradient'});
    gradu = utemp{1};
    
    valu = zeros(1, size(msh_side.quad_weights,1), size(msh_side.quad_weights,2));
    for idim = 1:msh.rdim
        valu = valu + (gradu(idim,:,:) .* msh_side.normal(idim,:,:));
    end
    
    valu = reshape (valu, sp_side.ncomp, msh_side.nqn, msh_side.nel);
    w =msh_side.quad_weights .* msh_side.jacdet;
    errl2_elem = sum (reshape (sum ((valu).^2, 1), [msh_side.nqn, msh_side.nel]) .* w);
    errl2  = sqrt (sum (errl2_elem));


    norm_flux = norm_flux + errl2;
end

end