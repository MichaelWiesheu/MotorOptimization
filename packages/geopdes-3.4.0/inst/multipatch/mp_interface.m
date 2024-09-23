% MP_INTERFACE: create a global numbering of basis functions in multipatch geometries.
%   The interfaces must be conforming: same control points and same knot vector
%
%   [glob_num, glob_ndof] = mp_interface (interfaces, sp);
%
% INPUT:
%
%  interfaces: structure with the information of the interfaces between patches (see mp_geo_load)
%  sp:         object representing the space of discrete functions, with the grad-preserving transform (see sp_scalar)
%
% OUTPUT:
%
%  glob_num:  global numbering of the discrete basis functions
%  glob_ndof: total number of degrees of freedom
%
% Copyright (C) 2011, 2015 Rafael Vazquez
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

function [glob_num, glob_ndof] = mp_interface (interfaces, sp)

  npatch = numel (sp);

  ndim = numel (sp{1}.ndof_dir);

  if (~isempty (interfaces))
    glob_num = cell (npatch, 1);
    patch_intrfc = cell (npatch, 1);
    ttform   = cell (npatch, numel (interfaces));
    ppnum    = cell (numel (interfaces), 1);
    for iptc = 1:npatch
      glob_num{iptc} = zeros (sp{iptc}.ndof, 1);
      patch_intrfc{iptc} = union (find([interfaces.patch1] == iptc), ...
                                  find([interfaces.patch2] == iptc));
    end

    for intrfc = 1:numel(interfaces)
      iptc1  = interfaces(intrfc).patch1;
      iptc2  = interfaces(intrfc).patch2;
      iside1 = interfaces(intrfc).side1;
      iside2 = interfaces(intrfc).side2;
      ttform{iptc1, intrfc} = sp{iptc1}.boundary(iside1).dofs;
      
      if (ndim == 1)
        ttform{iptc2, intrfc} = sp{iptc2}.boundary(iside2).dofs;
      elseif (ndim == 2)
        if (interfaces(intrfc).ornt == 1)
          ttform{iptc2, intrfc} = sp{iptc2}.boundary(iside2).dofs;
        else
          ttform{iptc2, intrfc} = fliplr (sp{iptc2}.boundary(iside2).dofs);
        end
      elseif (ndim == 3)
        nghbr_dofs = reshape (sp{iptc2}.boundary(iside2).dofs, ...
                      sp{iptc2}.boundary(iside2).ndof_dir);
        if (interfaces(intrfc).flag == -1)
          nghbr_dofs = nghbr_dofs';
        end
        if (interfaces(intrfc).ornt1 == -1)
          nghbr_dofs = flipud (nghbr_dofs);
        end
        if (interfaces(intrfc).ornt2 == -1)
          nghbr_dofs = fliplr (nghbr_dofs);
        end
        ttform{iptc2, intrfc} = nghbr_dofs(:)';
      end
      ppnum{intrfc} = zeros (1, sp{iptc1}.boundary(iside1).ndof);
    end

    glob_ndof = 0;
% We start with the dofs that do not belong to any interface
    for iptc = 1:npatch
      non_intrfc_dofs = setdiff(1:sp{iptc}.ndof, [ttform{iptc,:}]);
      glob_num{iptc}(non_intrfc_dofs) = glob_ndof + (1:numel(non_intrfc_dofs));
      glob_ndof = glob_ndof + numel (non_intrfc_dofs);
    end

% Then we set the interfaces
    for intrfc = 1:numel(interfaces)
      iptc     = interfaces(intrfc).patch1;
      new_dofs = find (glob_num{iptc}(ttform{iptc, intrfc}) == 0);
      new_dofs_number = glob_ndof + (1:numel (new_dofs));
      glob_num{iptc}(ttform{iptc, intrfc}(new_dofs)) = new_dofs_number;
      ppnum{intrfc}(new_dofs) = new_dofs_number;

      [glob_num, ppnum] = set_same_interface (iptc, intrfc, ttform, new_dofs,...
				   interfaces, glob_num, ppnum, patch_intrfc);
      [glob_num, ppnum] = set_same_patch (iptc, intrfc, ttform, new_dofs, ...
                                   interfaces, glob_num, ppnum, patch_intrfc);
      glob_ndof = glob_ndof + numel (new_dofs);

    end

  else
    glob_ndof = 0;
    glob_num = cell (npatch, 1);
    for iptc = 1:npatch
      glob_num{iptc} = glob_ndof + (1:sp{iptc}.ndof);
      glob_ndof = glob_ndof + sp{iptc}.ndof;
    end
  end    

end

function [glob_num, ppnum] = set_same_patch (iptc, intrfc, ttform, new_dofs, ...
                                 interfaces, glob_num, ppnum, patch_intrfc);

  intrfc_dofs = ttform{iptc, intrfc}(new_dofs);
  for ii = setdiff (patch_intrfc{iptc}, intrfc)
    [common_dofs, pos1, pos2] = intersect (ttform{iptc, ii}, intrfc_dofs);
    not_set = find (ppnum{ii}(pos1) == 0);
    if (~isempty (not_set))
      ppnum{ii}(pos1(not_set)) = ppnum{intrfc}(new_dofs(pos2(not_set)));
      [glob_num, ppnum] = set_same_interface (iptc, ii, ttform, ...
                   pos1(not_set), interfaces, glob_num, ppnum, patch_intrfc);
    end
  end
end

function [glob_num, ppnum] = set_same_interface (iptc, intrfc, ttform, ...
                         new_dofs, interfaces, glob_num, ppnum, patch_intrfc)
  intrfc_dofs = ttform{iptc, intrfc}(new_dofs);
  iptc2 = setdiff ([interfaces(intrfc).patch1 interfaces(intrfc).patch2], iptc);
  glob_num{iptc2}(ttform{iptc2, intrfc}(new_dofs)) = ...
                            glob_num{iptc}(intrfc_dofs);
  [glob_num, ppnum] = set_same_patch (iptc2, intrfc, ttform, new_dofs, ...
			     interfaces, glob_num, ppnum, patch_intrfc);

end
