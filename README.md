# Code snapshot to the publication "Combined Parameter and Shape Optimization of Electric Machines with Isogeometric Analysis"
- Published with Springer on 21.09.2024 
- Online link: https://link.springer.com/article/10.1007/s11081-024-09925-0

This repository is created for the purpose of reproducing of the results shown in the paper "Combined Parameter and Shape Optimization of Electric Machines with Isogeometric Analysis".
The aim is to promote the FAIR (findable, accessible, interoperable, reusable) standards and to help researchers to understand the methodology. Note, that the code is made for trying and prototyping new algorithms, usage is at ones own risk! 

# Prerequisites

In order to run the code you need
- Matlab (tested with version 2024a)
- Matlab Optimization Toolbox
- Matlab Curve Fitting Toolbox
- Matlab Global Optimization Toolbox

No other packages need to be downloaded at this point. To provide a fully working snapshot, we already include here the code snapshots for the packages

- GeoPDEs (https://rafavzqz.github.io/geopdes/)
- The NURBS toolbox (https://de.mathworks.com/matlabcentral/fileexchange/26390-nurbs-toolbox-by-d-m-spink)

to which we give credit.


# How to run the code

In the results folder, there are several scrips that can be run. These scripts reproduce the results from the publication. For the optimization, run the scripts:

- opt_param.m - runs the parameter optimization
- opt_shape.m - runs the shape optimization
- opt_param_then_shape.m - runs the sequential parameter and shape optimization
- opt_param_shape.m - runs the combined parameter and shape optimization

NOTE: We have observed that the number of iterations for the optimizations can sometimes be slightly different depending on the hardware. 
This is due to the high number of numerical solves, which can introduce deviations on a numerical level, that are passed to the next iteration. The overall convergence is not affected. 
We included also the results for our optimization runs (see the "OptimizationHistory.mat" files) for comparison.

For the evaluation of the results, run the scripts:

- eval_init.m - evaluation of the initial geometry without optimization
- eval_param.m - evaluation of the parameter optimization
- eval_shape.m - evaluation of the shape optimization
- eval_param_then_shape.m - evaluation of the sequential parameter and shape optimization
- eval_param_shape.m - evaluation of the combined parameter and shape optimization
- eval_opt_torques. - evaluation of the torque profiles generated from the optimization
- eval_JMAG_GA.m - evaluation of the results from the JMAG genetic algorithm

# License

These files are free to use for research purposes. You can edit it or modify it as your needs under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

These files are distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

While the code is under GPL 3.0, the geometrical motor design is taken from the JMAG-RT Motor Model Library and is subject to JSOL's terms of use https://www.jmag-international.com/modellibrary. Material data is taken from FEMM which is published under the Aladdin Free Public License https://www.femm.info/wiki/license.

# Authors: 
Michael Wiesheu, Theodor Komann, Melina Merkel, Sebastian Schöps, Stefan Ulbrich, Idoia Cortes Garcia

# Support:
michael.wiesheu@tu-darmstadt.de
komann@mathematik.tu-darmstadt.de
