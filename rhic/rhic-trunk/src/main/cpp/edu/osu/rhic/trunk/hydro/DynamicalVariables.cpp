/*
 * DynamicalVariables.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */
#include <stdlib.h>

#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"
#include "edu/osu/rhic/harness/lattice/LatticeParameters.h"
#include "edu/osu/rhic/trunk/hydro/EnergyMomentumTensor.h"
#include "edu/osu/rhic/trunk/hydro/FullyDiscreteKurganovTadmorScheme.h" // for ghost cells 

namespace rhic {
	dynamical_variables::dynamical_variables(unsigned con_law, pi_comp, vmu_comp, pimunu_comp) {
		number_conservation_laws = con_law;
		number_pi_components = pi_comp;
		number_propagated_vmu_components = vmu_comp;
		number_propagated_pimunu_components = pimunu_comp;
		number_dissipative_currents = number_pi_components + num_propagated_vmu_components
		                              + number_propagated_pimunu_components;
		number_conserved_variables = number_conservation_laws + number_dissipative_currents; 
		is_ideal = (number_dissipative_currents == 0);
	}

	inline int dynamical_variables::col_maj_linear_index(int i, int j, int k, int nx, int ny) {
		return i + nx * (j + ny * k);
	}

	void dynamical_variables::set_conserved_vars(lattice_params &params)
	{

		int nx = params::numLatticePointsX;
		int ny = params::numLatticePointsY;
		int nz = params::numLatticePointsRapidity;
		int ncx = params::numComputationalLatticePointsX;
		int ncy = params::numComputationalLatticePointsY;

		for (int k = N_GHOST_CELLS_M; k < nz+N_GHOST_CELLS_M; ++k) {
			for (int j = N_GHOST_CELLS_M; j < ny+N_GHOST_CELLS_M; ++j) {
				for (int i = N_GHOST_CELLS_M; i < nx+N_GHOST_CELLS_M; ++i) {
					int s = col_maj_linear_index(i, j, k, ncx, ncy);

					precision_t ut_s = u[t][s];
					precision_t ux_s = u[x][s];
					precision_t uy_s = u[y][s];
					precision_t un_s = u[n][s];
					precision_t e_s = e[s];
					precision_t p_s = p[s];
					precision_t pitt_s = q.pi[t, t][s];
					precision_t pitx_s = q.pi[t, x][s];
					precision_t pity_s = q.pi[t, y][s];
					precision_t pitn_s = q.pi[t, n][s];
					precision_t Pi_s = q.Pi[s];
				
					q.t[t][s] = Ttt(e_s, p_s + Pi_s, ut_s, pitt_s);
					q.t[x][s] = Ttx(e_s, p_s + Pi_s, ut_s, ux_s, pitx_s);
					q.t[y][s] = Tty(e_s, p_s + Pi_s, ut_s, uy_s, pity_s);
					q.t[n][s] = Ttn(e_s, p_s + Pi_s, ut_s, un_s, pitn_s);
				}
			}
		}
	}

	void dynamical_variables::set_ghost_cells(lattice_parameters &params)
	{
		setGhostCellsKernelI(params);
		setGhostCellsKernelJ(params);
		setGhostCellsKernelK(params);
	}

	void dynamical_variables::set_ghost_cell_vars()
	{
		e[s] = e[sBC];
		p[s] = p[sBC];
		u[t][s] = u[t][sBC];
		u[x][s] = u[x][sBC];
		u[y][s] = u[y][sBC];
		u[n][s] = u[n][sBC];	
		q.t[t][s] = q.t[t][sBC];
		q.t[x][s] = q.t[x][sBC];
		q.t[y][s] = q.t[y][sBC];
		q.t[n][s] = q.t[n][sBC];
		// set \pi^\mu\nu ghost cells if evolved
		q.pitt[s] = q.pi[t,t][sBC];
		q.pi[t,x][s] = q.pi[t,x][sBC];
		q.pi[t,y][s] = q.pi[t,y][sBC];
		q.pi[t,n][s] = q.pi[t,n][sBC];
		q.pi[x,x][s] = q.pi[x,x][sBC];
		q.pi[x,y][s] = q.pi[x,y][sBC];
		q.pi[x,n][s] = q.pi[x,n][sBC];
		q.pi[y,y][s] = q.pi[y,y][sBC];
		q.pi[y,n][s] = q.pi[y,n][sBC];
		q.pi[n,n][s] = q.pi[n,n][sBC];	
		// set \Pi ghost cells if evolved
		q.Pi[s] = q.Pi[sBC];	
	}

	void dynamical_variables::setGhostCellsKernelI(lattice_parameters &params) {

		int nx, ncx, ncy, ncz;
		nx = params::numLatticePointsX;
		ncx = params::numComputationalLatticePointsX;
		ncy = params::numComputationalLatticePointsY;
		ncz = params::numComputationalLatticePointsRapidity;

		int iBC, s, sBC;
		for(int j = 2; j < ncy; ++j) {
			for(int k = 2; k < ncz; ++k) {
				iBC = 2;
				for (int i = 0; i <= 1; ++i) {
					s = columnMajorLinearIndex(i, j, k, ncx, ncy);	
					sBC = columnMajorLinearIndex(iBC, j, k, ncx, ncy);
					set_ghost_cell_vars(s, sBC);
				}
				iBC = nx + 1;
				for (int i = nx + 2; i <= nx + 3; ++i) {
					s = columnMajorLinearIndex(i, j, k, ncx, ncy);
					sBC = columnMajorLinearIndex(iBC, j, k, ncx, ncy);
					set_ghost_cell_vars(s, sBC);
				}
			}
		}
	}

	void dynamical_variables::setGhostCellsKernelJ(lattice_parameters &params) {

		int ny, ncx, ncy, ncz;
		ny = params::numLatticePointsY;
		ncx = params::numComputationalLatticePointsX;
		ncy = params::numComputationalLatticePointsY;
		ncz = params::numComputationalLatticePointsRapidity;

		int jBC, s, sBC;
		for(int i = 2; i < ncx; ++i) {
			for(int k = 2; k < ncz; ++k) {
				jBC = 2;
				for (int j = 0; j <= 1; ++j) {
					s = columnMajorLinearIndex(i, j, k, ncx, ncy);
					sBC = columnMajorLinearIndex(i, jBC, k, ncx, ncy);	
					set_ghost_cell_vars(s, sBC);
				}
				jBC = ny + 1;
				for (int j = ny + 2; j <= ny + 3; ++j) {
					s = columnMajorLinearIndex(i, j, k, ncx, ncy);
					sBC = columnMajorLinearIndex(i, jBC, k, ncx, ncy);
					set_ghost_cell_vars(s, sBC);		
				}
			}
		}
	}

	void dynamical_variables::setGhostCellsKernelK(lattice_parameters &params) {

		int nz, ncx, ncy;
		nz = params::numLatticePointsRapidity;
		ncx = params::numComputationalLatticePointsX;
		ncy = params::numComputationalLatticePointsY;

		int kBC, s, sBC;
		for(int i = 2; i < ncx; ++i) {
			for(int j = 2; j < ncy; ++j) {
				kBC = 2;
				for (int k = 0; k <= 1; ++k) {
					s = columnMajorLinearIndex(i, j, k, ncx, ncy);
					sBC = columnMajorLinearIndex(i, j, kBC, ncx, ncy);
					set_ghost_cell_vars(s, sBC);
				}
				kBC = nz + 1;
				for (int k = nz + 2; k <= nz + 3; ++k) {
					s = columnMajorLinearIndex(i, j, k, ncx, ncy);
					sBC = columnMajorLinearIndex(i, j, kBC, ncx, ncy);
					set_ghost_cell_vars(s,sBC);
				}
			}
		}
	}

	void dynamical_variables::set_current_conserved_vars()
	{
		std::swap(q, Q);
	}

	void dynamical_variables::swap_fluid_velocity(fluid_velocity operand1, fluid_velocity operand2)
	{
		std::swap(operand1, operand2);
	}

}
