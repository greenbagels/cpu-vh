/*
 * DynamicalVariables.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef DYNAMICALVARIABLES_H_
#define DYNAMICALVARIABLES_H_

#include <array>
#include <vector>
#include "LatticeParameters.h"
#include "util.h"

namespace rhic
{
	using precision_t = double;
	
	class dynamical_variables
	{
		public:
			dynamical_variables(unsigned con_law, pi_comp, vmu_comp, pimunu_comp);
			int col_maj_linear_index(int i, int j, int k, int nx, int ny);
			void set_conserved_vars(lattice_params &params);
			void set_current_conserved_vars();
			void swap_fluid_velocity(fluid_velocity &, fluid_velocity &);
			void set_ghost_cells(lattice_params &params);
			void setGhostCellsKernelI(lattice_params &params);
			void setGhostCellsKernelJ(lattice_params &params);
			void setGhostCellsKernelK(lattice_params &params);

		private:
			template<typename T>
			struct conserved_variables {
				std::array<std::vector<T>,4> t;
				util::diagonal_matrix<std::vector<T>> pi(4,4);
				std::vector<T> Pi;
			}

			enum coord {
				t,
				x,
				y,
				z
			};


			using fluid_velocity = std::array<std::vector<precision_t>, 4>;

			unsigned number_conservation_laws;
			unsigned number_pi_components;
			unsigned number_propagated_vmu_components;
			unsigned number_propagated_pimunu_components;
			unsigned number_dissipative_currents;
			unsigned number_conserved_variables;
			bool is_ideal;

			conserved_variables<precision_t> q, Q, qS;
			fluid_velocity u, up, uS, uSS;
			std::vector<precision_t> e, p;
	}
}
#endif /* DYNAMICALVARIABLES_H_ */
