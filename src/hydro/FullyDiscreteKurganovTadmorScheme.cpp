/*
 * FullyDiscreteKurganovTadmorScheme.cpp
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#include <stdlib.h>
#include <stdio.h> // for printf
#include <math.h>

#include "../hydro/FullyDiscreteKurganovTadmorScheme.h"
#include "../lattice/LatticeParameters.h"
#include "../hydro/DynamicalVariables.h"
#include "../muscl/SemiDiscreteKurganovTadmorScheme.h"
#include "../muscl/HalfSiteExtrapolation.h"
#include "../hydro/FluxFunctions.h"
#include "../hydro/SpectralRadius.h"
#include "../hydro/SourceTerms.h"
#include "../hydro/EnergyMomentumTensor.h"
#include "../hydro/HydroParameters.h"

#include "../util/FiniteDifference.h" //temp

#include <omp.h>
//#define REGULATE_BULK //define to regulate bulk pressure according to inv reynolds #, otherwise bulk is not regulated

/**************************************************************************************************************************************************\
void setNeighborCells(const PRECISION * const __restrict__ data,
PRECISION * const __restrict__ I, PRECISION * const __restrict__ J, PRECISION * const __restrict__ K, PRECISION * const __restrict__ Q,
int s, unsigned int n, int ptr, int simm, int sim, int sip, int sipp, int sjmm, int sjm, int sjp, int sjpp, int skmm, int skm, int skp, int skpp) {
			PRECISION data_ns = data[s];
			// I
			*(I+ptr) = *(data+simm);
			*(I+ptr+1) = *(data+sim);
			*(I+ptr+2) = data_ns;
			*(I+ptr+3) = *(data+sip);
			*(I+ptr+4) = *(data+sipp);
			// J
			*(J+ptr) = *(data+sjmm);
			*(J+ptr+1) = *(data+sjm);
			*(J+ptr+2) = data_ns;
			*(J+ptr+3) = *(data+sjp);
			*(J+ptr+4) = *(data+sjpp);
			// K
			*(K+ptr) = *(data+skmm);
			*(K+ptr+1) = *(data+skm);
			*(K+ptr+2) = data_ns;
			*(K+ptr+3) = *(data+skp);
			*(K+ptr+4) = *(data+skpp);
			// Q
			*(Q + n) = data_ns;
}

void eulerStepKernel(PRECISION t,
const CONSERVED_VARIABLES * const __restrict__ currrentVars, CONSERVED_VARIABLES * const __restrict__ updatedVars,
const PRECISION * const __restrict__ e, const PRECISION * const __restrict__ p,
const FLUID_VELOCITY * const __restrict__ u, const FLUID_VELOCITY * const __restrict__ up,
int ncx, int ncy, int ncz, PRECISION dt, PRECISION dx, PRECISION dz, PRECISION etabar
) {
	PRECISION I[5 * NUMBER_CONSERVED_VARIABLES], J[5* NUMBER_CONSERVED_VARIABLES], K[5 * NUMBER_CONSERVED_VARIABLES];
	PRECISION hpx[NUMBER_CONSERVED_VARIABLES], hmx[NUMBER_CONSERVED_VARIABLES],
		hpy[NUMBER_CONSERVED_VARIABLES], hmy[NUMBER_CONSERVED_VARIABLES],
		hpz[NUMBER_CONSERVED_VARIABLES],	hmz[NUMBER_CONSERVED_VARIABLES];
	PRECISION Q[NUMBER_CONSERVED_VARIABLES], S[NUMBER_CONSERVED_VARIABLES];

	for(int i = 2; i < ncx-2; ++i) {
		for(int j = 2; j < ncy-2; ++j) {
			for(int k = 2; k < ncz-2; ++k) {

				int s = columnMajorLinearIndex(i, j, k, ncx, ncy);

				// calculate neighbor cell indices;
				int sim = s-1;
				int simm = sim-1;
				int sip = s+1;
				int sipp = sip+1;

				int sjm = s-ncx;
				int sjmm = sjm-ncx;
				int sjp = s+ncx;
				int sjpp = sjp+ncx;

				int skm = s-ncx*ncy;
				int skmm = skm-ncx*ncy;
				int skp = s+ncx*ncy;
				int skpp = skp+ncx*ncy;

				int ptr = 0;
				setNeighborCells(currrentVars->ttt,I,J,K,Q,s,0,ptr,simm,sim,sip,sipp,sjmm,sjm,sjp,sjpp,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCells(currrentVars->ttx,I,J,K,Q,s,1,ptr,simm,sim,sip,sipp,sjmm,sjm,sjp,sjpp,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCells(currrentVars->tty,I,J,K,Q,s,2,ptr,simm,sim,sip,sipp,sjmm,sjm,sjp,sjpp,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCells(currrentVars->ttn,I,J,K,Q,s,3,ptr,simm,sim,sip,sipp,sjmm,sjm,sjp,sjpp,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCells(currrentVars->pitt,I,J,K,Q,s,4,ptr,simm,sim,sip,sipp,sjmm,sjm,sjp,sjpp,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCells(currrentVars->pitx,I,J,K,Q,s,5,ptr,simm,sim,sip,sipp,sjmm,sjm,sjp,sjpp,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCells(currrentVars->pity,I,J,K,Q,s,6,ptr,simm,sim,sip,sipp,sjmm,sjm,sjp,sjpp,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCells(currrentVars->pitn,I,J,K,Q,s,7,ptr,simm,sim,sip,sipp,sjmm,sjm,sjp,sjpp,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCells(currrentVars->pixx,I,J,K,Q,s,8,ptr,simm,sim,sip,sipp,sjmm,sjm,sjp,sjpp,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCells(currrentVars->pixy,I,J,K,Q,s,9,ptr,simm,sim,sip,sipp,sjmm,sjm,sjp,sjpp,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCells(currrentVars->pixn,I,J,K,Q,s,10,ptr,simm,sim,sip,sipp,sjmm,sjm,sjp,sjpp,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCells(currrentVars->piyy,I,J,K,Q,s,11,ptr,simm,sim,sip,sipp,sjmm,sjm,sjp,sjpp,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCells(currrentVars->piyn,I,J,K,Q,s,12,ptr,simm,sim,sip,sipp,sjmm,sjm,sjp,sjpp,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCells(currrentVars->pinn,I,J,K,Q,s,13,ptr,simm,sim,sip,sipp,sjmm,sjm,sjp,sjpp,skmm,skm,skp,skpp);

				flux(I, hpx, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward, &spectralRadiusX, &Fx, t);
				flux(I, hmx, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards, &spectralRadiusX, &Fx, t);
				flux(J, hpy, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward, &spectralRadiusY, &Fy, t);
				flux(J, hmy, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards, &spectralRadiusY, &Fy, t);
				flux(K, hpz, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward, &spectralRadiusZ, &Fz, t);
				flux(K, hmz, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards, &spectralRadiusZ, &Fz, t);

				loadSourceTerms(I, J, K, Q, S, u->ut, u->ux, u->uy, u->un, up->ut[s], up->ux[s], up->uy[s], up->un[s], t, e[s], p,
					i, j, k, s, ncx, ncy, ncz, dt, dx, dz, etabar);

				PRECISION result[NUMBER_CONSERVED_VARIABLES];
				for (int n = 0; n < NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) = *(Q+n) + dt * ( *(S+n) - ( *(hpx+n) - *(hmx+n) + *(hpy+n) - *(hmy+n) ) / dx - ( *(hpz+n) - *(hmz+n) )/dz );
				}

				updatedVars->ttt[s] = result[0];
				updatedVars->ttx[s] = result[1];
				updatedVars->tty[s] = result[2];
				updatedVars->ttn[s] = result[3];
				updatedVars->pitt[s] = result[4];
				updatedVars->pitx[s] = result[5];
				updatedVars->pity[s] = result[6];
				updatedVars->pitn[s] = result[7];
				updatedVars->pixx[s] = result[8];
				updatedVars->pixy[s] = result[9];
				updatedVars->pixn[s] = result[10];
				updatedVars->piyy[s] = result[11];
				updatedVars->piyn[s] = result[12];
				updatedVars->pinn[s] = result[13];
			}
		}
	}
}
/**************************************************************************************************************************************************/

/**************************************************************************************************************************************************/
void setNeighborCellsJK2(const PRECISION * const __restrict__ in, PRECISION * const __restrict__ out,
int s, int ptr, int smm, int sm, int sp, int spp
) {
	PRECISION data_ns = in[s];
	*(out + ptr		) = in[smm];
	*(out + ptr + 1) = in[sm];
	*(out + ptr + 2) = data_ns;
	*(out + ptr + 3) = in[sp];
	*(out + ptr + 4) = in[spp];
}

void eulerStepKernelSource(PRECISION t,
const CONSERVED_VARIABLES * const __restrict__ currrentVars, CONSERVED_VARIABLES * const __restrict__ updatedVars,
const PRECISION * const __restrict__ e, const PRECISION * const __restrict__ p,
const FLUID_VELOCITY * const __restrict__ u, const FLUID_VELOCITY * const __restrict__ up,
int ncx, int ncy, int ncz, PRECISION dt, PRECISION dx, PRECISION dy, PRECISION dz, PRECISION etabar
) {
	//#pragma omp parallel for simd collapse(3)
	#pragma omp parallel for collapse(3)
	for(int i = 2; i < ncx-2; ++i) {
		for(int j = 2; j < ncy-2; ++j) {
			for(int k = 2; k < ncz-2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, ncx, ncy);
				PRECISION Q[NUMBER_CONSERVED_VARIABLES];
				PRECISION S[NUMBER_CONSERVED_VARIABLES];

				Q[0] = currrentVars->ttt[s];
				Q[1] = currrentVars->ttx[s];
				Q[2] = currrentVars->tty[s];
				Q[3] = currrentVars->ttn[s];
#ifdef PIMUNU
				Q[4] = currrentVars->pitt[s];
				Q[5] = currrentVars->pitx[s];
				Q[6] = currrentVars->pity[s];
				Q[7] = currrentVars->pitn[s];
				Q[8] = currrentVars->pixx[s];
				Q[9] = currrentVars->pixy[s];
				Q[10] = currrentVars->pixn[s];
				Q[11] = currrentVars->piyy[s];
				Q[12] = currrentVars->piyn[s];
				Q[13] = currrentVars->pinn[s];
#endif
#ifdef PI
				Q[14] = currrentVars->Pi[s];
#endif

				loadSourceTerms2(Q, S, u, up->ut[s], up->ux[s], up->uy[s], up->un[s], t, e[s], p, s, ncx, ncy, ncz, etabar, dt, dx, dy, dz);

				PRECISION result[NUMBER_CONSERVED_VARIABLES];
				for (unsigned int n = 0; n < NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) = *(Q+n) + dt * ( *(S+n) );
				}

				updatedVars->ttt[s] = result[0];
				updatedVars->ttx[s] = result[1];
				updatedVars->tty[s] = result[2];
				updatedVars->ttn[s] = result[3];
#ifdef PIMUNU
				updatedVars->pitt[s] = result[4];
				updatedVars->pitx[s] = result[5];
				updatedVars->pity[s] = result[6];
				updatedVars->pitn[s] = result[7];
				updatedVars->pixx[s] = result[8];
				updatedVars->pixy[s] = result[9];
				updatedVars->pixn[s] = result[10];
				updatedVars->piyy[s] = result[11];
				updatedVars->piyn[s] = result[12];
				updatedVars->pinn[s] = result[13];
#endif
#ifdef PI
				updatedVars->Pi[s] = result[14];
#endif
			}
		}
	}
}

void eulerStepKernelX(PRECISION t,
const CONSERVED_VARIABLES * const __restrict__ currrentVars, CONSERVED_VARIABLES * const __restrict__ updatedVars,
const FLUID_VELOCITY * const __restrict__ u, const PRECISION * const __restrict__ e,
int ncx, int ncy, int ncz, PRECISION dt, PRECISION dx
) {
	//#pragma omp parallel for simd collapse(3)
	#pragma omp parallel for collapse(3)
	for(int i = 2; i < ncx-2; ++i) {
		for(int j = 2; j < ncy-2; ++j) {
			for(int k = 2; k < ncz-2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, ncx, ncy);
				PRECISION I[5 * NUMBER_CONSERVED_VARIABLES];
				PRECISION H[NUMBER_CONSERVED_VARIABLES];

				// calculate neighbor cell indices;
				int sim = s-1;
				int simm = sim-1;
				int sip = s+1;
				int sipp = sip+1;

				int ptr=0;
				setNeighborCellsJK2(currrentVars->ttt,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->ttx,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->tty,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->ttn,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
#ifdef PIMUNU
				setNeighborCellsJK2(currrentVars->pitt,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pitx,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pity,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pitn,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixx,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixy,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixn,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->piyy,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->piyn,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pinn,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
#endif
#ifdef PI
				setNeighborCellsJK2(currrentVars->Pi,I,s,ptr,simm,sim,sip,sipp);
#endif

				PRECISION result[NUMBER_CONSERVED_VARIABLES];
				flux(I, H, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward, &spectralRadiusX, &Fx, t, e[s]);
				for (unsigned int n = 0; n < NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) = - *(H+n);
				}
				flux(I, H, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards, &spectralRadiusX, &Fx, t, e[s]);
				for (unsigned int n = 0; n < NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) += *(H+n);
					*(result+n) /= dx;
				}
#ifndef IDEAL
				loadSourceTermsX(I, H, u, s, dx);
				for (unsigned int n = 0; n < 4; ++n) {
					*(result+n) += *(H+n);
					*(result+n) *= dt;
				}
#else
				for (unsigned int n = 0; n < 4; ++n) {
					*(result+n) *= dt;
				}
#endif
				for (unsigned int n = 4; n < NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) *= dt;
				}

				updatedVars->ttt[s] += result[0];
				updatedVars->ttx[s] += result[1];
				updatedVars->tty[s] += result[2];
				updatedVars->ttn[s] += result[3];
#ifdef PIMUNU
				updatedVars->pitt[s] += result[4];
				updatedVars->pitx[s] += result[5];
				updatedVars->pity[s] += result[6];
				updatedVars->pitn[s] += result[7];
				updatedVars->pixx[s] += result[8];
				updatedVars->pixy[s] += result[9];
				updatedVars->pixn[s] += result[10];
				updatedVars->piyy[s] += result[11];
				updatedVars->piyn[s] += result[12];
				updatedVars->pinn[s] += result[13];
#endif
#ifdef PI
				updatedVars->Pi[s] += result[14];
#endif
			}
		}
	}
}

void eulerStepKernelY(PRECISION t,
const CONSERVED_VARIABLES * const __restrict__ currrentVars, CONSERVED_VARIABLES * const __restrict__ updatedVars,
const FLUID_VELOCITY * const __restrict__ u, const PRECISION * const __restrict__ e,
int ncx, int ncy, int ncz, PRECISION dt, PRECISION dy
) {
	//#pragma omp parallel for simd collapse(3)
	#pragma omp parallel for collapse(3)
	for(int i = 2; i < ncx-2; ++i) {
		for(int j = 2; j < ncy-2; ++j) {
			for(int k = 2; k < ncz-2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, ncx, ncy);
				PRECISION J[5* NUMBER_CONSERVED_VARIABLES];
				PRECISION H[NUMBER_CONSERVED_VARIABLES];

				// calculate neighbor cell indices;
				int sjm = s-ncx;
				int sjmm = sjm-ncx;
				int sjp = s+ncx;
				int sjpp = sjp+ncx;

				int ptr=0;
				setNeighborCellsJK2(currrentVars->ttt,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->ttx,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->tty,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->ttn,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
#ifdef PIMUNU
				setNeighborCellsJK2(currrentVars->pitt,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pitx,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pity,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pitn,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixx,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixy,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixn,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->piyy,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->piyn,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pinn,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
#endif
#ifdef PI
				setNeighborCellsJK2(currrentVars->Pi,J,s,ptr,sjmm,sjm,sjp,sjpp);
#endif

				PRECISION result[NUMBER_CONSERVED_VARIABLES];
				flux(J, H, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward, &spectralRadiusY, &Fy, t, e[s]);
				for (unsigned int n = 0; n < NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) = - *(H+n);
				}
				flux(J, H, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards, &spectralRadiusY, &Fy, t, e[s]);
				for (unsigned int n = 0; n < NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) += *(H+n);
					*(result+n) /= dy;
				}
#ifndef IDEAL
				loadSourceTermsY(J, H, u, s, dy);
				for (unsigned int n = 0; n < 4; ++n) {
					*(result+n) += *(H+n);
					*(result+n) *= dt;
				}
#else
				for (unsigned int n = 0; n < 4; ++n) {
					*(result+n) *= dt;
				}
#endif
				for (unsigned int n = 4; n < NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) *= dt;
				}

				updatedVars->ttt[s] += result[0];
				updatedVars->ttx[s] += result[1];
				updatedVars->tty[s] += result[2];
				updatedVars->ttn[s] += result[3];
#ifdef PIMUNU
				updatedVars->pitt[s] += result[4];
				updatedVars->pitx[s] += result[5];
				updatedVars->pity[s] += result[6];
				updatedVars->pitn[s] += result[7];
				updatedVars->pixx[s] += result[8];
				updatedVars->pixy[s] += result[9];
				updatedVars->pixn[s] += result[10];
				updatedVars->piyy[s] += result[11];
				updatedVars->piyn[s] += result[12];
				updatedVars->pinn[s] += result[13];
#endif
#ifdef PI
				updatedVars->Pi[s] += result[14];
#endif
			}
		}
	}
}

void eulerStepKernelZ(PRECISION t,
const CONSERVED_VARIABLES * const __restrict__ currrentVars, CONSERVED_VARIABLES * const __restrict__ updatedVars,
const FLUID_VELOCITY * const __restrict__ u, const PRECISION * const __restrict__ e,
int ncx, int ncy, int ncz, PRECISION dt, PRECISION dz
) {
	//#pragma omp parallel for simd collapse(3)
	#pragma omp parallel for collapse(3)
	for(int i = 2; i < ncx-2; ++i) {
		for(int j = 2; j < ncy-2; ++j) {
			for(int k = 2; k < ncz-2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, ncx, ncy);
				PRECISION K[5 * NUMBER_CONSERVED_VARIABLES];
				PRECISION H[NUMBER_CONSERVED_VARIABLES];

				// calculate neighbor cell indices;
				int stride = ncx * ncy;
				int skm = s-stride;
				int skmm = skm-stride;
				int skp = s+stride;
				int skpp = skp+stride;

				int ptr=0;
				setNeighborCellsJK2(currrentVars->ttt,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->ttx,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->tty,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->ttn,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
#ifdef PIMUNU
				setNeighborCellsJK2(currrentVars->pitt,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pitx,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pity,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pitn,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixx,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixy,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixn,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->piyy,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->piyn,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pinn,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
#endif
#ifdef PI
				setNeighborCellsJK2(currrentVars->Pi,K,s,ptr,skmm,skm,skp,skpp);
#endif

				PRECISION result[NUMBER_CONSERVED_VARIABLES];
				flux(K, H, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward, &spectralRadiusZ, &Fz, t, e[s]);
				for (unsigned int n = 0; n < NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) = -*(H+n);
				}
				flux(K, H, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards, &spectralRadiusZ, &Fz, t, e[s]);
				for (unsigned int n = 0; n < NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) += *(H+n);
					*(result+n) /= dz;
				}
#ifndef IDEAL
				loadSourceTermsZ(K, H, u, s, t, dz);
				for (unsigned int n = 0; n < 4; ++n) {
					*(result+n) += *(H+n);
					*(result+n) *= dt;
				}
#else
				for (unsigned int n = 0; n < 4; ++n) {
					*(result+n) *= dt;
				}
#endif
				for (unsigned int n = 4; n < NUMBER_CONSERVED_VARIABLES; ++n) {
					*(result+n) *= dt;
				}

				updatedVars->ttt[s] += result[0];
				updatedVars->ttx[s] += result[1];
				updatedVars->tty[s] += result[2];
				updatedVars->ttn[s] += result[3];
#ifdef PIMUNU
				updatedVars->pitt[s] += result[4];
				updatedVars->pitx[s] += result[5];
				updatedVars->pity[s] += result[6];
				updatedVars->pitn[s] += result[7];
				updatedVars->pixx[s] += result[8];
				updatedVars->pixy[s] += result[9];
				updatedVars->pixn[s] += result[10];
				updatedVars->piyy[s] += result[11];
				updatedVars->piyn[s] += result[12];
				updatedVars->pinn[s] += result[13];
#endif
#ifdef PI
				updatedVars->Pi[s] += result[14];
#endif
			}
		}
	}
}
/**************************************************************************************************************************************************\

/**************************************************************************************************************************************************/
void convexCombinationEulerStepKernel(const CONSERVED_VARIABLES * const __restrict__ q, CONSERVED_VARIABLES * const __restrict__ Q,
int ncx, int ncy, int ncz
) {
	//#pragma omp parallel for simd collapse(3)
	#pragma omp parallel for collapse(3)
	for(int i = 2; i < ncx-2; ++i) {
		for(int j = 2; j < ncy-2; ++j) {
			for(int k = 2; k < ncz-2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, ncx, ncy);
				Q->ttt[s] += q->ttt[s];
				Q->ttt[s] /= 2;
				Q->ttx[s] += q->ttx[s];
				Q->ttx[s] /= 2;
				Q->tty[s] += q->tty[s];
				Q->tty[s] /= 2;
				Q->ttn[s] += q->ttn[s];
				Q->ttn[s] /= 2;
				#ifdef PIMUNU
				Q->pitt[s] += q->pitt[s];
				Q->pitt[s] /= 2;
				Q->pitx[s] += q->pitx[s];
				Q->pitx[s] /= 2;
				Q->pity[s] += q->pity[s];
				Q->pity[s] /= 2;
				Q->pitn[s] += q->pitn[s];
				Q->pitn[s] /= 2;
				Q->pixx[s] += q->pixx[s];
				Q->pixx[s] /= 2;
				Q->pixy[s] += q->pixy[s];
				Q->pixy[s] /= 2;
				Q->pixn[s] += q->pixn[s];
				Q->pixn[s] /= 2;
				Q->piyy[s] += q->piyy[s];
				Q->piyy[s] /= 2;
				Q->piyn[s] += q->piyn[s];
				Q->piyn[s] /= 2;
				Q->pinn[s] += q->pinn[s];
				Q->pinn[s] /= 2;
				#endif
				#ifdef PI
				Q->Pi[s] += q->Pi[s];
				Q->Pi[s] /= 2.0;
				#endif
			}
		}
	}
}

/**************************************************************************************************************************************************/
#ifndef IDEAL
void
regulateDissipativeCurrents(PRECISION t,
const CONSERVED_VARIABLES * const __restrict__ currentVars,
const PRECISION * const __restrict__ e, const PRECISION * const __restrict__ p,
const FLUID_VELOCITY * const __restrict__ u,
int ncx, int ncy, int ncz
) {
	//#pragma omp parallel for simd collapse(3)
	#pragma omp parallel for collapse(3)
	for(int i = 2; i < ncx-2; ++i) {
		for(int j = 2; j < ncy-2; ++j) {
			for(int k = 2; k < ncz-2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, ncx, ncy);

#ifdef PIMUNU
				PRECISION pitt = currentVars->pitt[s];
				PRECISION pitx = currentVars->pitx[s];
				PRECISION pity = currentVars->pity[s];
				PRECISION pitn = currentVars->pitn[s];
				PRECISION pixx = currentVars->pixx[s];
				PRECISION pixy = currentVars->pixy[s];
				PRECISION pixn = currentVars->pixn[s];
				PRECISION piyy = currentVars->piyy[s];
				PRECISION piyn = currentVars->piyn[s];
				PRECISION pinn = currentVars->pinn[s];
#else
				PRECISION pitt = 0.0;
				PRECISION pitx = 0.0;
				PRECISION pity = 0.0;
				PRECISION pitn = 0.0;
				PRECISION pixx = 0.0;
				PRECISION pixy = 0.0;
				PRECISION pixn = 0.0;
				PRECISION piyy = 0.0;
				PRECISION piyn = 0.0;
				PRECISION pinn = 0.0;
#endif
#ifdef PI
				PRECISION Pi = currentVars->Pi[s];
#else
				PRECISION Pi = 0;
#endif

				PRECISION ut = u->ut[s];
				PRECISION ux = u->ux[s];
				PRECISION uy = u->uy[s];
				PRECISION un = u->un[s];

				PRECISION xi0 = (PRECISION)(1.0);
				PRECISION rhomax = (PRECISION)(10.0);
				//PRECISION xi0 = (PRECISION)(0.1);
				//PRECISION rhomax = (PRECISION)(0.8);
				PRECISION t2 = t*t;
		//		PRECISION pipi = pitt*pitt-2*(pitx*pitx+pity*pity-pixy*pixy+t2*(pitn*pitn-pixn*pixn-piyn*piyn))+pixx*pixx+piyy*piyy+pinn*pinn*t2*t2;
				PRECISION pipi = pitt*pitt-2*pitx*pitx-2*pity*pity+pixx*pixx+2*pixy*pixy+piyy*piyy-2*pitn*pitn*t2+2*pixn*pixn*t2+2*piyn*piyn*t2+pinn*pinn*t2*t2;
				if(isnan(pipi)==1) printf("found pipi Nan\n");
				//PRECISION spipi = sqrt(fabs(pipi+3*Pi*Pi)); //change this to sqrt(fabs(pipi)) to remove bulk pressure regulation of shear stress
				PRECISION spipi = sqrt(fabs(pipi)); //change this to sqrt(fabs(pipi)) to remove bulk pressure regulation of shear stress
				if(isnan(spipi)==1) printf("found spipi Nan\n");
				PRECISION pimumu = pitt - pixx - piyy - pinn*t*t;
				PRECISION piu0 = -(pitn*t2*un) + pitt*ut - pitx*ux - pity*uy;
				PRECISION piu1 = -(pixn*t2*un) + pitx*ut - pixx*ux - pixy*uy;
				PRECISION piu2 = -(piyn*t2*un) + pity*ut - pixy*ux - piyy*uy;
				PRECISION piu3 = -(pinn*t2*un) + pitn*ut - pixn*ux - piyn*uy;

				PRECISION a1 = spipi/rhomax/sqrtf(e[s]*e[s]+3*p[s]*p[s]);
				if(isnan(a1)==1) printf("found a1 Nan\n");
				PRECISION a2 = pimumu/xi0/rhomax/spipi;
				PRECISION a3 = piu0/xi0/rhomax/spipi;
				PRECISION a4 = piu1/xi0/rhomax/spipi;
				PRECISION a5 = piu2/xi0/rhomax/spipi;
				PRECISION a6 = piu3/xi0/rhomax/spipi;
				PRECISION a12 = fmax(a1,a2);
				PRECISION a34 = fmax(a3,a4);
				PRECISION a56 = fmax(a5,a6);
				PRECISION a3456 = fmax(a34,a56);
				PRECISION rho = fmax(a12,a3456);

				if(isnan(rho)==1) printf("found rho Nan\n");
				PRECISION fac = tanh(rho)/rho;
				if(fabs(rho)<1.e-7) fac = 1;
				if(isnan(fac)==1) printf("found fac Nan\n");

				//regulate the shear stress
				#ifdef PIMUNU
				currentVars->pitt[s] *= fac;
				currentVars->pitx[s] *= fac;
				currentVars->pity[s] *= fac;
				currentVars->pitn[s] *= fac;
				currentVars->pixx[s] *= fac;
				currentVars->pixy[s] *= fac;
				currentVars->pixn[s] *= fac;
				currentVars->piyy[s] *= fac;
				currentVars->piyn[s] *= fac;
				currentVars->pinn[s] *= fac;
				#endif

				//regulate the bulk pressure according to it's inverse reynolds #
				#ifdef REGULATE_BULK
				PRECISION rhoBulk = abs(Pi) / sqrtf(e[s]*e[s]+3*p[s]*p[s]);
				if(isnan(rhoBulk) == 1) printf("found rhoBulk Nan\n");
                                PRECISION facBulk = tanh(rhoBulk) / rhoBulk;
                                if(fabs(rhoBulk) < 1.e-7) facBulk = 1.0;
                                if(isnan(facBulk) == 1) printf("found facBulk Nan\n");

				//regulate bulk pressure
				#ifdef PI
				currentVars->Pi[s] *= facBulk;
				#endif

				#endif

			}
		}
	}
}
#endif
/**************************************************************************************************************************************************/
#ifndef IDEAL
#define REGULATE_DISSIPATIVE_CURRENTS
#endif

void
rungeKutta2(PRECISION t, PRECISION dt, CONSERVED_VARIABLES * __restrict__ q, CONSERVED_VARIABLES * __restrict__ Q,
void * latticeParams, void * hydroParams
) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;
	int ncx = lattice->numComputationalLatticePointsX;
	int ncy = lattice->numComputationalLatticePointsY;
	int ncz = lattice->numComputationalLatticePointsRapidity;

	PRECISION dx = (PRECISION)(lattice->latticeSpacingX);
	PRECISION dy = (PRECISION)(lattice->latticeSpacingY);
	PRECISION dz = (PRECISION)(lattice->latticeSpacingRapidity);

	PRECISION etabar = (PRECISION)(hydro->shearViscosityToEntropyDensity);

	//===================================================
	// STEP 1:
	//===================================================
	eulerStepKernelSource(t, q, qS, e, p, u, up, ncx, ncy, ncz, dt, dx, dy, dz, etabar);
	eulerStepKernelX(t, q, qS, u, e, ncx, ncy, ncz, dt, dx);
	eulerStepKernelY(t, q, qS, u, e, ncx, ncy, ncz, dt, dy);
	eulerStepKernelZ(t, q, qS, u, e, ncx, ncy, ncz, dt, dz);

	t+=dt;

	setInferredVariablesKernel(qS, e, p, uS, t, latticeParams);

#ifdef REGULATE_DISSIPATIVE_CURRENTS
	regulateDissipativeCurrents(t, qS, e, p, uS, ncx, ncy, ncz);
#endif

	setGhostCells(qS, e, p, uS, latticeParams);

	//===================================================
	// STEP 2:
	//===================================================
	eulerStepKernelSource(t, qS, Q, e, p, uS, u, ncx, ncy, ncz, dt, dx, dy, dz, etabar);
	eulerStepKernelX(t, qS, Q, uS, e, ncx, ncy, ncz, dt, dx);
	eulerStepKernelY(t, qS, Q, uS, e, ncx, ncy, ncz, dt, dy);
	eulerStepKernelZ(t, qS, Q, uS, e, ncx, ncy, ncz, dt, dz);

	convexCombinationEulerStepKernel(q, Q, ncx, ncy, ncz);

	swapFluidVelocity(&up, &u);
	setInferredVariablesKernel(Q, e, p, u, t, latticeParams);

#ifdef REGULATE_DISSIPATIVE_CURRENTS
	regulateDissipativeCurrents(t, Q, e, p, u, ncx, ncy, ncz);
#endif

	setGhostCells(Q, e, p, u, latticeParams);
}
