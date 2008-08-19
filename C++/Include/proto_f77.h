/*
 *  Prototypes of Fortran 77 external functions
 *
 * (system dependent)
 *
 */

/*
 *   Copyright (c) 2002  Eric Gourgoulhon
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
 *
 *   LORENE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with LORENE; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifndef __PROTO_F77_H_
#define __PROTO_F77_H_

/*
 * $Id$
 * $Log$
 * Revision 1.4  2008/08/19 06:41:59  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
 * Revision 1.3  2004/12/29 12:27:35  j_novak
 * permute is now a Itbl* which array is sent directly to the LAPACK routines.
 * It is now possible to solve a general system (i.e. even if the Matrice
 * is not in a banded form).
 *
 * Revision 1.2  2002/09/09 13:54:20  e_gourgoulhon
 *
 * Change of name of the Fortran subroutines
 * 	poisson2d -> poiss2d
 * 	poisson2di -> poiss2di
 * to avoid any clash with Map::poisson2d and Map::poisson2di
 *
 * Revision 1.1  2002/09/09 13:00:39  e_gourgoulhon
 * Modification of declaration of Fortran 77 prototypes for
 * a better portability (in particular on IBM AIX systems):
 * All Fortran subroutine names are now written F77_* and are
 * defined in the new file C++/Include/proto_f77.h.
 *
 *
 *
 * $Header$
 *
 */

// FFT
// ---
#ifdef __aix
        #define F77_fft991 fft991
        #define F77_fax fax
        #define F77_fftrig fftrig
#else
        #define F77_fft991 fft991_
        #define F77_fax fax_
        #define F77_fftrig fftrig_
#endif

extern "C" {
void F77_fft991(double [], double [], double [], int [], int *, int *, int *, int *, int *) ;
void F77_fax(int *, int*, int *) ;
void F77_fftrig(double [], int *, int *) ;
}

// Prototypage des routines BLAS et LAPACK
// ---------------------------------------
#ifdef __aix
        #define F77_dswap dswap
        #define F77_dgbtrf dgbtrf
        #define F77_dgbtrs dgbtrs
        #define F77_dgetrf dgetrf
        #define F77_dgetrs dgetrs
        #define F77_dgeev dgeev
        #define F77_dgesv dgesv
#else
        #define F77_dswap dswap_
        #define F77_dgbtrf dgbtrf_
        #define F77_dgbtrs dgbtrs_
        #define F77_dgetrf dgetrf_
        #define F77_dgetrs dgetrs_
        #define F77_dgeev dgeev_
        #define F77_dgesv dgesv_
#endif

extern "C" {
void F77_dswap(int*, double [], int*, double [], int*) ;
void F77_dgbtrf(int*, int*, int*, int*, double[], int*, int[], int *) ;
void F77_dgbtrs(const char*, int*, int*, int*, int*,
			double[], int*, int[], double [], int*, int *) ;
void F77_dgetrf(int*, int*, double[], int*, int[], int *) ;
void F77_dgetrs(const char*, int*, int*, double[], int*, int[], double [], int*, int* ) ;
void F77_dgeev(const char*, const char*, int*, double[], int*, double[], double[],
			double[], int*, double[], int*, double[], int*, int*) ;
void F77_dgesv(int *, int *, double [], int *, int [], double[], int *, int *) ;
}


// Fortran routines for 2-D computations
// -------------------------------------
#ifdef __aix
        #define F77_integrale2d integrale2d
        #define F77_poiss2d poiss2d
        #define F77_poiss2di poiss2di
#else
        #define F77_integrale2d integrale2d_
        #define F77_poiss2d poiss2d_
        #define F77_poiss2di poiss2di_
#endif

extern "C" {
void F77_integrale2d(int[], int*, int*, int*,  double[], double[],
			     double*) ;
void F77_poiss2d(int[], int*, int*, int*, int[], double[], double[],
			   double[], double*, double[]) ;
void F77_poiss2di(int[], int*, int*, int*, int[], double[],
			    double[], double[]) ;
}

// Miscellaneous
// -------------
#ifdef __aix
        #define F77_huntm huntm
        #define F77_insmts insmts
#else
        #define F77_huntm huntm_
        #define F77_insmts insmts_
#endif

extern "C" {
void F77_huntm(double[], int*, double*, int*) ;
void F77_insmts(int [], int *, double [], double [], double [], double []) ;
}

#endif
