/*
 *  Methods to impose the Dirac gauge: divergence-free condition.
 *
 *    (see file sym_tensor.h for documentation).
 *
 */

/*
 *   Copyright (c) 2006  Jerome Novak
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

char sol_Dirac_A_poisson_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2008/08/27 09:01:27  jl_cornou
 * Methods for solving Dirac systems for divergence free vectors
 *
 * Revision 1.2  2006/10/24 13:03:19  j_novak
 * New methods for the solution of the tensor wave equation. Perhaps, first
 * operational version...
 *
 * Revision 1.1  2006/09/05 15:38:45  j_novak
 * The fuctions sol_Dirac... are in a seperate file, with new parameters to
 * control the boundary conditions.
 *
 *
 * $Header$
 *
 */


// C++ headers
#include "headcpp.h"

// C headers
#include <stdlib.h>

// Lorene headers
#include "metric.h"
#include "nbr_spx.h"
#include "utilitaires.h"

// C headers
#include <assert.h>
#include <math.h>

// Lorene headers
#include "tensor.h"
#include "diff.h"
#include "proto.h"
#include "param.h"

//----------------------------------------------------------------------------------
//
//                               sol_Dirac_A
//
//----------------------------------------------------------------------------------

void Vector_divfree::sol_Dirac_A_poisson(const Scalar& aaa, Scalar& tilde_vr, Scalar& tilde_eta,
				   const Param* par_bc) const {


    Scalar source1 = -aaa.lapang();
    Scalar rvr = source1.poisson_tau();
    //rvr = rvr - rvr.val_grid_point(0,0,0,0);
    Scalar source2 = aaa.dsdr();
    source2.mult_r_dzpuis(2);
    source2 += 3*aaa;
    Scalar reta = source2.poisson_tau();
    //reta = reta - reta.val_grid_point(0,0,0,0);
    rvr.div_r();
    tilde_vr = rvr ;
    reta.div_r();
    tilde_eta = reta ;
     

} 
