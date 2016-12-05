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

 

/*
 * $Id$
 * $Log$
 * Revision 1.5  2016/12/05 16:18:18  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:53:45  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:20  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2009/10/23 13:18:46  j_novak
 * Minor modifications
 *
 * Revision 1.1  2008/08/27 09:01:27  jl_cornou
 * Methods for solving Dirac systems for divergence free vectors
 *
 * $Header$
 *
 */


// C headers
#include <cstdlib>
#include <cassert>
#include <cmath>

// Lorene headers
#include "metric.h"
#include "diff.h"
#include "proto.h"
#include "param.h"

//----------------------------------------------------------------------------------
//
//                               sol_Dirac_A
//
//----------------------------------------------------------------------------------

namespace Lorene {
void Vector_divfree::sol_Dirac_A_poisson(const Scalar& aaa, Scalar& tilde_vr, Scalar& tilde_eta,
				   const Param* ) const {


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
}
