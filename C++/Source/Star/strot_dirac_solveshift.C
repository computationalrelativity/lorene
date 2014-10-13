/*
 *  Solution of the shift equation for rotating stars 
 *  in Dirac gauge.
 *
 *    (see file star_rot_dirac.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005 Lap-Ming Lin & Jerome Novak
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

char strot_dirac_solveshift_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:53:40  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2005/02/17 17:31:12  f_limousin
 * Change the name of some quantities to be consistent with other classes
 * (for instance nnn is changed to nn, shift to beta, beta to lnq...)
 *
 * Revision 1.1  2005/01/31 08:51:48  j_novak
 * New files for rotating stars in Dirac gauge (still under developement).
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "star_rot_dirac.h"
#include "unites.h"

namespace Lorene {
void Star_rot_Dirac::solve_shift(Vector& beta_new) const {

    using namespace Unites ;
    
    const Metric_flat& mets = mp.flat_met_spher() ;
    
    const Vector& dln_psi = ln_psi.derive_cov(mets) ; // D_i ln(Psi)
    const Vector& dnn = nn.derive_cov(mets) ;         // D_i N

    Vector source_beta = 2.* contract(aa, 1, 
                                   dnn - 6.*nn * dln_psi, 0) ;

    source_beta += 2.* nn * ( 2.*qpig* psi4 * j_euler 
                        - contract(tgamma.connect().get_delta(), 1, 2, 
                                   aa, 0, 1) ) ;
            
    Vector vtmp = contract(hh, 0, 1, 
                           beta.derive_cov(mets).derive_cov(mets), 1, 2)
                + 0.3333333333333333*
                  contract(hh, 1, beta.divergence(mets).derive_cov(mets), 0) ; 
    vtmp.inc_dzpuis() ; // dzpuis: 3 -> 4
                    
    source_beta -= vtmp ; 
    source_beta.set(1).set_dzpuis(4) ; //## these components are null
    source_beta.set(2).set_dzpuis(4) ; //## in axial symmetry

    Vector_divfree sou_beta_df = source_beta.div_free(mets) ;

    beta_new = sou_beta_df.poisson() ;
    beta_new.set(1) = 0 ; //## these components are null
    
    beta_new.set(2) = 0 ; //## in axial symmetry



}
}
