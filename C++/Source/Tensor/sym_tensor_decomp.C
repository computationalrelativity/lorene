/*
 *  Methods transverse( ) and longit_pot( ) of class Sym_tensor
 *
 *    (see file sym_tensor.h for documentation).
 *
 */

/*
 *   Copyright (c) 2003-2004  Eric Gourgoulhon & Jerome Novak
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

char sym_tensor_decomp_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.12  2008/12/03 10:20:00  j_novak
 * Modified output.
 *
 * Revision 1.11  2004/06/17 06:56:42  e_gourgoulhon
 * Simplified code for method transverse (use of Vector::ope_killing).
 * Slight modif. of output in method longit_pot.
 *
 * Revision 1.10  2004/06/14 20:45:41  e_gourgoulhon
 * Added argument method_poisson to methods longit_pot and transverse.
 *
 * Revision 1.9  2004/05/25 15:05:57  f_limousin
 * Add parameters in argument of the functions transverse, longit_pot
 * for the case of Map_et.
 *
 * Revision 1.8  2004/03/30 14:01:19  j_novak
 * Copy constructors and operator= now copy the "derived" members.
 *
 * Revision 1.7  2004/03/30 08:01:16  j_novak
 * Better treatment of dzpuis in mutators.
 *
 * Revision 1.6  2004/03/29 16:13:07  j_novak
 * New methods set_longit_trans and set_tt_trace .
 *
 * Revision 1.5  2004/02/09 12:56:27  e_gourgoulhon
 * Method longit_pot: added test of the vector Poisson equation.
 *
 * Revision 1.4  2004/02/02 09:18:11  e_gourgoulhon
 * Method longit_pot: treatment of case divergence dzpuis = 5.
 *
 * Revision 1.3  2003/12/10 10:17:54  e_gourgoulhon
 * First operational version.
 *
 * Revision 1.2  2003/11/27 16:01:47  e_gourgoulhon
 * First implmentation.
 *
 * Revision 1.1  2003/11/26 21:57:03  e_gourgoulhon
 * First version; not ready yet.
 *
 *
 * $Header$
 *
 */


// C headers
#include <stdlib.h>

// Lorene headers
#include "tensor.h"
#include "metric.h"
#include "param.h"

void Sym_tensor::set_longit_trans(const Vector& v_pot, 
				  const Sym_tensor_trans& ht ) {

  assert ( v_pot.get_index_type(0) == CON ) ;

  const Metric& metre = ht.get_met_div() ;

  *this = ht + v_pot.ope_killing(metre) ; // this has dzpuis = 2
  dec_dzpuis(2) ; // which is decreased so to add *this to a flat metric

  del_deriv() ;

  set_dependance(metre) ;
  int jp = get_place_met(metre) ;
  assert ((jp>=0) && (jp<N_MET_MAX)) ;

  p_transverse[jp] = new Sym_tensor_trans(ht) ;
  p_longit_pot[jp] = new Vector( v_pot ) ;
  
}

const Sym_tensor_trans& Sym_tensor::transverse(const Metric& metre, Param* par,
            int method_poisson) const {

    set_dependance(metre) ;
    int jp = get_place_met(metre) ;
    assert ((jp>=0) && (jp<N_MET_MAX)) ;

    if (p_transverse[jp] == 0x0) { // a new computation is necessary

        assert( (type_indice(0) == CON) && (type_indice(1) == CON) ) ; 

        for (int ic=0; ic<n_comp; ic++) {
            assert(cmp[ic]->check_dzpuis(4)) ;  // dzpuis=4 is assumed
        }

        const Vector& ww = longit_pot(metre, par, method_poisson) ;	    

        Sym_tensor lww = ww.ope_killing(metre) ;  // D^i W^j + D^j W^i        
        
        lww.inc_dzpuis(2) ;                     
        
        p_transverse[jp] = new Sym_tensor_trans(*mp, *triad, metre) ;
        
        *(p_transverse[jp]) = *this - lww ; 

    }

    return *p_transverse[jp] ;
    

}


const Vector& Sym_tensor::longit_pot(const Metric& metre, Param* par,
    int method_poisson) const {

    set_dependance(metre) ;
    int jp = get_place_met(metre) ;
    assert ((jp>=0) && (jp<N_MET_MAX)) ;

    if (p_longit_pot[jp] == 0x0) {  // a new computation is necessary
        
        const Metric_flat* metf = dynamic_cast<const Metric_flat*>(&metre) ; 
        if (metf == 0x0) {
            cout << "Sym_tensor::longit_pot : the case of a non flat metric"
             << endl <<"  is not treated yet !" << endl ; 
            abort() ; 
        }
        
        Vector hhh = divergence(metre) ; 


        // If dpzuis is 5, it should be decreased to 4 for the Poisson equation:
        bool dzp5 = false ; 
        for (int i=1; i<=3; i++) {
            dzp5 = dzp5 || hhh(i).check_dzpuis(5) ;
        }
        if (dzp5) hhh.dec_dzpuis() ; 
                
	if (dynamic_cast<const Map_af*>(mp) != 0x0) 
	    p_longit_pot[jp] = new Vector( hhh.poisson(double(1), 
                                                       method_poisson) ) ; 
        else
	    p_longit_pot[jp] = new Vector( hhh.poisson(double(1), *par, 
                                                       method_poisson) ) ; 


        // Test of resolution of the vector Poisson equation:
        const Vector& ww = *(p_longit_pot[jp]) ; 

        hhh.dec_dzpuis() ; 

        Vector lapw = ww.derive_con(metre).divergence(metre) 
                        + (ww.divergence(metre)).derive_con(metre) ;
#ifndef NDEBUG
        cout << "## Sym_tensor::longit_pot : test of Poisson : \n" ;
        cout << 
        "  Max absolute error in each domain on the vector Poisson equation:\n" ;   
        maxabs(lapw - hhh) ; 

	int nz = mp->get_mg()->get_nzone() ;	    // total number of domains
		
	cout << "  Relative error in each domain on the vector Poisson equation:\n" ;
	for (int i=1; i<=3; i++){
	    cout << "   Comp. " << i << " :  " ;
	    for (int l=0; l<nz; l++){
		cout << diffrel(lapw(i),hhh(i) )(l) << " " ;
	    }
	    cout << endl ;
	}
	cout << endl ;
#endif
    }

    return *p_longit_pot[jp] ;
    

}


