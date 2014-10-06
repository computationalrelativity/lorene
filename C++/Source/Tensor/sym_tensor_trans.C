/*
 *  Methods of class Sym_tensor_trans
 *
 *   (see file sym_tensor.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003-2004 Eric Gourgoulhon & Jerome Novak
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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


char sym_tensor_trans_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.17  2014/10/06 15:13:19  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.16  2008/12/03 10:20:00  j_novak
 * Modified output.
 *
 * Revision 1.15  2006/09/15 08:48:01  j_novak
 * Suppression of some messages in the optimized version.
 *
 * Revision 1.14  2005/01/03 12:37:08  j_novak
 * Sym_tensor_trans::trace_from_det_one : modified the test on hijtt to
 * be compatible with older compilers.
 *
 * Revision 1.13  2005/01/03 08:36:36  f_limousin
 * Come back to the previous version.
 *
 * Revision 1.12  2005/01/03 08:13:26  f_limousin
 * The first argument of the function trace_from_det_one(...) is now
 * a Sym_tensor_trans instead of a Sym_tensor_tt (because of a
 * compilation error with some compilators).
 *
 * Revision 1.11  2004/12/28 14:21:48  j_novak
 * Added the method Sym_tensor_trans::trace_from_det_one
 *
 * Revision 1.10  2004/05/25 15:07:12  f_limousin
 * Add parameters in argument of the function tt_part for the case
 * of a Map_et.
 *
 * Revision 1.9  2004/03/30 14:01:19  j_novak
 * Copy constructors and operator= now copy the "derived" members.
 *
 * Revision 1.8  2004/03/30 08:01:16  j_novak
 * Better treatment of dzpuis in mutators.
 *
 * Revision 1.7  2004/03/29 16:13:07  j_novak
 * New methods set_longit_trans and set_tt_trace .
 *
 * Revision 1.6  2004/03/03 13:22:14  j_novak
 * The case where dzpuis = 0 is treated in tt_part().
 *
 * Revision 1.5  2004/02/18 18:48:39  e_gourgoulhon
 * Method trace() renamed the_trace() to avoid any confusion with
 * the new method Tensor::trace().
 *
 * Revision 1.4  2004/02/09 12:57:13  e_gourgoulhon
 * First implementation of method tt_part().
 *
 * Revision 1.3  2004/01/04 20:52:45  e_gourgoulhon
 * Added assignement (operator=) to a Tensor_sym.
 *
 * Revision 1.2  2003/10/28 21:24:52  e_gourgoulhon
 * Added new methods trace() and tt_part().
 *
 * Revision 1.1  2003/10/27 10:50:54  e_gourgoulhon
 * First version.
 *
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>

// Headers Lorene
#include "metric.h"
#include "param.h"

			//--------------//
			// Constructors //
			//--------------//

// Standard constructor 
// --------------------
Sym_tensor_trans::Sym_tensor_trans(const Map& map, const Base_vect& triad_i,
		const Metric& met) 
	: Sym_tensor(map, CON, triad_i),
	  met_div(&met) {
				
	set_der_0x0() ;

}

// Copy constructor
// ----------------
Sym_tensor_trans::Sym_tensor_trans (const Sym_tensor_trans& source)
	: Sym_tensor(source), 
	  met_div(source.met_div) {
    
	set_der_0x0() ;

	if (source.p_trace != 0x0) p_trace = new Scalar( *(source.p_trace) ) ; 
	if (source.p_tt != 0x0) p_tt = new Sym_tensor_tt( *(source.p_tt) ) ; 
	
}   


// Constructor from a file
// -----------------------
Sym_tensor_trans::Sym_tensor_trans(const Map& mapping, const Base_vect& triad_i, 
	const Metric& met, FILE* fd) 
	: Sym_tensor(mapping, triad_i, fd),
	  met_div(&met) {

	set_der_0x0() ;
}

			//--------------//
			//  Destructor  //
			//--------------//

Sym_tensor_trans::~Sym_tensor_trans() {

  Sym_tensor_trans::del_deriv() ;	// in order not to follow the virtual aspect
  									// of del_deriv()

}



			//-------------------//
			// Memory managment  //
			//-------------------//

void Sym_tensor_trans::del_deriv() const {

	if (p_trace != 0x0) delete p_trace ; 
	if (p_tt != 0x0) delete p_tt ; 
	
	set_der_0x0() ;
	Sym_tensor::del_deriv() ;

}

void Sym_tensor_trans::set_der_0x0() const {

	p_trace = 0x0 ; 
	p_tt = 0x0 ; 

}


			//--------------//
			//  Assignment  //
			//--------------//

void Sym_tensor_trans::operator=(const Sym_tensor_trans& source) {
    
    // Assignment of quantities common to all the derived classes of Sym_tensor
	Sym_tensor::operator=(source) ; 

    // Assignment of proper quantities of class Sym_tensor_trans
	assert(met_div == source.met_div) ; 
	
	del_deriv() ; 
	
	if (source.p_trace != 0x0) p_trace = new Scalar( *(source.p_trace) ) ; 
	if (source.p_tt != 0x0) p_tt = new Sym_tensor_tt( *(source.p_tt) ) ; 
	
}


void Sym_tensor_trans::operator=(const Sym_tensor& source) {
    
    // Assignment of quantities common to all the derived classes of Sym_tensor
	Sym_tensor::operator=(source) ; 

  	// The metric which was set by the constructor is kept
	
	del_deriv() ; 	
}


void Sym_tensor_trans::operator=(const Tensor_sym& source) {
    
    // Assignment of quantities common to all the derived classes of Sym_tensor
    Sym_tensor::operator=(source) ; 

    // The metric which was set by the constructor is kept
	
    del_deriv() ; 	
}



void Sym_tensor_trans::operator=(const Tensor& source) {
    
    // Assignment of quantities common to all the derived classes of Sym_tensor
	Sym_tensor::operator=(source) ; 

  	// The metric which was set by the constructor is kept
	
	del_deriv() ; 	
}

void Sym_tensor_trans::set_tt_trace(const Sym_tensor_tt& htt, 
				    const Scalar& htrace, Param* par ) {
  
  assert (met_div == &htt.get_met_div() ) ;

  assert (htrace.check_dzpuis(4)) ;

  Scalar pot (*mp) ;
  if (dynamic_cast<const Map_af*>(mp) != 0x0) {
      pot = htrace.poisson() ;
  }
  else {
      pot = 0. ;
      htrace.poisson(*par, pot) ;
  }

  Sym_tensor tmp = (pot.derive_con(*met_div)).derive_con(*met_div) ; 
  tmp.dec_dzpuis() ;  
        
  *this = htrace * met_div->con() ; 
  dec_dzpuis(2) ; // this has now dzpuis = 2
  *this = htt + 0.5 * ( *this - tmp ) ;
  
  del_deriv() ;

  p_trace = new Scalar( htrace ) ;
  p_tt = new Sym_tensor_tt( htt ) ;

}


			//-----------------------------//
			//    Computational methods    //
			//-----------------------------//

const Scalar& Sym_tensor_trans::the_trace() const {

	if (p_trace == 0x0) {   // a new computation is necessary
        
		assert( (type_indice(0)==CON) && (type_indice(1)==CON) ) ; 
        p_trace = new Scalar( trace(*met_div) ) ; 
		
	}
	
	return *p_trace ; 

}


const Sym_tensor_tt& Sym_tensor_trans::tt_part(Param* par) const {
  
  if (p_tt == 0x0) {   // a new computation is necessary

    int dzp = the_trace().get_dzpuis() ;

    assert((dzp == 0) || (dzp == 4)) ;

    p_tt = new Sym_tensor_tt(*mp, *triad, *met_div) ; 

    Scalar pot (*mp) ;        
    if (dynamic_cast<const Map_af*>(mp) != 0x0) {
	pot =  the_trace().poisson() ; 
    }
    else {
	pot = 0. ;
	the_trace().poisson(*par, pot) ; 
    }

    Sym_tensor tmp = (pot.derive_con(*met_div)).derive_con(*met_div) ; 
    (dzp == 4) ? tmp.inc_dzpuis() : tmp.dec_dzpuis(3) ;  //## to be improved ?
        
    *p_tt = *this - 0.5 * ( the_trace() * met_div->con() - tmp ) ; 
		
  }
	
  return *p_tt ; 

}


void Sym_tensor_trans::trace_from_det_one(const Sym_tensor_tt& hijtt, 
					   double precis, int it_max) {
    
#ifndef NDEBUG
  const Sym_tensor_trans* ptmp = 
      dynamic_cast<const Sym_tensor_trans*>(&hijtt) ;
  assert (ptmp != 0x0) ;
  assert (ptmp != this) ;
  for (int i=0; i<hijtt.get_n_comp(); i++)
      assert(hijtt.cmp[i]->check_dzpuis(2)) ;
#endif
    assert( (precis > 0.) && (it_max > 0) ) ;
    assert (met_div == &hijtt.get_met_div() ) ;
    
    Sym_tensor_trans& hij = *this ;
    hij = hijtt ; //initialization

    // The trace h = f_{ij} h^{ij} :
    Scalar htrace(*mp) ;
        
    // Value of h at previous step of the iterative procedure below :
    Scalar htrace_prev(*mp) ;
    htrace_prev.set_etat_zero() ;   // initialisation to zero
    
    for (int it=0; it<=it_max; it++) {
      
        // Trace h from the condition det(f^{ij} + h^{ij}) = det f^{ij} :
      
        htrace = hij(1,1) * hij(2,3) * hij(2,3) 
	    + hij(2,2) * hij(1,3) * hij(1,3) + hij(3,3) * hij(1,2) * hij(1,2)
	    - 2.* hij(1,2) * hij(1,3) * hij(2,3) 
            - hij(1,1) * hij(2,2) * hij(3,3) ;
        
        htrace.dec_dzpuis(2) ; // dzpuis: 6 --> 4
        
	htrace += hij(1,2) * hij(1,2) + hij(1,3) * hij(1,3) 
                    + hij(2,3) * hij(2,3) - hij(1,1) * hij(2,2) 
                    - hij(1,1) * hij(3,3) - hij(2,2) * hij(3,3) ;

        // New value of hij from htrace and hijtt (obtained by solving 
        //    the Poisson equation for Phi) : 

        set_tt_trace(hijtt, htrace) ; 

        double diff = max(max(abs(htrace - htrace_prev))) ;
#ifndef NDEBUG
        cout << "Sym_tensor_trans::trace_from_det_one : " 
	     << "iteration : " << it << " convergence on trace(h): " << diff << endl ;
#endif
        if (diff < precis) break ;
        else htrace_prev = htrace ;

        if (it == it_max) {
            cout << "Sym_tensor_trans::trace_from_det_one : convergence not reached \n" ;
            cout << "  for the required accuracy (" << precis << ") ! " << endl ;
            abort() ;
        }
    }

}



