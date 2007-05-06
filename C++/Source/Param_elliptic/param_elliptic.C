/*
 *   Copyright (c) 2004 Philippe Grandclement
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

char param_elliptic_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.19  2007/05/06 10:48:14  p_grandclement
 * Modification of a few operators for the vorton project
 *
 * Revision 1.18  2007/04/24 09:04:13  p_grandclement
 * Addition of an operator for the vortons
 *
 * Revision 1.17  2005/05/12 09:49:44  j_novak
 * Temptative treatment of the case where the source is null in the CED (putting
 * dzpuis to 4). May be a bad idea...
 *
 * Revision 1.16  2005/02/15 15:43:17  j_novak
 * First version of the block inversion for the vector Poisson equation (method 6).
 *
 * Revision 1.15  2004/12/23 16:30:16  j_novak
 * New files and class for the solution of the rr component of the tensor Poisson
 * equation.
 *
 * Revision 1.14  2004/08/24 09:14:49  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.13  2004/06/22 13:46:52  j_novak
 * Deplacement du conte++ dans set_pois_vect_r
 *
 * Revision 1.12  2004/06/22 08:49:59  p_grandclement
 * Addition of everything needed for using the logarithmic mapping
 *
 * Revision 1.11  2004/06/14 15:07:12  j_novak
 * New methods for the construction of the elliptic operator appearing in
 * the vector Poisson equation (acting on eta).
 *
 * Revision 1.10  2004/05/17 15:50:54  j_novak
 * Removed unused nr variables
 *
 * Revision 1.9  2004/05/17 15:42:22  j_novak
 * The method 1 of vector Poisson eq. solves directly for F^r.
 * Some bugs were corrected in the operator poisson_vect.
 *
 * Revision 1.8  2004/05/14 08:51:02  p_grandclement
 * *** empty log message ***
 *
 * Revision 1.7  2004/05/10 15:28:22  j_novak
 * First version of functions for the solution of the r-component of the
 * vector Poisson equation.
 *
 * Revision 1.6  2004/03/05 09:18:49  p_grandclement
 * Addition of operator sec_order_r2
 *
 * Revision 1.5  2004/01/15 09:15:39  p_grandclement
 * Modification and addition of the Helmholtz operators
 *
 * Revision 1.4  2004/01/07 14:36:38  p_grandclement
 * Modif mineure in Param_elliptic.set_variable
 *
 * Revision 1.3  2003/12/11 16:11:38  e_gourgoulhon
 * Changed #include <iostream.h> to #include "headcpp.h".
 *
 * Revision 1.2  2003/12/11 15:57:27  p_grandclement
 * include stdlib.h encore ...
 *
 * Revision 1.1  2003/12/11 14:48:51  p_grandclement
 * Addition of ALL (and that is a lot !) the files needed for the general elliptic solver ... UNDER DEVELOPEMENT...
 *
 * 
 * $Header$
 *
 */

#include "headcpp.h"

#include <math.h>
#include <stdlib.h>

#include "map.h"
#include "ope_elementary.h"
#include "param_elliptic.h"
#include "base_val.h" 
#include "scalar.h"
#include "proto.h"


const Map_radial& Param_elliptic::get_mp() const {

  switch (type_map) {
  case MAP_AFF:
    return *mp_af ;
    break ;
  case MAP_LOG:
    return *mp_log ;
    break ;
  default:
    cout << "Unknown mapping in Param_elliptic" << endl ;
    abort() ;
    return *mp_af ;
  }
}

double Param_elliptic::get_alpha(int l) const {
  
 switch (type_map) {
  case MAP_AFF:
    return mp_af->get_alpha()[l] ;
    break ;
  case MAP_LOG:
    return mp_log->get_alpha(l) ;
    break ;
  default:
    cout << "Unknown mapping in Param_elliptic" << endl ;
    abort() ;
    return 1 ;
  }
}

double Param_elliptic::get_beta(int l) const {
  
 switch (type_map) {
  case MAP_AFF:
    return mp_af->get_beta()[l] ;
    break ;
  case MAP_LOG:
    return mp_log->get_beta(l) ;
    break ;
  default:
    cout << "Unknown mapping in Param_elliptic" << endl ;
    abort() ;
    return 1 ;
  }
}

int Param_elliptic::get_type(int l) const {
  
 switch (type_map) {
  case MAP_AFF:
    return AFFINE ;
    break ;
  case MAP_LOG:
    return mp_log->get_type(l) ;
    break ;
  default:
    cout << "Unknown mapping in Param_elliptic" << endl ;
    abort() ;
    return 1 ;
  }
}

// Construction (By default it is set to Poisson with appropriate dzpuis...)
Param_elliptic::Param_elliptic(const Scalar& so) : var_F(so.get_mp()), var_G(so.get_mp()), 
						   done_F (so.get_mp().get_mg()->get_nzone(), 
							   so.get_mp().get_mg()->get_np(0) + 1,  
							   so.get_mp().get_mg()->get_nt(0)), 
						   done_G (so.get_mp().get_mg()->get_nzone()), 
						   val_F_plus (so.get_mp().get_mg()->get_nzone(), 
							   so.get_mp().get_mg()->get_np(0) + 1,  
							   so.get_mp().get_mg()->get_nt(0)), 
						   val_F_minus (so.get_mp().get_mg()->get_nzone(), 
							   so.get_mp().get_mg()->get_np(0) + 1,  
							    so.get_mp().get_mg()->get_nt(0)),
						   val_dF_plus (so.get_mp().get_mg()->get_nzone(), 
							   so.get_mp().get_mg()->get_np(0) + 1,  
							   so.get_mp().get_mg()->get_nt(0)), 
						   val_dF_minus (so.get_mp().get_mg()->get_nzone(), 
							   so.get_mp().get_mg()->get_np(0) + 1,  
							    so.get_mp().get_mg()->get_nt(0)),
						   val_G_plus (so.get_mp().get_mg()->get_nzone()), 
						   val_G_minus (so.get_mp().get_mg()->get_nzone()),
						   val_dG_plus (so.get_mp().get_mg()->get_nzone()), 
						   val_dG_minus (so.get_mp().get_mg()->get_nzone())
						   
						   
{

  // On passe en Ylm
  Scalar auxi(so) ;
  auxi.set_spectral_va().coef() ;
  auxi.set_spectral_va().ylm() ;

  Base_val base (auxi.get_spectral_va().base) ;
  int dzpuis = (auxi.dz_nonzero() ? auxi.get_dzpuis() : 4) ; //## to be modified??
  
  // Right now, only applicable with affine mapping
  const Map_af* map_affine = dynamic_cast <const Map_af*> (&so.get_mp()) ;
  const Map_log* map_log = dynamic_cast <const Map_log*> (&so.get_mp()) ;

  
  if ((map_affine == 0x0) && (map_log == 0x0)) {
    cout << "Param_elliptic not yet defined on this type of mapping" << endl ;
    abort() ;
  }
  else  {
    
    if (map_affine != 0x0) {
      type_map = MAP_AFF ;
      mp_af = map_affine ;
      mp_log = 0x0 ;
    }
    if (map_log != 0x0) {
      type_map = MAP_LOG ;
      mp_af = 0x0 ;
      mp_log = map_log ;
    }
    int nz = get_mp().get_mg()->get_nzone() ;
    int nbr = 0 ;
    for (int l=0 ; l<nz ; l++)
      nbr += (get_mp().get_mg()->get_np(l)+1)*
	get_mp().get_mg()->get_nt(l) ;
    
    operateurs = new Ope_elementary* [nbr] ;
    
    for (int l=0 ; l<nbr ; l++)
      operateurs[l] = 0x0 ;

      int nr ;
      int base_r, m_quant, l_quant ;
      
      int conte = 0 ;
      for (int l=0 ; l<nz ; l++) {
	
	nr = get_mp().get_mg()->get_nr(l) ;
	
	for (int k=0 ; k<get_mp().get_mg()->get_np(l)+1 ; k++)
	  for (int j=0 ; j<get_mp().get_mg()->get_nt(l) ; j++) {
	    
	    so.get_spectral_va().base.give_quant_numbers 
	      (l, k, j, m_quant, l_quant, base_r) ;
 
	    switch (type_map) {
	    case MAP_AFF: 
	      operateurs[conte] = new 
		Ope_poisson(nr, base_r, get_alpha(l), 
			    get_beta(l), l_quant, dzpuis) ;
	      break ;
	    case MAP_LOG:
	      if (mp_log->get_type(l) == AFFINE)
		operateurs[conte] = new 
		  Ope_poisson(nr, base_r, get_alpha(l), 
			      get_beta(l), l_quant, dzpuis) ;
	      else 
		operateurs[conte] = new 
		  Ope_sec_order (nr, base_r, get_alpha(l), get_beta(l), 
				 1., 2. , l_quant) ;
	      break ;
	    default :
	      cout << "Unknown mapping in Param_elliptic::Param_elliptic" 
		   << endl ; 
	      
	    }
	    conte ++ ;
	  }
      }
  }
  
  // STD VARIABLE CHANGE :
  var_F.annule_hard() ;
  var_F.set_spectral_va().set_base (so.get_spectral_va().get_base()) ;
  
  var_G.set_etat_qcq() ;
  var_G = 1 ;
  var_G.std_spectral_base() ;
  
  done_F.annule_hard() ;
  done_G.annule_hard() ;
  
  val_F_plus.set_etat_qcq() ;
  val_F_minus.set_etat_qcq() ;
  val_dF_plus.set_etat_qcq() ;
  val_dF_minus.set_etat_qcq() ;
  
  val_G_plus.set_etat_qcq() ;
  val_G_minus.set_etat_qcq() ;
  val_dG_plus.set_etat_qcq() ;
  val_dG_minus.set_etat_qcq() ;
}

Param_elliptic::~Param_elliptic () {
  
  int nbr = 0 ;
  for (int l=0 ; l<get_mp().get_mg()->get_nzone() ; l++) 
    nbr += (get_mp().get_mg()->get_np(l)+1)*get_mp().get_mg()->get_nt(l) ;

  for (int i=0 ; i<nbr ; i++)
    if (operateurs[i] != 0x0)
      delete operateurs[i] ;
  
  delete [] operateurs ;
}

void Param_elliptic::inc_l_quant (int zone) {
  
  int np, nt ;

  int conte = 0 ;
  for (int l=0 ; l<get_mp().get_mg()->get_nzone() ; l++) {
    
    np = get_mp().get_mg()->get_np(l) ;
    nt = get_mp().get_mg()->get_nt(l) ;

    for (int k=0 ; k<np+1 ; k++)
      for (int j=0 ; j<nt ; j++) { 
	if ((operateurs[conte] != 0x0) && (l==zone))
	  operateurs[conte]->inc_l_quant() ;
	conte ++ ;
      }
  }
}


void Param_elliptic::set_helmholtz_minus (int zone, double masse, Scalar& source) {
  
  source.set_spectral_va().coef() ;
  source.set_spectral_va().ylm() ;

  int nz = get_mp().get_mg()->get_nzone() ;
  int nr ;
  int conte = 0 ;
  int m_quant, base_r_1d, l_quant ;

  switch (type_map) {

  case MAP_AFF:
 
    for (int l=0 ; l<nz ; l++) {
      
      nr = get_mp().get_mg()->get_nr(l) ;
      
      for (int k=0 ; k<get_mp().get_mg()->get_np(l)+1 ; k++)
	for (int j=0 ; j<get_mp().get_mg()->get_nt(l) ; j++) {
	  if (l==zone) {
	    if (operateurs[conte] != 0x0) {	    
	      delete operateurs[conte] ;
	      source.get_spectral_va().base.give_quant_numbers 
		(l, k, j, m_quant, l_quant, base_r_1d) ;
	      operateurs[conte] = new Ope_helmholtz_minus (nr, base_r_1d, l_quant, get_alpha(l), 
							   get_beta(l) , masse) ;
	    }
	  }
	  conte ++ ;
	}
    }
    break ;
    
  case MAP_LOG :
    if (mp_log->get_type(zone) != AFFINE) {
      cout << "Operator not define with LOG mapping..." << endl ;
      abort() ;
    }
    else {
      for (int l=0 ; l<nz ; l++) {
	
	nr = get_mp().get_mg()->get_nr(l) ;
	
	for (int k=0 ; k<get_mp().get_mg()->get_np(l)+1 ; k++)
	  for (int j=0 ; j<get_mp().get_mg()->get_nt(l) ; j++) {
	    if (l==zone) {
	      if (operateurs[conte] != 0x0) {	    
		delete operateurs[conte] ;
		source.get_spectral_va().base.give_quant_numbers 
		  (l, k, j, m_quant, l_quant, base_r_1d) ;
		operateurs[conte] = new Ope_helmholtz_minus (nr, base_r_1d, l_quant,
							     get_alpha(l), get_beta(l), masse) ;
	      }
	    }
	    conte ++ ;
	  }
      }
    }
    break ;
    
  default :
    cout << "Unkown mapping in set_helmhotz_minus" << endl ;
    abort() ;
    break ;
  }
}

void Param_elliptic::set_helmholtz_plus (int zone, double masse, Scalar& source) {
  
  source.set_spectral_va().coef() ;
  source.set_spectral_va().ylm() ;

  int nz = get_mp().get_mg()->get_nzone() ;
  int nr ;
  int conte = 0 ;
  int m_quant, base_r_1d, l_quant ;

  switch (type_map) {

  case MAP_AFF:
 
    for (int l=0 ; l<nz ; l++) {
      
      nr = get_mp().get_mg()->get_nr(l) ;
      
      for (int k=0 ; k<get_mp().get_mg()->get_np(l)+1 ; k++)
	for (int j=0 ; j<get_mp().get_mg()->get_nt(l) ; j++) {
	  if (l==zone) {
	    if (operateurs[conte] != 0x0) {	
		int old_base = operateurs[conte]->get_base_r() ;
	    // PROVISOIRE, DANS LE NOYAU SEUL LE CAS SPHERIQUE EST IMPLEMENTE
	    if ((old_base != R_CHEBI)) {    
	      delete operateurs[conte] ;
	      source.get_spectral_va().base.give_quant_numbers 
		(l, k, j, m_quant, l_quant, base_r_1d) ;
	      operateurs[conte] = new Ope_helmholtz_plus (nr, base_r_1d, l_quant, get_alpha(l), 
							   get_beta(l) , masse) ;
	    }
		}
	  }
	  conte ++ ;
	}
    }
    break ;
    
  case MAP_LOG :
    if (mp_log->get_type(zone) != AFFINE) {
      cout << "Operator not define with LOG mapping..." << endl ;
      abort() ;
    }
    else {
      for (int l=0 ; l<nz ; l++) {
	
	nr = get_mp().get_mg()->get_nr(l) ;
	
	for (int k=0 ; k<get_mp().get_mg()->get_np(l)+1 ; k++)
	  for (int j=0 ; j<get_mp().get_mg()->get_nt(l) ; j++) {
	    if (l==zone) {
	      if (operateurs[conte] != 0x0) {
		int old_base = operateurs[conte]->get_base_r() ;
	    // PROVISOIRE, DANS LE NOYAU SEUL LE CAS SPHERIQUE EST IMPLEMENTE
	    	if ((old_base != R_CHEBI)) {	    
			delete operateurs[conte] ;
			source.get_spectral_va().base.give_quant_numbers 
		  		(l, k, j, m_quant, l_quant, base_r_1d) ;
			operateurs[conte] = new Ope_helmholtz_plus (nr, base_r_1d, l_quant,
							     get_alpha(l), get_beta(l), masse) ;
	      }
		}
	    }
	    conte ++ ;
	  }
      }
    }
    break ;
    
  default :
    cout << "Unkown mapping in set_helmhotz_minus" << endl ;
    abort() ;
    break ;
  }
}

void Param_elliptic::set_poisson_vect_r(int zone, bool only_l_zero) {
 
  if (type_map != MAP_AFF) {
    cout << "set_poisson_vect_r only defined for an affine mapping..." << endl ;
    abort() ;
  }
  else {

    int nz = get_mp().get_mg()->get_nzone() ;
    
    int nr ;
    
    int conte = 0 ;
    for (int l=0 ; l<nz ; l++) {
      
      nr = get_mp().get_mg()->get_nr(l) ;
      
      for (int k=0 ; k<get_mp().get_mg()->get_np(l)+1 ; k++)
	for (int j=0 ; j<get_mp().get_mg()->get_nt(l) ; j++) {
	  if ((operateurs[conte] != 0x0) && (l==zone)) {
	    int old_base = operateurs[conte]->get_base_r() ;
	    Ope_poisson* op_pois = 
	      dynamic_cast<Ope_poisson*>(operateurs[conte]) ;
	    assert (op_pois !=0x0) ;
	    int lq_old = op_pois->get_lquant() ;
	    int dzp = op_pois->get_dzpuis() ;
	    
	    delete operateurs[conte] ;
	    if (type_map == MAP_AFF) {
		if ((!only_l_zero)||(lq_old == 0)) {
		    operateurs[conte] = new Ope_pois_vect_r(nr, old_base,get_alpha(l), 
							    get_beta(l), lq_old, dzp) ;
		}
		else {
		    operateurs[conte] = 0x0 ;
		}
	    }
	    else
	      operateurs[conte] = 0x0 ;
	  }
      conte ++ ;
	}
    }
  }
}

void Param_elliptic::set_poisson_vect_eta(int zone) {
 
  int nz = get_mp().get_mg()->get_nzone() ;
  
  int conte = 0 ;
  for (int l=0 ; l<nz ; l++) {
    
    bool ced = (get_mp().get_mg()->get_type_r(l) == UNSURR ) ;
    for (int k=0 ; k<get_mp().get_mg()->get_np(l)+1 ; k++)
      for (int j=0 ; j<get_mp().get_mg()->get_nt(l) ; j++) {
	if ((operateurs[conte] != 0x0) && (l==zone)) {
	  Ope_poisson* op_pois = 
	    dynamic_cast<Ope_poisson*>(operateurs[conte]) ;
	  assert (op_pois !=0x0) ;
	  int lq_old = op_pois->get_lquant() ;
	  if (lq_old == 0) {
	    delete operateurs[conte] ;
	    operateurs[conte] = 0x0 ;
	  }
	  else 
	    ced ? op_pois->inc_l_quant() : op_pois->dec_l_quant() ;
	}
	conte ++ ;
      }
  }
}

void Param_elliptic::set_poisson_tens_rr(int zone) {
 
    if (type_map != MAP_AFF) {
	cout << "set_poisson_tens_rr only defined for an affine mapping..." 
	     << endl ;
	abort() ;
    }
    else {
	
	int nz = get_mp().get_mg()->get_nzone() ;
	
	int nr ;
	
	int conte = 0 ;
	for (int l=0 ; l<nz ; l++) {
	    
	    nr = get_mp().get_mg()->get_nr(l) ;
	    
	    for (int k=0 ; k<get_mp().get_mg()->get_np(l)+1 ; k++)
		for (int j=0 ; j<get_mp().get_mg()->get_nt(l) ; j++) {
		    if ((operateurs[conte] != 0x0) && (l==zone)) {
			int old_base = operateurs[conte]->get_base_r() ;
			Ope_poisson* op_pois = 
			    dynamic_cast<Ope_poisson*>(operateurs[conte]) ;
			assert (op_pois !=0x0) ;
			int lq_old = op_pois->get_lquant() ;
			int dzp = op_pois->get_dzpuis() ;
			
			delete operateurs[conte] ;
			if (lq_old >= 2) {
			    operateurs[conte] = new Ope_pois_tens_rr
				(nr, old_base,get_alpha(l), get_beta(l), lq_old, dzp) ;
			}
			else
			    operateurs[conte] = 0x0 ;
		    }
		    conte ++ ;
		}
	}
    }
}

void Param_elliptic::set_sec_order_r2 (int zone, double a, double b, double c){
 
   
  int nz = get_mp().get_mg()->get_nzone() ;
  int nr ;
  int conte = 0 ;

  switch (type_map) {

  case MAP_AFF:
 
    for (int l=0 ; l<nz ; l++) {
      
      nr = get_mp().get_mg()->get_nr(l) ;
      
      for (int k=0 ; k<get_mp().get_mg()->get_np(l)+1 ; k++)
	for (int j=0 ; j<get_mp().get_mg()->get_nt(l) ; j++) {
	  if ((operateurs[conte] != 0x0) && (l==zone)) {
	    int old_base = operateurs[conte]->get_base_r() ;
	    // PROVISOIRE, DANS LE NOYAU SEUL LE CAS SPHERIQUE EST IMPLEMENTE
	    if ((old_base != R_CHEBI)) {
	      delete operateurs[conte] ;
	      operateurs[conte] = new Ope_sec_order_r2 (nr, old_base, 
							   get_alpha(l), 
							   get_beta(l), a, b, c) ;
	    }
	  }
	  conte ++ ;
	}
    }
    break ;
    
  case MAP_LOG :
    if (mp_log->get_type(zone) != AFFINE) {
      cout << "Operator not define with LOG mapping..." << endl ;
      abort() ;
    }
    else {
      for (int l=0 ; l<nz ; l++) {
	
	nr = get_mp().get_mg()->get_nr(l) ;
	
	for (int k=0 ; k<get_mp().get_mg()->get_np(l)+1 ; k++)
	  for (int j=0 ; j<get_mp().get_mg()->get_nt(l) ; j++) {
	    if ((operateurs[conte] != 0x0) && (l==zone)) {
	      int old_base = operateurs[conte]->get_base_r() ;
	      // PROVISOIRE, DANS LE NOYAU SEUL LE CAS SPHERIQUE EST IMPLEMENTE
	      if ((old_base != R_CHEBI)) {
		delete operateurs[conte] ;
		operateurs[conte] = new Ope_sec_order_r2 (nr, old_base, 
							     get_alpha(l), 
							     get_beta(l), a, b, c) ;
	      }
	    }
	    conte ++ ;
	  }
      }
    }
    break ;
    
  default :
    cout << "Unkown mapping in set_sec_order_r2" << endl ;
    abort() ;
    break ;
  }
}

void Param_elliptic::set_sec_order (int zone, double a, double b, double c){
 
  if ((type_map == MAP_AFF) || (mp_log->get_type(zone) == AFFINE)) {
    cout << "set_sec_order only defined for a log mapping" << endl ;
    abort() ;
  }
  else {
 
    int nz = get_mp().get_mg()->get_nzone() ;
    
    int nr ;
    
    int conte = 0 ;
    for (int l=0 ; l<nz ; l++) {
      
      nr = get_mp().get_mg()->get_nr(l) ;
      
      for (int k=0 ; k<get_mp().get_mg()->get_np(l)+1 ; k++)
	for (int j=0 ; j<get_mp().get_mg()->get_nt(l) ; j++) {
	  if ((operateurs[conte] != 0x0) && (l==zone)) {

	    int old_base = operateurs[conte]->get_base_r() ;
	    // PROVISOIRE, DANS LE NOYAU SEUL LE CAS SPHERIQUE EST IMPLEMENTE
	    if (old_base != R_CHEBI) {
	      delete operateurs[conte] ;
	      operateurs[conte] = new Ope_sec_order (nr, old_base, 
						     get_alpha(l), 
						     get_beta(l), a, b, c) ;
	    }
	  }
	  conte ++ ;
	}
    }
  }
}

void Param_elliptic::set_ope_vorton (int zone, Scalar& source) {
  
  source.set_spectral_va().coef() ;
  source.set_spectral_va().ylm() ;

  int nz = get_mp().get_mg()->get_nzone() ;
  int nr ;
  int conte = 0 ;
  int m_quant, base_r_1d, l_quant ;
  int dzpuis = source.get_dzpuis() ;

  switch (type_map) {

  case MAP_AFF:
 
    for (int l=0 ; l<nz ; l++) {
      int dz = (l==nz-1) ? dzpuis : 0 ;
      nr = get_mp().get_mg()->get_nr(l) ;
      
      for (int k=0 ; k<get_mp().get_mg()->get_np(l)+1 ; k++)
	for (int j=0 ; j<get_mp().get_mg()->get_nt(l) ; j++) {
	  if (l==zone) {
	    if (operateurs[conte] != 0x0) {	    
	      delete operateurs[conte] ;
	      source.get_spectral_va().base.give_quant_numbers 
		(l, k, j, m_quant, l_quant, base_r_1d) ;
	      operateurs[conte] = new Ope_vorton (nr, base_r_1d, get_alpha(l), 
							   get_beta(l), l_quant, dz) ;
	    }
	  }
	  conte ++ ;
	}
    }
    break ;
    
  case MAP_LOG :
    if (mp_log->get_type(zone) != AFFINE) {
      cout << "Operator not define with LOG mapping..." << endl ;
      abort() ;
    }
    else {
      for (int l=0 ; l<nz ; l++) {
	int dz = (l==nz-1) ? dzpuis : 0 ;
	nr = get_mp().get_mg()->get_nr(l) ;
	
	for (int k=0 ; k<get_mp().get_mg()->get_np(l)+1 ; k++)
	  for (int j=0 ; j<get_mp().get_mg()->get_nt(l) ; j++) {
	    if (l==zone) {
	      if (operateurs[conte] != 0x0) {	    
		delete operateurs[conte] ;
		source.get_spectral_va().base.give_quant_numbers 
		  (l, k, j, m_quant, l_quant, base_r_1d) ;
		operateurs[conte] = new Ope_vorton (nr, base_r_1d,
							     get_alpha(l), get_beta(l), l_quant, dz) ;
	      }
	    }
	    conte ++ ;
	  }
      }
    }
    break ;
    
  default :
    cout << "Unkown mapping in set_ope_vorton" << endl ;
    abort() ;
    break ;
  }
}

void Param_elliptic::set_variable_F (const Scalar& so) {

  assert (so.get_etat() != ETATNONDEF) ;
  assert (so.get_mp() == get_mp()) ;
  
  var_F = so ;
  done_F.annule_hard() ;
}

void Param_elliptic::set_variable_G (const Scalar& so) {

  assert (so.get_etat() != ETATNONDEF) ;
  assert (so.get_mp() == get_mp()) ;
  
  var_G = so ;
  done_G.annule_hard() ;
}
