/*
 *   Copyright (c) 2003 Philippe Grandclement
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

#include "ope_elementary.h"
#include "param_elliptic.h"
#include "change_var.h"
#include "base_val.h" 
#include "scalar.h"

// Construction (By default it is set to Poisson with appropriate dzpuis...)
Param_elliptic::Param_elliptic(const Scalar& so) {

  // On passe en Ylm
  Scalar auxi(so) ;
  auxi.set_spectral_va().coef() ;
  auxi.set_spectral_va().ylm() ;

  Base_val base (auxi.get_spectral_va().base) ;
  int dzpuis = auxi.get_dzpuis() ;
  
  // Right now, only applicable with affine mapping
  const Map_af* map_affine = dynamic_cast <const Map_af*> (&so.get_mp()) ;
  mp = map_affine ;

  if (map_affine == 0x0) {
    cout << "Param_elliptic only defined for affine mapping" << endl ;
    abort() ;
  }
  else {
  
    int nz = mp->get_mg()->get_nzone() ;
    int nbr = 0 ;
    for (int l=0 ; l<nz ; l++)
      nbr += (mp->get_mg()->get_np(l)+1)*
	mp->get_mg()->get_nt(l) ;
    
    operateurs = new Ope_elementary* [nbr] ;
    variables = new Change_var* [nbr] ;
    
    for (int l=0 ; l<nbr ; l++) {
      operateurs[l] = 0x0 ;
      variables[l] = 0x0 ;
    }
    
    int nr ;
    double alpha, beta ;
    int base_r, m_quant, l_quant ;
    
    int conte = 0 ;
    for (int l=0 ; l<nz ; l++) {
      
      nr = mp->get_mg()->get_nr(l) ;
      
      alpha = mp->get_alpha()[l] ;
      beta = mp->get_beta()[l] ;
      
      for (int k=0 ; k<mp->get_mg()->get_np(l)+1 ; k++)
	for (int j=0 ; j<mp->get_mg()->get_nt(l) ; j++) {
	  if (nullite_plm(j, mp->get_mg()->get_nt(l), k, 
			  mp->get_mg()->get_np(l), base)) {
	    
	    donne_lm(nz, l, j, k, base, m_quant, l_quant, base_r) ;
	    operateurs[conte] = new 
	      Ope_poisson(nr, base_r, alpha, beta, l_quant, dzpuis) ;
	    variables[conte] = new Change_var(STD) ;
	  }
	  conte ++ ;
	}
    }
  }
}

Param_elliptic::~Param_elliptic () {
  
  int nbr = 0 ;
  for (int l=0 ; l<mp->get_mg()->get_nzone() ; l++) 
    nbr += (mp->get_mg()->get_np(l)+1)*mp->get_mg()->get_nt(l) ;

  for (int i=0 ; i<nbr ; i++) {
    if (operateurs[i] != 0x0)
      delete operateurs[i] ;
    if (variables[i] != 0x0)
      delete variables[i] ;
  }

  delete [] operateurs ;
  delete [] variables ;
}

void Param_elliptic::inc_l_quant (int zone) {
  
  int np, nt, nr ;

  int conte = 0 ;
  for (int l=0 ; l<mp->get_mg()->get_nzone() ; l++) {
    
    np = mp->get_mg()->get_np(l) ;
    nt = mp->get_mg()->get_nt(l) ;
    nr = mp->get_mg()->get_nr(l) ;

    for (int k=0 ; k<np+1 ; k++)
      for (int j=0 ; j<nt ; j++) {
	if ((operateurs[conte] != 0x0) && (l==zone))
	  operateurs[conte]->inc_l_quant() ;
	conte ++ ;
      }
  }
}


void Param_elliptic::set_helmholtz_minus (int zone, double masse) {
 
  int nz = mp->get_mg()->get_nzone() ;

  int nr ;
  double alpha, beta ;

  int conte = 0 ;
  for (int l=0 ; l<nz ; l++) {

    nr = mp->get_mg()->get_nr(l) ;

    alpha = mp->get_alpha()[l] ;
    beta = mp->get_beta()[l] ;
  
    for (int k=0 ; k<mp->get_mg()->get_np(l)+1 ; k++)
      for (int j=0 ; j<mp->get_mg()->get_nt(l) ; j++) {
	if ((operateurs[conte] != 0x0) && (l==zone)) {
	  int old_base = operateurs[conte]->get_base_r() ;
	  // PROVISOIRE, DANS LE NOYAU SEUL LE CAS SPHERIQUE EST IMPLEMENTE
	  if (old_base != R_CHEBI) {
	    delete operateurs[conte] ;
	    operateurs[conte] = new Ope_helmholtz_minus (nr, old_base, 
							 alpha, beta, masse) ;
	  }
	}
	conte ++ ;
      }
  }
}

void Param_elliptic::set_helmholtz_plus (int zone, double masse) {
 
  int nz = mp->get_mg()->get_nzone() ;
  
  int nr ;
  double alpha, beta ;

  int conte = 0 ;
  for (int l=0 ; l<nz ; l++) {

    nr = mp->get_mg()->get_nr(l) ;

    alpha = mp->get_alpha()[l] ;
    beta = mp->get_beta()[l] ;
  
    for (int k=0 ; k<mp->get_mg()->get_np(l)+1 ; k++)
      for (int j=0 ; j<mp->get_mg()->get_nt(l) ; j++) {
	if ((operateurs[conte] != 0x0) && (l==zone)) {
	  int old_base = operateurs[conte]->get_base_r() ;
	  // PROVISOIRE, DANS LE NOYAU SEUL LE CAS SPHERIQUE EST IMPLEMENTE
	  if (old_base != R_CHEBI) {
	    delete operateurs[conte] ;
	    operateurs[conte] = new Ope_helmholtz_plus (nr, old_base, 
							alpha, beta,  masse) ;
	  }
	}
	conte ++ ;
      }
  }
}

void Param_elliptic::set_poisson_vect_r(int zone) {
 
  int nz = mp->get_mg()->get_nzone() ;
  
  int nr ;
  double alpha, beta ;
  
  int conte = 0 ;
  for (int l=0 ; l<nz ; l++) {
    
    nr = mp->get_mg()->get_nr(l) ;
    
    alpha = mp->get_alpha()[l] ;
    beta = mp->get_beta()[l] ;
    
    for (int k=0 ; k<mp->get_mg()->get_np(l)+1 ; k++)
      for (int j=0 ; j<mp->get_mg()->get_nt(l) ; j++) {
	if ((operateurs[conte] != 0x0) && (l==zone)) {
	  int old_base = operateurs[conte]->get_base_r() ;
	  Ope_poisson* op_pois = 
	    dynamic_cast<Ope_poisson*>(operateurs[conte]) ;
	  assert (op_pois !=0x0) ;
	  int lq_old = op_pois->get_lquant() ;
	  int dzp = op_pois->get_dzpuis() ;

	  delete operateurs[conte] ;
	  operateurs[conte] = new Ope_pois_vect_r(nr, old_base,alpha, 
						  beta, lq_old, dzp) ;
	}
      }
    conte ++ ;
  }
}


void Param_elliptic::set_sec_order_r2 (int zone, double a, double b, double c){
 
  int nz = mp->get_mg()->get_nzone() ;
  
  int nr ;
  double alpha, beta ;

  int conte = 0 ;
  for (int l=0 ; l<nz ; l++) {

    nr = mp->get_mg()->get_nr(l) ;

    alpha = mp->get_alpha()[l] ;
    beta = mp->get_beta()[l] ;
  
    for (int k=0 ; k<mp->get_mg()->get_np(l)+1 ; k++)
      for (int j=0 ; j<mp->get_mg()->get_nt(l) ; j++) {
	if ((operateurs[conte] != 0x0) && (l==zone)) {
	  int old_base = operateurs[conte]->get_base_r() ;
	  // PROVISOIRE, DANS LE NOYAU SEUL LE CAS SPHERIQUE EST IMPLEMENTE
	  if (old_base != R_CHEBI) {
	    delete operateurs[conte] ;
	    operateurs[conte] = new Ope_sec_order_r2 (nr, old_base, 
						      alpha, beta, a, b, c) ;
	  }
	}
	conte ++ ;
      }
  }
}

void Param_elliptic::set_variable (int type_variable, int zone) {
 
  int np, nt, nr ;
  
  int conte = 0 ;
  for (int l=0 ; l<mp->get_mg()->get_nzone() ; l++) {
    
    np = mp->get_mg()->get_np(l) ;
    nt = mp->get_mg()->get_nt(l) ;
    nr = mp->get_mg()->get_nr(l) ;
    
    for (int k=0 ; k<np+1 ; k++)
      for (int j=0 ; j<nt ; j++) {
	if ((variables[conte] != 0x0) && (l==zone) && (k==0) && (j==0)) {
	  delete variables[conte] ;
	  variables[conte] = new Change_var(type_variable) ;
	}
	conte ++ ;
      }
  }
}

void Param_elliptic::set_variable (int type_variable, double mult, int zone) {
 
  int np, nt, nr ;
  
  int conte = 0 ;
  for (int l=0 ; l<mp->get_mg()->get_nzone() ; l++) {
    
    np = mp->get_mg()->get_np(l) ;
    nt = mp->get_mg()->get_nt(l) ;
    nr = mp->get_mg()->get_nr(l) ;
    
    for (int k=0 ; k<np+1 ; k++)
      for (int j=0 ; j<nt ; j++) {
	if ((variables[conte] != 0x0) && (l==zone) && (k==0) && (j==0)) {
	  delete variables[conte] ;
	  variables[conte] = new Change_var(type_variable, mult) ;
	}
	conte ++ ;
      }
  }
}

void Param_elliptic::set_variable (int type_variable, double mult, double add, int zone) {
 
  int np, nt, nr ;
  
  int conte = 0 ;
  for (int l=0 ; l<mp->get_mg()->get_nzone() ; l++) {
    
    np = mp->get_mg()->get_np(l) ;
    nt = mp->get_mg()->get_nt(l) ;
    nr = mp->get_mg()->get_nr(l) ;
    
    for (int k=0 ; k<np+1 ; k++)
      for (int j=0 ; j<nt ; j++) {
	if ((variables[conte] != 0x0) && (l==zone) && (k==0) && (j==0)) {
	  delete variables[conte] ;
	  variables[conte] = new Change_var(type_variable, mult, add) ;
	}
	conte ++ ;
      }
  }
}
