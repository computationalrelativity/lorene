/*
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *   Copyright (c) 2002 Jerome Novak
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


 

/*
 * $Id$
 * $Log$
 * Revision 1.9  2016/12/05 16:18:17  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.8  2014/10/13 08:53:43  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.7  2014/10/06 15:13:19  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.6  2003/03/03 19:41:34  f_limousin
 * Suppression of an assert on a metric associated with a tensor.
 *
 * Revision 1.5  2002/10/16 14:37:15  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.4  2002/09/10 13:44:17  j_novak
 * The method "manipule" of one indice has been removed for Tenseur_sym objects
 * (the result cannot be a Tenseur_sym).
 * The method "sans_trace" now computes the traceless part of a Tenseur (or
 * Tenseur_sym) of valence 2.
 *
 * Revision 1.3  2002/09/06 14:49:26  j_novak
 * Added method lie_derive for Tenseur and Tenseur_sym.
 * Corrected various errors for derive_cov and arithmetic.
 *
 * Revision 1.2  2002/08/07 16:14:12  j_novak
 * class Tenseur can now also handle tensor densities, this should be transparent to older codes
 *
 * Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
 * LORENE
 *
 * Revision 2.2  2000/02/09  19:32:22  eric
 * MODIF IMPORTANTE: la triade de decomposition est desormais passee en
 * argument des constructeurs.
 *
 * Revision 2.1  2000/01/11  11:15:08  eric
 * Gestion de la base vectorielle (triad).
 *
 * Revision 2.0  1999/12/02  17:19:02  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>
#include <cassert>
#include <cmath>

// Headers Lorene
#include "tenseur.h"
#include "metrique.h"

namespace Lorene {
Tenseur_sym operator*(const Tenseur& t1, const Tenseur_sym& t2) {
   
    assert ((t1.get_etat() != ETATNONDEF) && (t2.etat != ETATNONDEF)) ;
    assert (t1.get_mp() == t2.mp) ;
    
    int val_res = t1.get_valence() + t2.valence ;
    double poids_res = t1.get_poids() + t2.poids ;
    poids_res = (fabs(poids_res) < 1.e-10 ? 0. : poids_res) ;
    const Metrique* met_res = 0x0 ;
    if (poids_res != 0.) {
      // assert((t1.get_metric() != 0x0) || (t2.metric != 0x0)) ;
      if (t1.get_metric() != 0x0) met_res = t1.get_metric() ;
      else met_res = t2.metric ;
    }
   
    Itbl tipe (val_res) ;
    tipe.set_etat_qcq() ;
    for (int i=0 ; i<t1.get_valence() ; i++)
	tipe.set(i) = t1.get_type_indice(i) ;
    for (int i=0 ; i<t2.valence ; i++)
	tipe.set(i+t1.get_valence()) = t2.type_indice(i) ;
	

    if ( t1.get_valence() != 0 ) {
	    assert ( *(t1.get_triad()) == *(t2.get_triad()) ) ;
    }

    Tenseur_sym res(*t1.get_mp(), val_res, tipe, *(t2.get_triad()),
		    met_res, poids_res) ;
	

    if ((t1.get_etat() == ETATZERO) || (t2.etat == ETATZERO))
	res.set_etat_zero() ;
    else {
	res.set_etat_qcq() ;
	Itbl jeux_indice_t1 (t1.get_valence()) ;
	jeux_indice_t1.set_etat_qcq() ;
	Itbl jeux_indice_t2 (t2.valence) ;
	jeux_indice_t2.set_etat_qcq() ;
	    
	for (int i=0 ; i<res.n_comp ; i++) {
	    Itbl jeux_indice_res(res.donne_indices(i)) ;
	    for (int j=0 ; j<t1.get_valence() ; j++)
		jeux_indice_t1.set(j) = jeux_indice_res(j) ;
	    for (int j=0 ; j<t2.valence ; j++)
		jeux_indice_t2.set(j) = jeux_indice_res(j+t1.get_valence()) ;
		
	    res.set(jeux_indice_res) = t1(jeux_indice_t1)*t2(jeux_indice_t2) ;
	}
    }
    return res ;
}

Tenseur_sym manipule (const Tenseur_sym& t1, const Metrique& met) {
    
    Tenseur* auxi ;
    Tenseur* auxi_old = new Tenseur(t1) ;
    
    for (int i=0 ; i<t1.valence ; i++) {
	auxi = new Tenseur(manipule(*auxi_old, met, i)) ;
	delete auxi_old ;
	auxi_old = new Tenseur(*auxi) ;
	delete auxi ;
    }
    
    Tenseur_sym result(*auxi_old) ;
    delete auxi_old ;
    return result ;
}
  
Tenseur_sym lie_derive (const Tenseur_sym& t, const Tenseur& x, 
			const Metrique* met)
{
  assert ( (t.get_etat() != ETATNONDEF) && (x.get_etat() != ETATNONDEF) ) ;
  assert(x.get_valence() == 1) ;
  assert(x.get_type_indice(0) == CON) ;
  assert(x.get_poids() == 0.) ;
  assert(t.get_mp() == x.get_mp()) ;
 
  int val = t.get_valence() ;
  double poids = t.get_poids() ;
  Itbl tipe (val+1) ;
  tipe.set_etat_qcq() ;
  tipe.set(0) = COV ;
  Itbl tipx(2) ;
  tipx.set_etat_qcq() ;
  tipx.set(0) = COV ;
  tipx.set(1) = CON ;
  for (int i=0 ; i<val ; i++)
    tipe.set(i+1) = t.get_type_indice(i) ;
  Tenseur_sym dt(*t.get_mp(), val+1, tipe, *t.get_triad(), t.get_metric(), 
		 poids) ;
  Tenseur dx(*x.get_mp(), 2, tipx, x.get_triad()) ; 
  if (met == 0x0) {
    dx = x.gradient() ;
    dt = t.gradient() ;
  }
  else {
    dx = x.derive_cov(*met) ;
    dt = t.derive_cov(*met) ;
  }
  Tenseur_sym resu(contract(x,0,dt,0)) ;
  Tenseur* auxi ;
  if ((val!=0)&&(t.get_etat()!=ETATZERO)&&(x.get_etat()!=ETATZERO)) {
    assert(t.get_triad()->identify() == x.get_triad()->identify()) ;

    for (int i=0 ; i<val ; i++) {
      if (t.get_type_indice(i) == COV) {
	auxi = new Tenseur(contract(t,i,dx,1)) ;
	
	Itbl indices_aux(val) ;
	indices_aux.set_etat_qcq() ;
	for (int j=0 ; j<resu.get_n_comp() ; j++) {
	  
	  Itbl indices (resu.donne_indices(j)) ;
	  indices_aux.set(val-1) = indices(i) ;
	  for (int idx=0 ; idx<val-1 ; idx++)
	    if (idx<i)
	      indices_aux.set(idx) = indices(idx) ;
	    else
	      indices_aux.set(idx) = indices(idx+1) ;
	  
	  resu.set(indices) += (*auxi)(indices_aux) ;
	}
      }   
      else {
	auxi = new Tenseur(contract(t,i,dx,0)) ;
	
	Itbl indices_aux(val) ;
	indices_aux.set_etat_qcq() ;
	
	//On range comme il faut :
	for (int j=0 ; j<resu.get_n_comp() ; j++) {
	  
	  Itbl indices (resu.donne_indices(j)) ;
	  indices_aux.set(val-1) = indices(i) ;
	  for (int idx=0 ; idx<val-1 ; idx++)
	    if (idx<i)
	      indices_aux.set(idx) = indices(idx) ;
	    else
	      indices_aux.set(idx) = indices(idx+1) ;
	  resu.set(indices) -= (*auxi)(indices_aux) ;
	}
      }
      delete auxi ;
    }
  }
  if ((poids != 0.)&&(t.get_etat()!=ETATZERO)&&(x.get_etat()!=ETATZERO)) 
    resu = resu + poids*contract(dx,0,1)*t ;

  return resu ;
}

Tenseur_sym sans_trace(const Tenseur_sym& t, const Metrique& metre) 
{
  assert(t.get_etat() != ETATNONDEF) ;
  assert(metre.get_etat() != ETATNONDEF) ;
  assert(t.get_valence() == 2) ;

  Tenseur_sym resu(t) ;
  if (resu.get_etat() == ETATZERO) return resu ;
  assert(resu.get_etat() == ETATQCQ) ;

  int t0 = t.get_type_indice(0) ;
  int t1 = t.get_type_indice(1) ;
  Itbl mix(2) ;
  mix.set_etat_qcq() ;
  mix.set(0) = (t0 == t1 ? -t0 : t0) ;
  mix.set(1) = t1 ;

  Tenseur tmp(*t.get_mp(), 2, mix, *t.get_triad(), t.get_metric(), 
	      t.get_poids()) ;
  if (t0 == t1)
    tmp = manipule(t, metre, 0) ;
  else
    tmp = t ;

  Tenseur trace(contract(tmp, 0, 1)) ;

  if (t0 == t1) {
	switch (t0) {
	case COV : {
	  resu = resu - 1./3.*trace * metre.cov() ;
	  break ;
	}
	case CON : {
	  resu = resu - 1./3.*trace * metre.con() ;	
	  break ;
	}
	default :
	  cout << "Erreur bizarre dans sans_trace!" << endl ;
	  abort() ;
	  break ;
	}
  }
  else {
    Tenseur_sym delta(*t.get_mp(), 2, mix, *t.get_triad(), 
		      t.get_metric(), t.get_poids()) ;
    delta.set_etat_qcq() ;
    for (int i=0; i<3; i++) 
      for (int j=i; j<3; j++)
	delta.set(i,j) = (i==j ? 1 : 0) ;
    resu = resu - trace/3. * delta ;
  }
  resu.set_std_base() ;
  return resu ;
}
}
