/*
 *   Copyright (c) 2004 Philippe Grandclement
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


char valeur_val_propre_1d_C[] = "$Header$" ;



/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/06 15:13:24  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2004/08/24 09:14:52  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * 
 * $Header$
 *
 */

// headers du C
#include <cassert>
#include <cstdlib>

// headers Lorene

#include "type_parite.h"
#include "valeur.h"
#include "matrice.h"
#include "proto.h"


//****************************************************************
//                           CAS PAIR
//****************************************************************

void rotate_propre_pair (Valeur& so, bool sens) {
  
  so.coef() ;
  so.set_etat_cf_qcq() ;

  static int nt_courant = 0 ;
  static Matrice* passage = 0x0 ;
  static Matrice* inv = 0x0 ;
  
  int nt = so.mg->get_nt(0) ;
 
  
  if (nt != nt_courant) {
    // On doit calculer les matrices... Pas de la tarte...
    if (passage != 0x0) 
      delete passage ;
    if (inv != 0x0)
      delete inv ;

    nt_courant = nt ;
    
    Matrice ope (nt, nt) ;
    ope.set_etat_qcq() ;
    for (int i=0 ; i<nt ; i++)
      for (int j=0 ; j<nt ; j++)
	ope.set(i,j) = 0 ;
    
    double c_courant = 0 ;
    for (int j=0 ; j<nt ; j++) {
      ope.set(0, j) = 2*j ;
      for (int i=1 ; i<j ; i++)
	ope.set(i,j) = 4*j ;
      ope.set(j,j) = c_courant ;
      c_courant -= 8*j + 2 ;
    }

    passage = new Matrice (ope.vect_propre()) ;
    passage->set_band(nt-1,0) ;
    passage->set_lu() ;
    // Un peu nul pour calculer l'inverse mais bon...
    inv = new Matrice (nt, nt) ;
    inv->set_etat_qcq() ;
      
    Tbl auxi (nt) ;
    auxi.set_etat_qcq() ;

    for (int i=0 ; i<nt ; i++) {
      for (int j=0 ; j<nt ; j++)
	auxi.set(j) = 0 ;
      auxi.set(i) = 1 ;
      Tbl sortie (passage->inverse(auxi)) ;
      for (int j=0 ; j<nt ; j++)
	inv->set(j,i) = sortie(j) ;
    }    
  }

  // Fin du calcul des matrices...
  for (int l=0 ; l<so.mg->get_nzone() ; l++)
    if (so.c_cf->t[l]->get_etat() != ETATZERO)
      for (int k=0 ; k<so.mg->get_np(l) ; k++)
	for (int i=0 ; i<so.mg->get_nr(l) ; i++) {
	  
	  Tbl coefs(nt) ;
	  coefs.set_etat_qcq() ;
	  for (int j=0 ; j<nt ; j++)
	    coefs.set(j) = (*so.c_cf)(l,k,j,i) ;
	  Tbl prod (nt) ;
	  prod.set_etat_qcq() ;
	  for (int j=0 ; j<nt ; j++)
	    prod.set(j) = 0 ;
	  
	  if (sens) {
	    //calcul direct
	    for (int j=0 ; j<nt ; j++)
	      for (int jb=0 ; jb<nt ; jb++)
		prod.set(j) += (*inv)(j, jb) * coefs(jb) ;
	  }
	  else {
	    //calcul inverse
	    for (int j=0 ; j<nt ; j++)
	      for (int jb=0 ; jb<nt ; jb++)
		prod.set(j) += (*passage)(j, jb) * coefs(jb) ;
	  }
	  
	  for (int j=0 ; j<nt ; j++)
	    so.c_cf->set(l,k,j,i) = prod(j) ;
	  
	}
  
  // On met les bonnes bases :
  int base_tet = so.base.get_base_t(0) ;
  Base_val new_base (so.base) ;
  
  if (sens) {
    switch (base_tet) {
    case T_COS_P:
      new_base.set_base_t(T_CL_COS_P) ;
      break ;
    case T_SIN_P:
      new_base.set_base_t(T_CL_SIN_P) ;
      break ;
    default:
      cout << "Problem in rotate_propre_pair" << endl ;
      abort() ;
      break ;
    }
  }
  else {
    switch (base_tet) {
    case T_CL_COS_P:
      new_base.set_base_t(T_COS_P) ;
      break ;
    case T_CL_SIN_P:
      new_base.set_base_t(T_SIN_P) ;
      break ;
    default:
      cout << "Problem in rotate_propre_pair" << endl ;
      abort() ;
      break ;
    }
  }

  so.c_cf->base = new_base ;
  so.base = new_base ; 
}

//****************************************************************
//                           CAS IMPAIR
//****************************************************************

void rotate_propre_impair (Valeur& so, bool sens) {  
  so.coef() ;
  so.set_etat_cf_qcq() ;

  static int nt_courant = 0 ;
  static Matrice* passage = 0x0 ;
  static Matrice* inv = 0x0 ;
  
  int nt = so.mg->get_nt(0) ;
 
  
  if (nt != nt_courant) {
    // On doit calculer les matrices... Pas de la tarte...
    if (passage != 0x0) 
      delete passage ;
    if (inv != 0x0)
      delete inv ;

    nt_courant = nt ;
    
    Matrice ope (nt, nt) ;
    ope.set_etat_qcq() ;
    for (int i=0 ; i<nt ; i++)
      for (int j=0 ; j<nt ; j++)
	ope.set(i,j) = 0 ;

    double c_courant = 0 ; 
    for (int j=0 ; j<nt ; j++) {
      for (int i=0 ; i<j ; i++)
	ope.set(i,j) = 2+4*j ;
      ope.set(j,j) = c_courant ;
      c_courant -= 8*j + 6 ;
    }

    passage = new Matrice (ope.vect_propre()) ;
    passage->set_band(nt-1,0) ;
    passage->set_lu() ;
    // Un peu nul pour calculer l'inverse mais bon...
    inv = new Matrice (nt, nt) ;
    inv->set_etat_qcq() ;
      
    Tbl auxi (nt) ;
    auxi.set_etat_qcq() ;

    for (int i=0 ; i<nt ; i++) {
      for (int j=0 ; j<nt ; j++)
	auxi.set(j) = 0 ;
      auxi.set(i) = 1 ;
      Tbl sortie (passage->inverse(auxi)) ;
      for (int j=0 ; j<nt ; j++)
	inv->set(j,i) = sortie(j) ;
    }
  }   
  
  // Fin du calcul des matrices...

  for (int l=0 ; l<so.mg->get_nzone() ; l++)
    if (so.c_cf->t[l]->get_etat() != ETATZERO)
      for (int k=0 ; k<so.mg->get_np(l) ; k++)
	for (int i=0 ; i<so.mg->get_nr(l) ; i++) {
	  
	  Tbl coefs(nt) ;
	  coefs.set_etat_qcq() ;
	  for (int j=0 ; j<nt ; j++)
	    coefs.set(j) = (*so.c_cf)(l,k,j,i) ;
	  Tbl prod (nt) ;
	  prod.set_etat_qcq() ;
	  for (int j=0 ; j<nt ; j++)
	    prod.set(j) = 0 ;
	  
	  if (sens) {
	    //calcul direct
	    for (int j=0 ; j<nt ; j++)
	      for (int jb=0 ; jb<nt ; jb++)
		prod.set(j) += (*inv)(j, jb) * coefs(jb) ;
	  }
	  else {
	    //calcul inverse
	    for (int j=0 ; j<nt ; j++)
	      for (int jb=0 ; jb<nt ; jb++)
		prod.set(j) += (*passage)(j, jb) * coefs(jb) ;
	  }
	  
	  for (int j=0 ; j<nt ; j++)
	    so.c_cf->set(l,k,j,i) = prod(j) ;
	  
	}
  
  // On met les bonnes bases :
  int base_tet = so.base.get_base_t(0) ;
  Base_val new_base (so.base) ;
  
  if (sens) {
    switch (base_tet) {
    case T_COS_I:
      new_base.set_base_t(T_CL_COS_I) ;
      break ;
    case T_SIN_I:
      new_base.set_base_t(T_CL_SIN_I) ;
      break ;
    default:
      cout << "Problem in rotate_propre_impair" << endl ;
      abort() ;
      break ;
    }
  }
  else {
    switch (base_tet) {
    case T_CL_COS_I:
      new_base.set_base_t(T_COS_I) ;
      break ;
    case T_CL_SIN_I:
      new_base.set_base_t(T_SIN_I) ;
      break ;
    default:
      cout << "Problem in rotate_propre_impair" << endl ;
      abort() ;
      break ;
    }
  }

  so.c_cf->base = new_base ;
  so.base = new_base ;
}


void Valeur::val_propre_1d() {

  // On recupere la base en tet :
  int nz = mg->get_nzone() ;
  int base_tet = base.get_base_t(0) ;
  // On verifie que c'est bien la meme partout...
  for (int i=1 ; i<nz ; i++)
    assert (base_tet == base.get_base_t(i)) ;

  switch (base_tet) {
  case T_COS_P:
    rotate_propre_pair (*this, true) ;
    break ;
  case T_SIN_P:
    rotate_propre_pair (*this, true) ;
    break ;
  case T_COS_I:
    rotate_propre_impair (*this, true) ;
    break ;
  case T_SIN_I:
    rotate_propre_impair (*this, true) ;
    break ;
  case T_CL_COS_P:
    break ;
  case T_CL_SIN_P:
    break ;
  case T_CL_COS_I:
    break ;
  case T_CL_SIN_I:
    break ;
  default:
    cout << "Unknown basis in Valeur::val_propre_1d" << endl ;
    abort() ;
    break ;
  }
}
void Valeur::val_propre_1d_i() {

  // On recupere la base en tet :
  int nz = mg->get_nzone() ;
  int base_tet = base.get_base_t(0) ;
  // On verifie que c'est bien la meme partout...
  for (int i=1 ; i<nz ; i++)
    assert (base_tet == base.get_base_t(i)) ;

  switch (base_tet) {
  case T_CL_COS_P:
    rotate_propre_pair (*this, false) ;
    break ;
  case T_CL_SIN_P:
    rotate_propre_pair (*this, false) ;
    break ;
  case T_CL_COS_I:
    rotate_propre_impair (*this, false) ;
    break ;
  case T_CL_SIN_I:
    rotate_propre_impair (*this, false) ;
    break ; 
  case T_COS_P:
    break ;
  case T_SIN_P:
    break ;
  case T_COS_I:
    break ;
  case T_SIN_I:
    break ;
  default:
    cout << "Unknown basis in Valeur::val_propre_1d_i" << endl ;
    abort() ;
    break ;
  }
}


