/*
 *  Methods of class Map_log
 *
 *   (see file map.h for documentation)
 *
 */

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


char map_log_C[] = "$Header $" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/06/22 08:49:58  p_grandclement
 * Addition of everything needed for using the logarithmic mapping
 *
 * 
 * $Header$
 *
 */

// headers C
#include <math.h>

// headers Lorene
#include "itbl.h"
#include "tbl.h"
#include "coord.h"
#include "grilles.h"
#include "map.h"

			
			//---------------//
			// Constructeurs //
			//---------------//

// Constructor from a grid
// -----------------------
Map_log::Map_log (const Mg3d& mgrille, const Tbl& bornes, const Itbl& typevar) : 
  Map_radial(mgrille), alpha (mgrille.get_nzone()), beta (mgrille.get_nzone()), 
  type_var(typevar)
										
{
    // Les bornes
    int nzone = mg->get_nzone() ;
    
    alpha.set_etat_qcq() ;
    beta.set_etat_qcq() ;

    for (int l=0 ; l<nzone ; l++) {
      switch (type_var(l)) {
      case AFFINE: {
	switch (mg->get_type_r(l)) {
	case RARE:	{
	  alpha.set(l) = bornes(l+1) - bornes(l) ;
	  beta.set(l) = bornes(l) ;
	  break ; 
	}
	    
	case FIN:	{
	  alpha.set(l) = (bornes(l+1) - bornes(l)) * .5 ;
	  beta.set(l) = (bornes(l+1) + bornes(l)) * .5 ;
	  break ;
	}
	  
	case UNSURR: {
	  double umax = 1./bornes(l) ;
	  double umin = 1./bornes(l+1) ;
	  alpha.set(l) = (umin - umax) * .5 ;  // u est une fonction decroissante
	  beta.set(l) = (umin + umax) * .5 ;   //  de l'indice i en r
	  break ;
	}
	  
	default:	{
	  cout << "Map_log::Map_log: unkown type_r ! " << endl ;
	  abort () ;
	  break ;
	}
	    
	}
	break ;
      }
      case LOG:{
	switch (mg->get_type_r(l)) {
	case FIN:	{
	  alpha.set(l) = (log(bornes(l+1)) - log(bornes(l))) * .5 ;
	  beta.set(l) = (log(bornes(l+1)) + log(bornes(l))) * .5 ;
	  break ;
	}
	  
	default:	{
	  cout << "Map_log::Map_log: unkown type_r ! " << endl ;
	  abort () ;
	  break ;
	}
	}
	break ;
      }
      
      default: { 
	cout << "Map_log::Map_log: unkown type_r ! " << endl ;
	abort () ;
	break ;
      }
      }
    }

    set_coord() ;
}

Map_log::Map_log (const Map_log& so) : Map_radial (*so.mg), alpha(so.alpha), beta(so.beta), 
				       type_var(so.type_var) {
  set_coord() ;
}


Map_log::Map_log (const Mg3d& mgrille, FILE* fd) : Map_radial (mgrille, fd), 
						   alpha (mgrille.get_nzone()), 
						   beta (mgrille.get_nzone()), 
						   type_var (mgrille.get_nzone()) {
  Tbl alpha_lu (fd) ;
  Tbl beta_lu (fd) ;
  Itbl type_lu (fd) ;

  alpha = alpha_lu ;
  beta = beta_lu ;
  type_var = type_lu ;

  set_coord() ;
}

Map_log::~Map_log() {}

// Sauvegarde :
void Map_log::sauve(FILE* fd) const {

  Map_radial::sauve(fd) ;
  alpha.sauve(fd) ;
  beta.sauve(fd) ;
  type_var.sauve(fd) ;
}

// Comparison operator :
bool Map_log::operator==(const Map& mpi) const {
  
  // Precision of the comparison
  double precis = 1e-10 ;
  bool resu = true ;

  // Dynamic cast pour etre sur meme Map...
  const Map_log* mp0 = dynamic_cast<const Map_log*>(&mpi) ;
  if (mp0 == 0x0)
    resu = false ;
  else {
    if (*mg != *(mpi.get_mg()))
      resu = false ;
    
    if (fabs(ori_x-mpi.get_ori_x()) > precis) resu = false ;
    if (fabs(ori_y-mpi.get_ori_y()) > precis) resu = false ;
    if (fabs(ori_z-mpi.get_ori_z()) > precis)  resu = false ;

    if (bvect_spher != mpi.get_bvect_spher()) resu = false ;
    if (bvect_cart != mpi.get_bvect_cart()) resu = false ;

    if (diffrelmax (alpha, mp0->alpha) > precis) resu = false ;
    if (diffrelmax (beta, mp0->beta) > precis) resu = false ;
    if (diffrelmax(type_var, mp0->type_var) > precis) resu = false ;
  }

  return resu ;
}
			//------------//
			// Impression //
			//------------//

ostream & Map_log::operator>>(ostream & ost) const {

  ost << "Log mapping (class Map_log)" << endl ; 
  int nz = mg->get_nzone() ;
  for (int l=0; l<nz; l++) {
    ost << "     Domain #" << l << " ; Variable type " ;
    if (type_var(l) == AFFINE)
      ost << "affine : " ;
    if (type_var(l) == LOG)
      ost << "log : " ;
    ost << "alpha_l = " << alpha(l) 
	 << " , beta_l = " << beta(l) << endl ;  
  }

  
  ost << "            Coord r : " ; 
  for (int l=0; l<nz; l++) {
    int nrm1 = mg->get_nr(l) - 1 ;
    ost << " " << (+r)(l, 0, 0, nrm1) ; 
  }
    ost << endl ; 

    return ost ;
}    

// Affectation to a affine mapping :
void Map_log::operator=(const Map_af& mpi) {
    
  assert(mpi.get_mg() == mg) ; 
    
  set_ori( mpi.get_ori_x(), mpi.get_ori_y(), mpi.get_ori_z() ) ; 
    
  set_rot_phi( mpi.get_rot_phi() ) ; 
    
    // The members bvect_spher and bvect_cart are treated by the functions
    //  set_ori and set_rot_phi.
    
    for (int l=0; l<mg->get_nzone(); l++){
	alpha.set(l) = mpi.get_alpha()[l] ; 
	beta.set(l) = mpi.get_beta()[l] ; 
    }
    
    type_var = AFFINE ;

    reset_coord() ;
    set_coord() ;
}


            //-------------------------------------------------//
	    //  Assignement of the Coord building functions    //
	    //-------------------------------------------------//
	    
void Map_log::set_coord(){

    // ... Coord's introduced by the base class Map : 
    r.set(this, map_log_fait_r) ;
    tet.set(this, map_log_fait_tet) ;
    phi.set(this, map_log_fait_phi) ;
    sint.set(this, map_log_fait_sint) ;
    cost.set(this, map_log_fait_cost) ;
    sinp.set(this, map_log_fait_sinp) ;
    cosp.set(this, map_log_fait_cosp) ;

    x.set(this, map_log_fait_x) ;
    y.set(this, map_log_fait_y) ;
    z.set(this, map_log_fait_z) ;

    xa.set(this, map_log_fait_xa) ;
    ya.set(this, map_log_fait_ya) ;
    za.set(this, map_log_fait_za) ;
    
    // ... Coord's introduced by the base class Map_radial : 
    xsr.set(this, map_log_fait_xsr) ;
    dxdr.set(this, map_log_fait_dxdr) ;
    drdt.set(this, map_log_fait_drdt) ;
    stdrdp.set(this, map_log_fait_stdrdp) ;
    srdrdt.set(this, map_log_fait_srdrdt) ;
    srstdrdp.set(this, map_log_fait_srstdrdp) ;
    sr2drdt.set(this, map_log_fait_sr2drdt) ;
    sr2stdrdp.set(this, map_log_fait_sr2stdrdp) ;
    d2rdx2.set(this, map_log_fait_d2rdx2) ;
    lapr_tp.set(this, map_log_fait_lapr_tp) ;
    d2rdtdx.set(this, map_log_fait_d2rdtdx) ;
    sstd2rdpdx.set(this, map_log_fait_sstd2rdpdx) ;
    sr2d2rdt2.set(this, map_log_fait_sr2d2rdt2) ;
    
    // ... Coord's introduced by the base Map_log itself
    dxdlnr.set(this, map_log_fait_dxdlnr) ;
}
