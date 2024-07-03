/*
 *   Copyright (c) 2024 Jerome Novak
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
 * Methods of classes Map_af and Map_star to interpolate between one another.
 * (see files map.h and map_star.h) 
 */ 

#include "tensor.h"
#include "param.h"
#include "proto.h"
#include "utilitaires.h"

using namespace Lorene ;

    //---------------------------------------------------------------
    // Static variables : arrays of coefficient/summation functions
    //---------------------------------------------------------------
namespace
{
  bool first_call = true ;
  
  void (*coef_r[MAX_BASE])(const int*, const int*, double*, const int*, double*) ;
  double (*som_r_1d[MAX_BASE])(const double*, int, double) ;
}

// Local functions prototypes
//----------------------------
void do_first_call_initializations(
  void (*coef_r0[MAX_BASE])(const int*, const int*, double*, const int*, double*),
  double (*som_r_1d0[MAX_BASE])(const double*, int, double) ) ;
bool check_grids(const Map_af&, const Map_star&, int&) ;
bool check_grids(const Map_star&, const Map_af&, int&) ;
void pasprevu_r(const int*, const int*, double*, const int*, double*) ;
void base_non_def_r(const int*, const int*, double*, const int*, double*) ;



            //---------------------------------------------
            //  Spectral summation from Map_af to Map_star
            //---------------------------------------------

Scalar Map_star::interpolate_from_map_af(const Scalar& f_a) const
{
  
  //## Only valid for TYPE_T = SYM & TYPE_P = SYM for the moment (to be improved)
  assert((mg->get_type_t() == SYM) && (mg->get_type_p() == SYM) ) ; 

  // First call operations
  if (first_call) {
    do_first_call_initializations(coef_r, som_r_1d) ;
    first_call = false ;
  }

  // Check whether the input Scalar is defined on a Map_af
  const Map* p_mp = &f_a.get_mp() ;
  const Map_af* p_mpa = dynamic_cast<const Map_af*>(p_mp) ;
  if (p_mpa == 0x0) 
    throw(invalid_argument("Map_star::interpolate_from_map_af: the input Scalar field is not of type Map_af.")) ;

  const Map_af& mpa = *p_mpa ;
  const Mg3d* mg_a = mpa.get_mg() ;

  
  // Checks that both grids are compatible... (## add symmetries!!)
  int ndom_a = 0 ;
  if (!check_grids(mpa, *this, ndom_a)) {
    throw(invalid_argument("Map_star::interpolate_from_map_af: both mappings are not compatible")) ;
  }
  int np0 = mg_a->get_np(0) ;
  int nt0 = mg_a->get_nt(0) ;

  Scalar resu(*this) ;
  
  switch (f_a.get_etat()) {
    
  case ETATNONDEF: {
    throw(invalid_argument("Map_star::interpolate_from_map_af: the input Scalar field is not defined.")) ;
    break;
  }
    
  case ETATZERO:{
    resu.set_etat_zero() ;
    break ;
  }
    
  case ETATUN:{
    resu.set_etat_one() ;
    break ;
  }
    
  case ETATQCQ:{ // General case
    resu.allocate_all() ; 
    const Base_val& base = f_a.get_spectral_va().base ;

    //Intermediate array: coefficients in r, values in \theta & \varphi
    Mtbl interm(mg_a) ;
    interm.annule_hard() ;

    // Compute function values  # Improve with a transform in theta, phi only
    f_a.get_spectral_va().coef_i() ;

    // Makes the coefficient transform in the r variable only
    //-------------------------------------------------------
    for (int l=0; l<ndom_a; l++) {
      const Tbl* f = ((f_a.get_spectral_va().c)->t)[l] ;
      Tbl* cf = (interm.t)[l];
      if (f->get_etat() == ETATZERO) {
	cf->set_etat_zero() ;
	continue ; // nothing to do if the Tbl = 0  
      }	    
      int nr = f->get_dim(0) ;
      int deg[3] ;
      deg[0] = np0 ;
      deg[1] = nt0 ;
      deg[2] = nr ;
      *cf = *f ;
      
      // Takes information on the r-basis 
      int base_r = ( base.b[l] & MSQ_R ) >> TRA_R ;
      assert(base_r < MAX_BASE) ; 
      
      // Transformation in r:
      // --------------------
      if ( nr > 1 ) {
	assert( admissible_fft(nr-1) || (mg->get_colloc_r(l) != BASE_CHEB) );
	coef_r[base_r](deg, deg, (cf->t), deg, (cf->t)) ;
      }	
    } // end of loop on domains
    
    // Computing values on grid points associated to Map_star coordinates
    //-------------------------------------------------------------------
    int nzs = get_mg()->get_nzone() ;
    Param par_dummy ;
    
    for (int lz=0; lz<nzs; lz++) {
      for (int k=0; k<np0; k++) {
	for (int j=0; j<nt0; j++) {
	  int nrs = get_mg()->get_nr(lz) ;
	  for (int i=0; i<nrs; i++) {
	    double radius = (+r)(lz, k, j, i) ;
	    int l_a = 0 ;
	    double xi_a = 0. ;
	    mpa.val_lx_jk(radius, j, k, par_dummy, l_a, xi_a) ;
	    int base_r = ( base.b[l_a] & MSQ_R ) >> TRA_R ;
	    const double* cfr = &interm.set(l_a, k, j, 0) ;
	    int nra = mg_a->get_nr(l_a) ;

	    // Call to the 'som_r_1d' function that makes the spectral summation
	    resu.set_grid_point(lz, k, j, i) = som_r_1d[base_r](cfr, nra, xi_a) ;
	  }
	}
      }
    }
    break ;
  } // end of case ETATQCQ
  default: {
   throw(invalid_argument("Map_star::interpolate_from_map_af: the input Scalar field is ill-formed.")) ;
   break ;
  }
  } // end of switch
    
  return resu ;

}

            //---------------------------------------------
            //  Spectral summation from Map_star to Map_af
            //---------------------------------------------

Scalar Map_af::interpolate_from_map_star(const Scalar& f_s) const
{
  
  //## Only valid for TYPE_T = SYM & TYPE_P = SYM for the moment (to be improved)
  assert((mg->get_type_t() == SYM) && (mg->get_type_p() == SYM) ) ; 

  // First call operations
  if (first_call) {
    do_first_call_initializations(coef_r, som_r_1d) ;
    first_call = false ;
  }

  // Check whether the input Scalar is defined on a Map_star
  const Map* p_mp = &f_s.get_mp() ;
  const Map_star* p_mps = dynamic_cast<const Map_star*>(p_mp) ;
  if (p_mps == 0x0) 
    throw(invalid_argument("Map_af::interpolate_from_map_star: the input Scalar field is not of type Map_star.")) ;

  const Map_star& mps = *p_mps ;
  const Mg3d* mg_s = mps.get_mg() ;
  
  // Checks that both grids are compatible...
  int ndom_a = 0 ;
  if (!check_grids(mps, *this, ndom_a)) {
    throw(invalid_argument("Map_star::interpolate_from_map_star: both mappings are not compatible")) ;
  }
  int np0 = mg_s->get_np(0) ;
  int nt0 = mg_s->get_nt(0) ;

  Scalar resu(*this) ;
  
  switch (f_s.get_etat()) {
    
  case ETATNONDEF: {
    throw(invalid_argument("Map_af::interpolate_from_map_star: the input Scalar field is not defined.")) ;
    break;
  }
    
  case ETATZERO:{
    resu.set_etat_zero() ;
    break ;
  }
    
  case ETATUN:{
    resu.set_etat_one() ;
    break ;
  }
    
  case ETATQCQ:{ // General case
    resu.allocate_all() ; 
    const Base_val& base = f_s.get_spectral_va().base ;

    //Intermediate array: coefficients in r, values in \theta & \varphi
    Mtbl interm(mg_s) ;
    interm.annule_hard() ;

    // Compute function values  # Improve with a transform in theta, phi only
    f_s.get_spectral_va().coef_i() ;

    // Makes the coefficient transform in the r variable only
    //-------------------------------------------------------
    int nzs = mg_s->get_nzone() ;
    for (int l=0; l<nzs; l++) {
      const Tbl* f = ((f_s.get_spectral_va().c)->t)[l] ;
      Tbl* cf = (interm.t)[l];
      if (f->get_etat() == ETATZERO) {
	cf->set_etat_zero() ;
	continue ; // nothing to do if the Tbl = 0  
      }	    
      int nr = f->get_dim(0) ;
      int deg[3] ;
      deg[0] = np0 ;
      deg[1] = nt0 ;
      deg[2] = nr ;
      *cf = *f ;
      
      // Takes information on the r-basis 
      int base_r = ( base.b[l] & MSQ_R ) >> TRA_R ;
      assert(base_r < MAX_BASE) ; 
      
      // Transformation in r:
      // --------------------
      if ( nr > 1 ) {
	assert( admissible_fft(nr-1) || (mg->get_colloc_r(l) != BASE_CHEB) );
	coef_r[base_r](deg, deg, (cf->t), deg, (cf->t)) ;
      }	
    } // end of loop on domains
    
    // Computing values on grid points associated to Map_star coordinates
    //-------------------------------------------------------------------
    int nza = get_mg()->get_nzone() ;
    Param par_dummy ;
    const Coord& r_s = mps.r ;
    
    int nrs_max = mg_s->get_nr(nzs-1) ;
    for (int lz=0; lz<ndom_a; lz++) {
      for (int k=0; k<np0; k++) {
	for (int j=0; j<nt0; j++) {
	  double rmax_s = (+r_s)(nzs-1, k, j, nrs_max-1) ; 
	  int nra = get_mg()->get_nr(lz) ;
	  for (int i=0; i<nra; i++) {
	    double radius_a = (+r)(lz, k, j, i) ;
	    if (radius_a <= rmax_s) {
	      int l_s = 0 ;
	      double xi_s = 0. ;
	      mps.val_lx_jk(radius_a, j, k, par_dummy, l_s, xi_s) ;
	      int base_r = ( base.b[l_s] & MSQ_R ) >> TRA_R ;
	      const double* cfr = &interm.set(l_s, k, j, 0) ;
	      int nrs = mg_s->get_nr(l_s) ;

	    // Call to the 'som_r_1d' function that makes the spectral summation
	    resu.set_grid_point(lz, k, j, i) = som_r_1d[base_r](cfr, nrs, xi_s) ;
	    }
	    else // the point is outside the Map_star grid
	      resu.set_grid_point(lz, k, j, i) = 0. ;
	  }
	}
      }
    }
    resu.annule(ndom_a, nza-1) ;// Zero outside the domains containing the Map_star 
    break ;
  } // end of case ETATQCQ
  default: {
   throw(invalid_argument("Map_af::interpolate_from_map_star: the input Scalar field is ill-formed.")) ;
   break ;
  }
  } // end of switch
    
  return resu ;

}

                       //------------------------
                       //     Miscellaneous
                       //------------------------

void do_first_call_initializations(
  void (*coef_r0[MAX_BASE])
  (const int*, const int*, double*, const int*, double*),
  double (*som_r_1d0[MAX_BASE])(const double*, int, double) ) {

  for (int i=0; i<MAX_BASE; i++) {
      coef_r0[i] = pasprevu_r ;
      som_r_1d0[i] = som_r_1d_pas_prevu ;
  }
  coef_r0[NONDEF] = base_non_def_r ;
  coef_r0[R_CHEB >> TRA_R] = cfrcheb ;	    
  coef_r0[R_CHEBU >> TRA_R] = cfrcheb ;	    
  coef_r0[R_CHEBP >> TRA_R] = cfrchebp ;	    
  coef_r0[R_CHEBI >> TRA_R] = cfrchebi ;	    
  coef_r0[R_LEG >> TRA_R] = cfrleg ;	    
  coef_r0[R_LEGP >> TRA_R] = cfrlegp ;	    
  coef_r0[R_LEGI >> TRA_R] = cfrlegi ;	    
  coef_r0[R_JACO02 >> TRA_R] = cfrjaco02 ;

  som_r_1d0[R_CHEB >> TRA_R] = som_r_1d_cheb ;
  som_r_1d0[R_CHEBP >> TRA_R] = som_r_1d_chebp ;
  som_r_1d0[R_CHEBI >> TRA_R] = som_r_1d_chebi ;
  som_r_1d0[R_CHEBU >> TRA_R] = som_r_1d_cheb ; // no error here... 
  som_r_1d0[R_LEG >> TRA_R] = som_r_1d_leg ;
  som_r_1d0[R_LEGP >> TRA_R] = som_r_1d_legp ;
  som_r_1d0[R_LEGI >> TRA_R] = som_r_1d_legi ;
  som_r_1d0[R_JACO02 >> TRA_R] = som_r_1d_jaco02 ;
}

bool check_grids(const Map_af& mpa, const Map_star& mps, int& ndom) {

  const Mg3d* mga = mpa.get_mg() ;
  const Mg3d* mgs = mps.get_mg() ;
  int nza = mga->get_nzone() ;
  int nzs = mgs->get_nzone() ;

  // First, check that both grids have same symmetries.
  bool resu = ( mga->get_type_t() == mgs->get_type_t() )
    && ( mga->get_type_p() == mgs->get_type_p() ) ;

  // Then, check that the number of points in theta & phi are the same.
  int np0 = mga->get_np(0) ;
  int nt0 = mga->get_nt(0) ;
  for (int l=1; l<nza; l++) 
    resu = resu && (np0 == mga->get_np(l)) && (nt0 == mga->get_nt(l)) ;
  for (int l=0; l<nzs; l++) 
    resu = resu && (np0 == mgs->get_np(l)) && (nt0 == mgs->get_nt(l)) ;

  // Finally, looks at the maximal radius of Map_star ...
  const Coord& rr = mps.r ;
  double rmax = 0. ;
  int jmax = -1 ;
  int kmax = -1 ;
  int nr = mgs->get_nr(nzs-1) ;
  for (int k=0; k<np0; k++) {
    for (int j=0; j<nt0; j++) {
      if ((+rr)(nzs-1, k, j, nr-1) > rmax) {
	rmax = (+rr)(nzs-1, k, j, nr-1) ;
	jmax = j ;
	kmax = k ;
      }
    }
  }
  // ... and the corresponding domain for Map_af.
  double dummy = 0. ;
  Param par_dummy ;
  mpa.val_lx_jk(rmax, jmax, kmax, par_dummy, ndom, dummy) ;
  ndom++ ;
  
  return resu ;
}

bool check_grids(const Map_star& mps, const Map_af& mpa, int& ndom) {

  const Mg3d* mga = mpa.get_mg() ;
  const Mg3d* mgs = mps.get_mg() ;
  int nza = mga->get_nzone() ;
  int nzs = mgs->get_nzone() ;

  // First, check that both grids have same symmetries.
  bool resu = ( mga->get_type_t() == mgs->get_type_t() )
    && ( mga->get_type_p() == mgs->get_type_p() ) ;

  // Then, check that the number of points in theta & phi are the same.
  int np0 = mga->get_np(0) ;
  int nt0 = mga->get_nt(0) ;
  for (int l=1; l<nza; l++) 
    resu = resu && (np0 == mga->get_np(l)) && (nt0 == mga->get_nt(l)) ;
  for (int l=0; l<nzs; l++) 
    resu = resu && (np0 == mgs->get_np(l)) && (nt0 == mgs->get_nt(l)) ;

  // Finally, looks at the maximal radius of Map_star ...
  const Coord& rr = mps.r ;
  double rmax = 0. ;
  int jmax = -1 ;
  int kmax = -1 ;
  int nr = mgs->get_nr(nzs-1) ;
  for (int k=0; k<np0; k++) {
    for (int j=0; j<nt0; j++) {
      if ((+rr)(nzs-1, k, j, nr-1) > rmax) {
	rmax = (+rr)(nzs-1, k, j, nr-1) ;
	jmax = j ;
	kmax = k ;
      }
    }
  }
  // ... and the corresponding domain for Map_af.
  double dummy = 0. ;
  Param par_dummy ;
  mpa.val_lx_jk(rmax, jmax, kmax, par_dummy, ndom, dummy) ;
  ndom++ ;
  
  return resu ;
}

void pasprevu_r(const int*, const int*, double*, const int*, double*) {
    cout << "Valeur::coef: the required expansion basis in r " << endl ;
    cout << "  is not implemented !" << endl ;
    abort() ;
}

void base_non_def_r(const int*, const int*, double*, const int*, double*) {
    cout << "Valeur::coef: the expansion basis in r is undefined !" << endl ;
    abort() ;
}


