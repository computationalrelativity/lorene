/*
 * Methods to compare the values of tensors for a Bin_star and a Bin_ns_ncp
 *
 * (see file etoile.h for documentation)
 */

/*
 *   Copyright (c) 2003 Francois Limousin
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


char bin_ns_compare_C[] = "$Header$" ;

/*
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "bin_ns_ncp.h"
#include "eos.h"
#include "binaire.h"
#include "utilitaires.h"

void Bin_ns_ncp::compare(const Binaire& binaire) {

  // Parameters of the grid
  //-----------------------
  
  for (int n=1; n<3; n++) {

  double nzone = (*((binaire(n)).get_mp().get_mg())).get_nzone() ;
  double nr = (*((binaire(n)).get_mp().get_mg())).get_nr(0) ;
  double nt = (*((binaire(n)).get_mp().get_mg())).get_nt(0) ;
  double np = (*((binaire(n)).get_mp().get_mg())).get_np(0) ;


  // Declarations 
  //-------------

  // lapse

  double diff_logn_auto, diffp_logn_auto, i_logn_auto, j_logn_auto ;
  double k_logn_auto, l_logn_auto, moy_logn_auto;
  
  double diff_logn_comp, diffp_logn_comp, i_logn_comp, j_logn_comp ;
  double k_logn_comp, l_logn_comp, moy_logn_comp ;

  double diff_nnn, diffp_nnn, i_nnn, j_nnn ;
  double k_nnn, l_nnn, moy_nnn ;

  // shift

  double diff_shift_auto[3], diffp_shift_auto[3], i_shift_auto[3] ;
  double j_shift_auto[3], k_shift_auto[3], l_shift_auto[3], moy_shift_auto[3] ;

  double diff_shift_comp[3], diffp_shift_comp[3], i_shift_comp[3] ;
  double j_shift_comp[3], k_shift_comp[3], l_shift_comp[3], moy_shift_comp[3] ;

  double diff_shift[3], diffp_shift[3], i_shift[3], j_shift[3] ;
  double k_shift[3], l_shift[3], moy_shift[3] ;

  // a_car

  double diff_a_car, diffp_a_car, i_a_car, j_a_car ;
  double k_a_car, l_a_car, moy_a_car ; 

  // beta

  double diff_beta, diffp_beta, i_beta, j_beta;
  double k_beta, l_beta, moy_beta ;

  // hydro quantities

  double diff_ener, diffp_ener, i_ener, j_ener ;
  double k_ener, l_ener, moy_ener ;

  double diff_press, diffp_press, i_press, j_press ;
  double k_press, l_press, moy_press ;

  double diff_ener_euler, diffp_ener_euler, i_ener_euler, j_ener_euler ;
  double k_ener_euler, l_ener_euler, moy_ener_euler ;

  double diff_s_euler, diffp_s_euler, i_s_euler, j_s_euler ;
  double k_s_euler, l_s_euler, moy_s_euler ;

  double diff_gam_euler, diffp_gam_euler, i_gam_euler, j_gam_euler ;
  double k_gam_euler, l_gam_euler, moy_gam_euler ;

  double diff_u_euler[3], diffp_u_euler[3], i_u_euler[3], j_u_euler[3] ;
  double k_u_euler[3], l_u_euler[3], moy_u_euler[3] ;

  // miscallaneous

  double diff_tkij_auto[9], diffp_tkij_auto[9], i_tkij_auto[9] ;
  double j_tkij_auto[9], k_tkij_auto[9], l_tkij_auto[9], moy_tkij_auto[9] ;

  double diff_tkij_comp[9], diffp_tkij_comp[9], i_tkij_comp[9] ;
  double j_tkij_comp[9], k_tkij_comp[9], l_tkij_comp[9], moy_tkij_comp[9] ;

  double diff_kcar_auto, diffp_kcar_auto, i_kcar_auto, j_kcar_auto ;
  double k_kcar_auto, l_kcar_auto, moy_kcar_auto ;

  double diff_kcar_comp, diffp_kcar_comp, i_kcar_comp, j_kcar_comp ;
  double k_kcar_comp, l_kcar_comp, moy_kcar_comp ;

  double diff_mass_b, diff_mass_g, diff_ray_eq, diff_ray_pole ;
  double diff_omega, diff_x_axe ;

  int compt_logn_auto, compt_logn_comp, compt_nnn, compt_shift_auto[3] ;
  int compt_shift_comp[3], compt_shift[3], compt_a_car, compt_beta,compt_ener ;
  int compt_press, compt_ener_euler, compt_s_euler, compt_gam_euler ; 
  int compt_u_euler[3], compt_tkij_auto[9], compt_tkij_comp[9] ;
  int compt_kcar_auto, compt_kcar_comp ;

  // Quantities providing from star.
  //-------------------------------

  Tenseur logn_auto_star = binaire(n).get_logn_auto() ;
  Tenseur logn_comp_star = binaire(n).get_logn_comp() ;
  Tenseur nnn_star = binaire(n).get_nnn() ;

  Tenseur shift_auto_star = binaire(n).get_shift_auto() ;
  Tenseur shift_comp_star = binaire(n).get_shift_comp() ;
  Tenseur shift_star = binaire(n).get_shift();
  
  shift_auto_star.change_triad(((*et[n-1]).get_mp()).get_bvect_cart()) ;
  shift_comp_star.change_triad(((*et[n-1]).get_mp()).get_bvect_cart()) ;
  shift_star.change_triad(((*et[n-1]).get_mp()).get_bvect_cart()) ;


  Tenseur a_car_star = binaire(n).get_a_car() ;

  Tenseur beta_auto_star = binaire(n).get_beta_auto() ;
  Tenseur beta_comp_star = binaire(n).get_beta_comp() ;

  Tenseur ener_star = binaire(n).get_ener() ;
  Tenseur press_star = binaire(n).get_press() ;
  Tenseur ener_euler_star = binaire(n).get_ener_euler() ;
  Tenseur s_euler_star = binaire(n).get_s_euler() ;
  Tenseur gam_euler_star = binaire(n).get_gam_euler() ;

  Tenseur u_euler_star = binaire(n).get_u_euler() ;
  u_euler_star.change_triad(((*et[n-1]).get_mp()).get_bvect_cart()) ;
  
  Tenseur tkij_auto_star = binaire(n).get_tkij_auto() ;
  Tenseur tkij_comp_star = binaire(n).get_tkij_comp() ;

  tkij_auto_star.change_triad(((*et[n-1]).get_mp()).get_bvect_cart()) ;
  tkij_comp_star.change_triad(((*et[n-1]).get_mp()).get_bvect_cart()) ;

  Tenseur kcar_auto_star = binaire(n).get_akcar_auto() ;
  Tenseur kcar_comp_star = binaire(n).get_akcar_comp() ;

  double mass_b_star = binaire(n).mass_b() ;
  double mass_g_star = binaire(n).mass_g() ;
  double ray_eq_star = binaire(n).ray_eq() ;
  double ray_pole_star = binaire(n).ray_pole() ;
  double omega_star = binaire.get_omega() ;
  double x_axe_star = binaire.get_x_axe() ;

  // Initialisations
  //----------------

  diff_logn_auto = -1. ;
  moy_logn_auto = 0. ;
  compt_logn_auto = 0 ;
  diff_logn_comp = -1. ;
  moy_logn_comp = 0. ;
  compt_logn_comp = 0 ;
  diff_nnn = -1. ;
  moy_nnn = 0. ;
  compt_nnn = 0 ;

  for (int m=0; m<3; m++) {
    diff_shift_auto[m] = -1. ;
    moy_shift_auto[m] = 0. ;
    compt_shift_auto[m] = 0 ;
    diff_shift_comp[m] = -1. ;
    moy_shift_comp[m] = 0. ;
    compt_shift_comp[m] = 0 ;
    diff_shift[m] = -1. ;
    moy_shift[m] = 0. ;
    compt_shift[m] = 0 ;
  }

  diff_a_car = -1. ;
  moy_a_car = 0. ;
  compt_a_car = 0 ;

  diff_beta = -1. ;
  moy_beta = 0. ;
  compt_beta = 0 ;

  diff_ener = -1. ;
  moy_ener = 0. ;
  compt_ener = 0 ;
  diff_press = -1. ;
  moy_press = 0. ;
  compt_press = 0 ;
  diff_ener_euler = -1. ;
  moy_ener_euler = 0. ;
  compt_ener_euler = 0 ;
  diff_s_euler = -1. ;
  moy_s_euler = 0. ;
  compt_s_euler = 0 ;
  diff_gam_euler = -1. ;
  moy_gam_euler = 0. ;
  compt_gam_euler = 0 ;
  for (int m=0; m<3; m++) {
    diff_u_euler[m] = -1. ;
    moy_u_euler[m] = 0. ;
    compt_u_euler[m] = 0 ;
  }

  for (int m=0; m<9; m++) { 
    diff_tkij_auto[m] = -1. ;
    moy_tkij_auto[m] = 0. ;
    compt_tkij_auto[m] = 0 ;
    diff_tkij_comp[m] = -1. ;
    moy_tkij_comp[m] = 0. ;
    compt_tkij_comp[m] = 0 ;
  }
  
  diff_kcar_auto = -1 ;
  moy_kcar_auto = 0 ; 
  compt_kcar_auto = 0 ;
  diff_kcar_comp = -1 ;
  moy_kcar_comp = 0 ; 
  compt_kcar_comp = 0 ;
 
  // Begin of the comparison 
  //------------------------


  for (int i=0; i<nr; i++) {

    for (int j=0; j<nt; j++) {
      for (int k=0; k<np; k++) {
	for (int l=0; l<nzone; l++) {
	
	  // Lapse
	  //-------

	  // logn_auto
	  
	  if (fabs((*et[n-1]).get_logn_auto()()(l,k,j,i)) > 1.e-10) {
	    diffp_logn_auto = fabs(((*et[n-1]).get_logn_auto()()(l,k,j,i) 
	  -logn_auto_star()(l,k,j,i))/(*et[n-1]).get_logn_auto()()(l,k,j,i)) ;  
	    moy_logn_auto += diffp_logn_auto ; 
	    compt_logn_auto++ ;

	    if (diffp_logn_auto > diff_logn_auto) {
	      diff_logn_auto = diffp_logn_auto ;
	      l_logn_auto = l ;
	      k_logn_auto= k ;
	      j_logn_auto = j ;
	      i_logn_auto = i ;
	    }
	  }   

	  // logn_comp

	  if (fabs((*et[n-1]).get_logn_comp()()(l,k,j,i)) > 1.e-10) {
	    diffp_logn_comp = fabs(((*et[n-1]).get_logn_comp()()(l,k,j,i) 
  	   - logn_comp_star()(l,k,j,i))/(*et[n-1]).get_logn_comp()()(l,k,j,i)) ;
	    moy_logn_comp += diffp_logn_comp ; 
	    compt_logn_comp++ ;

	    if (diffp_logn_comp > diff_logn_comp) {
	      diff_logn_comp = diffp_logn_comp ;
	      l_logn_comp = l ;
	      k_logn_comp = k ;
	      j_logn_comp = j ;
	      i_logn_comp = i ;
	    }
	  } 

	  // nnn

	  if (fabs((*et[n-1]).get_nnn()()(l,k,j,i)) > 1.e-10) {
	    diffp_nnn = fabs(((*et[n-1]).get_nnn()()(l,k,j,i) 
	       	      - nnn_star()(l,k,j,i))/(*et[n-1]).get_nnn()()(l,k,j,i)) ;
	    moy_nnn += diffp_nnn ; 
	    compt_nnn++ ;
	      
	    if (diffp_nnn > diff_nnn) {
	      diff_nnn = diffp_nnn ;
	      l_nnn = l ;
	      k_nnn = k ;
	      j_nnn = j ;
	      i_nnn = i ;
	    }
	  }

	  // Shift
	  //------

	  // shift_auto 

	      for (int m=0; m<3; m++) {
	    if (fabs((*et[n-1]).get_shift_auto()(m)(l,k,j,i)) > 1.e-10) {
	    diffp_shift_auto[m] = fabs(((*et[n-1]).get_shift_auto()(m)(l,k,j,i) 
       - shift_auto_star(m)(l,k,j,i))/(*et[n-1]).get_shift_auto()(m)(l,k,j,i)) ;
	      moy_shift_auto[m] += diffp_shift_auto[m] ;
	      compt_shift_auto[m]++ ;

	      if (diffp_shift_auto[m] > diff_shift_auto[m]) {
		diff_shift_auto[m] = diffp_shift_auto[m] ;
		l_shift_auto[m] = l ;
		k_shift_auto[m] = k ;
		j_shift_auto[m] = j ;
		i_shift_auto[m] = i ;
	      }
	    }
	  }

	  // shift_comp 

	  for (int m=0; m<3; m++) {
	    if (fabs((*et[n-1]).get_shift_comp()(m)(l,k,j,i)) > 1.e-10) {
	    diffp_shift_comp[m] = fabs(((*et[n-1]).get_shift_comp()(m)(l,k,j,i) 
       - shift_comp_star(m)(l,k,j,i))/(*et[n-1]).get_shift_comp()(m)(l,k,j,i)) ;
	      moy_shift_comp[m] += diffp_shift_comp[m] ;
	      compt_shift_comp[m]++ ;

	      if (diffp_shift_comp[m] > diff_shift_comp[m]) {
		diff_shift_comp[m] = diffp_shift_comp[m] ;
		l_shift_comp[m] = l ;
		k_shift_comp[m] = k ;
		j_shift_comp[m] = j ;
		i_shift_comp[m] = i ;
	      }
	    }
	  }

	  // shift


	  for (int m=0; m<3; m++) {
	    if (fabs((*et[n-1]).get_shift()(m)(l,k,j,i)) > 1.e-10) {
	      diffp_shift[m] = fabs(((*et[n-1]).get_shift()(m)(l,k,j,i) 
	       - shift_star(m)(l,k,j,i))/(*et[n-1]).get_shift()(m)(l,k,j,i)) ;
	      moy_shift[m] += diffp_shift[m] ;
	      compt_shift[m]++ ;

	      if (diffp_shift[m] > diff_shift[m]) {
		diff_shift[m] = diffp_shift[m] ;
		l_shift[m] = l ;
		k_shift[m] = k ;
		j_shift[m] = j ;
		i_shift[m] = i ;
	      }
	    }
	  } 


	  // a_car
	  //-------

	  if (fabs((*et[n-1]).get_a_car()()(l,k,j,i)) > 1.e-10) {
	    diffp_a_car = fabs(((*et[n-1]).get_a_car()()(l,k,j,i) 
      		- a_car_star()(l,k,j,i))/(*et[n-1]).get_a_car()()(l,k,j,i)) ;
	    moy_a_car += diffp_a_car ; 
	    compt_a_car++ ;

	    if (diffp_a_car > diff_a_car) {
	      diff_a_car = diffp_a_car ;
	      l_a_car = l ;
	      k_a_car = k ;
	      j_a_car = j ;
	      i_a_car = i ;
	    }
	  }


	  // beta

	  Cmp beta = (log((*et[n-1]).get_a_car()()))/2. 
	    + log((*et[n-1]).get_nnn()()) ; 
	  if (fabs(beta(l,k,j,i)) > 1.e-10) {
	    diffp_beta = fabs((beta(l,k,j,i) - beta_auto_star()(l,k,j,i) 
			      - beta_comp_star()(l,k,j,i))/(beta(l,k,j,i))) ;
	    moy_beta += diffp_beta ; 
	    compt_beta++ ;
	      
	    if (diffp_beta > diff_beta) {
	      diff_beta = diffp_beta ;
	      l_beta = l ;
	      k_beta = k ;
	      j_beta = j ;
	      i_beta = i ;
	    }
	  }


	  // Hydro quantities
	  //-----------------

	  // ener
	  
	  if (fabs((*et[n-1]).get_ener()()(l,k,j,i)) > 1.e-10) {
	    diffp_ener = fabs(((*et[n-1]).get_ener()()(l,k,j,i) 
      	       - ener_star()(l,k,j,i))/(*et[n-1]).get_ener()()(l,k,j,i)) ;
	    moy_ener += diffp_ener ; 
	    compt_ener++ ;

	    if (diffp_ener > diff_ener) {
	      diff_ener = diffp_ener ;
	      l_ener = l ;
	      k_ener = k ;
	      j_ener = j ;
	      i_ener = i ;
	    }
	  }

	  // press

	  if (fabs((*et[n-1]).get_press()()(l,k,j,i)) > 1.e-10) {
	    diffp_press = fabs(((*et[n-1]).get_press()()(l,k,j,i) 
       	     	- press_star()(l,k,j,i))/(*et[n-1]).get_press()()(l,k,j,i)) ;
	    moy_press += diffp_press ; 
	    compt_press++ ;

	    if (diffp_press > diff_press) {
	      diff_press = diffp_press ;
	      l_press = l ;
	      k_press= k ;
	      j_press = j ;
	      i_press = i ;
	    }
	  }

	  // ener_euler

	  if (fabs((*et[n-1]).get_ener_euler()()(l,k,j,i)) > 1.e-10) {
	    diffp_ener_euler = fabs(((*et[n-1]).get_ener_euler()()(l,k,j,i) 
         - ener_euler_star()(l,k,j,i))/(*et[n-1]).get_ener_euler()()(l,k,j,i)) ;
	    moy_ener_euler += diffp_ener_euler ; 
	    compt_ener_euler++ ;

	    if (diffp_ener_euler > diff_ener_euler) {
	      diff_ener_euler = diffp_ener_euler;
	      l_ener_euler = l ;
	      k_ener_euler = k ;
	      j_ener_euler = j ;
	      i_ener_euler = i ;
	    }
	  }

	  // s_euler

	  if (fabs((*et[n-1]).get_s_euler()()(l,k,j,i)) > 1.e-10) {
	    diffp_s_euler = fabs(((*et[n-1]).get_s_euler()()(l,k,j,i) 
       	  - s_euler_star()(l,k,j,i))/(*et[n-1]).get_s_euler()()(l,k,j,i)) ;
	    moy_s_euler += diffp_s_euler ; 
	    compt_s_euler++ ;

	    if (diffp_s_euler > diff_s_euler) {
	      diff_s_euler = diffp_s_euler;
	      l_s_euler = l ;
	      k_s_euler = k ;
	      j_s_euler = j ;
	      i_s_euler = i ;
	    }
	  }

	  // gam_euler

	  if (fabs((*et[n-1]).get_gam_euler()()(l,k,j,i)) > 1.e-10) {
	    diffp_gam_euler = fabs(((*et[n-1]).get_gam_euler()()(l,k,j,i) 
	   - gam_euler_star()(l,k,j,i))/(*et[n-1]).get_gam_euler()()(l,k,j,i)) ;
	    moy_gam_euler += diffp_gam_euler ; 
	    compt_gam_euler++ ;

	    if (diffp_gam_euler > diff_gam_euler) {
	      diff_gam_euler = diffp_gam_euler ;
	      l_gam_euler = l ;
	      k_gam_euler = k ;
	      j_gam_euler = j ;
	      i_gam_euler = i ;
	    }
	  }

	  // u_euler

	  for (int m=0; m<3; m++) {
	    if (fabs((*et[n-1]).get_u_euler()(m)(l,k,j,i)) > 1.e-10) {
	      diffp_u_euler[m] = fabs(((*et[n-1]).get_u_euler()(m)(l,k,j,i) 
             - u_euler_star(m)(l,k,j,i))/(*et[n-1]).get_u_euler()(m)(l,k,j,i)) ;
	      moy_u_euler[m] += diffp_u_euler[m] ;
	      compt_u_euler[m]++ ;

	      if (diffp_u_euler[m] > diff_u_euler[m]) {
		diff_u_euler[m] = diffp_u_euler[m] ;
		l_u_euler[m] = l ;
		k_u_euler[m] = k ;
		j_u_euler[m] = j ;
		i_u_euler[m] = i ;
	      }
	    }
	  } 



	  // Miscallaneous
	  //---------------

	  // tkij_auto

	  for (int m=0; m<3; m++) {
	    for (int p=0; p<3; p++) {
	      if (fabs((*et[n-1]).get_tkij_auto()(m,p)(l,k,j,i)) > 1.e-10) {
		diffp_tkij_auto[3*m+p] = fabs(((*et[n-1]).get_tkij_auto()(m,p)(l,k,j,i) 
    - tkij_auto_star(m,p)(l,k,j,i))/(*et[n-1]).get_tkij_auto()(m,p)(l,k,j,i)) ;
	      moy_tkij_auto[3*m+p] += diffp_tkij_auto[3*m+p] ;
	      compt_tkij_auto[3*m+p]++ ;

		if (diffp_tkij_auto[3*m+p] > diff_tkij_auto[3*m+p]) {
		  diff_tkij_auto[3*m+p] = diffp_tkij_auto[3*m+p] ;
		  l_tkij_auto[3*m+p] = l ;
		  k_tkij_auto[3*m+p] = k ;
		  j_tkij_auto[3*m+p] = j ;
		  i_tkij_auto[3*m+p] = i ;
		}
	      }
	    } 
	  }

	  // tkij_comp

	  for (int m=0; m<3; m++) {
	    for (int p=0; p<3; p++) {
	      if (fabs((*et[n-1]).get_tkij_comp()(m,p)(l,k,j,i)) > 1.e-10) {
		diffp_tkij_comp[3*m+p] = fabs(((*et[n-1]).get_tkij_comp()(m,p)(l,k,j,i) 
     - tkij_comp_star(m,p)(l,k,j,i))/(*et[n-1]).get_tkij_comp()(m,p)(l,k,j,i)) ;
	      moy_tkij_comp[3*m+p] += diffp_tkij_comp[3*m+p] ;
	      compt_tkij_comp[3*m+p]++ ;

		if (diffp_tkij_comp[3*m+p] > diff_tkij_comp[3*m+p]) {
		  diff_tkij_comp[3*m+p] = diffp_tkij_comp[3*m+p] ;
		  l_tkij_comp[3*m+p] = l ;
		  k_tkij_comp[3*m+p] = k ;
		  j_tkij_comp[3*m+p] = j ;
		  i_tkij_comp[3*m+p] = i ;
		}
	      }
	    } 
	  }
	 
 
	  // kcar_auto

	  if (fabs((*et[n-1]).get_kcar_auto()()(l,k,j,i)) > 1.e-10) {
	    diffp_kcar_auto = fabs(((*et[n-1]).get_kcar_auto()()(l,k,j,i) 
	   - kcar_auto_star()(l,k,j,i))/(*et[n-1]).get_kcar_auto()()(l,k,j,i)) ;
	    moy_kcar_auto += diffp_kcar_auto ; 
	    compt_kcar_auto++ ;

	    if (diffp_kcar_auto > diff_kcar_auto) {
	      diff_kcar_auto = diffp_kcar_auto ;
	      l_kcar_auto = l ;
	      k_kcar_auto = k ;
	      j_kcar_auto = j ;
	      i_kcar_auto = i ;
	    }
	  }

	  // kcar_comp
	  
	  if (fabs((*et[n-1]).get_kcar_comp()()(l,k,j,i)) > 1.e-10) {
	    diffp_kcar_comp = fabs(((*et[n-1]).get_kcar_comp()()(l,k,j,i) 
	   - kcar_comp_star()(l,k,j,i))/(*et[n-1]).get_kcar_comp()()(l,k,j,i)) ;
	    moy_kcar_comp += diffp_kcar_comp ; 
	    compt_kcar_comp++ ;

	    if (diffp_kcar_comp > diff_kcar_comp) {
	      diff_kcar_comp = diffp_kcar_comp ;
	      l_kcar_comp = l ;
	      k_kcar_comp = k ;
	      j_kcar_comp = j ;
	      i_kcar_comp = i ;
	    }
	  }
	  


	}
      }
    }
  }

  moy_logn_auto /= compt_logn_auto ;
  moy_logn_comp /= compt_logn_comp ;
  moy_nnn /= compt_nnn ;
  
  for (int m=0; m<3; m++) {
    moy_shift_auto[m] /= compt_shift_auto[m] ;
    moy_shift_comp[m] /= compt_shift_comp[m] ;
    moy_shift[m] /= compt_shift[m] ;
    moy_u_euler[m] /= compt_u_euler[m] ;
  }
  
  moy_a_car /= compt_a_car ;
  moy_beta /= compt_beta ;
  moy_ener /= compt_ener ;
  moy_press /= compt_press ;
  moy_ener_euler /= compt_ener_euler ;
  moy_s_euler /= compt_s_euler ;
  moy_gam_euler /= compt_gam_euler ; 
  
  for (int m=0; m<9; m++) {
    moy_tkij_auto[m] /= compt_tkij_auto[m] ;
    moy_tkij_comp[m] /= compt_tkij_comp[m] ;
  }   

  moy_kcar_auto /= compt_kcar_auto ;
  moy_kcar_comp /= compt_kcar_comp ;


  diff_mass_b = fabs( ((*et[n-1]).mass_b() - mass_b_star) 
		      / (*et[n-1]).mass_b() ) ;
  diff_mass_g = fabs( ((*et[n-1]).mass_g() - mass_g_star) 
		      / (*et[n-1]).mass_g() ) ;
  diff_ray_eq = fabs( ((*et[n-1]).ray_eq() - ray_eq_star) 
		      / (*et[n-1]).ray_eq() ) ;
  diff_ray_pole = fabs( ((*et[n-1]).ray_pole() - ray_pole_star) 
			/ (*et[n-1]).ray_pole() ) ;
  diff_omega = fabs( (omega - omega_star)/omega ) ;
  if (x_axe == 0) {
    diff_x_axe = -1 ; 
  }
  else {
  diff_x_axe = fabs( (x_axe - x_axe_star)/x_axe ) ;
  }

	    
  // Outputs
  //-------------

  cout << endl << "---------------------------------------------------" 
       << endl ;

 
  
  // lapse
  //-------

  // logn_auto

  if (diff_logn_auto == -1) {
    cout << "logn_auto = 0" << endl ; 
  }
  else {
    cout << "max relative difference for logn_auto(" << l_logn_auto << "," 
	 << k_logn_auto << "," << j_logn_auto<< "," << i_logn_auto << ") = "
	 << diff_logn_auto << endl ;
    cout << "average difference for logn_auto : " << moy_logn_auto << endl ;
  }

  // logn_comp

  if (diff_logn_comp == -1) {
    cout << "logn_comp = 0" << endl ; 
  }
  else {
    cout << "max relative difference for logn_comp(" << l_logn_comp << "," 
	 << k_logn_comp << "," << j_logn_comp << "," << i_logn_comp << ") = "
	 << diff_logn_comp << endl ;
    cout << "average difference for logn_comp : " << moy_logn_comp << endl ;
  }
  
  // nnn

  if (diff_nnn == -1) {
    cout << "nnn = 0" << endl << endl ; 
  }
  else {
    cout << "max relative difference for nnn(" << l_nnn << "," 
	 << k_nnn << "," << j_nnn << "," << i_nnn << ") = "
	 << diff_nnn << endl  ;
    cout << "average difference for nnn : " << moy_nnn << endl << endl ;
  }
 
  // Shift
  //-------

  // shift_auto
  
  if (diff_shift_auto[0] == -1) {
    cout << "shift_auto_x = 0" << endl ; 
  }
  else {
    cout << "max relative difference for shift_auto_x(" << l_shift_auto[0] 
	 << "," << k_shift_auto[0] << "," << j_shift_auto[0] << "," 
	 << i_shift_auto[0] << ") = " << diff_shift_auto[0] << endl ;
     cout << "average difference for shift_auto_x : " << moy_shift_auto[0]
	  << endl ;
  }
  
  if (diff_shift_auto[1] == -1) {
    cout << "shift_auto_y = 0" << endl ; 
  }
  else {
    cout << "max relative difference for shift_auto_y(" << l_shift_auto[1] 
	 << "," << k_shift_auto[1] << "," << j_shift_auto[1] << "," 
	 << i_shift_auto[1] << ") = " << diff_shift_auto[1] << endl ;
     cout << "average difference for shift_auto_y : " << moy_shift_auto[1]
	  << endl ; 
 }
  
  if (diff_shift_auto[2] == -1) {
    cout << "shift_auto_z = 0" << endl << endl ; 
  }
  else {
    cout << "max relative difference for shift_auto_z(" << l_shift_auto[2] 
	 << "," << k_shift_auto[2] << "," << j_shift_auto[2] << "," 
	 << i_shift_auto[2] << ") = " << diff_shift_auto[2] << endl ;
      cout << "average difference for shift_auto_z : " << moy_shift_auto[2]
	  << endl << endl ; 
}
  
  // shift_comp

  if (diff_shift_comp[0] == -1) {
    cout << "shift_comp_x = 0" << endl ; 
  }
  else {
    cout << "max relative difference for shift_comp_x(" << l_shift[0] 
	 << "," << k_shift_comp[0] << "," << j_shift_comp[0] << "," 
	 << i_shift_comp[0] << ") = " << diff_shift_comp[0] << endl ;
      cout << "average difference for shift_comp_x : " << moy_shift_comp[0]
	  << endl ; 
}
  
  if (diff_shift_comp[1] == -1) {
    cout << "shift_comp_y = 0" << endl ; 
  }
  else {
    cout << "max relative difference for shift_comp_y(" << l_shift_comp[1] 
	 << "," << k_shift_comp[1] << "," << j_shift_comp[1] << "," 
	 << i_shift_comp[1] << ") = " << diff_shift_comp[1] << endl ;
      cout << "average difference for shift_comp_y : " << moy_shift_comp[1]
	  << endl ;
 }
  
  if (diff_shift_comp[2] == -1) {
    cout << "shift_comp_z = 0" << endl << endl ; 
  }
  else {
    cout << "max relative difference for shift_comp_z(" << l_shift_comp[2] 
	 << "," << k_shift_comp[2] << "," << j_shift_comp[2] << "," 
	 << i_shift_comp[2] << ") = " << diff_shift_comp[2] << endl ;
      cout << "average difference for shift_comp_z : " << moy_shift_comp[2]
	  << endl << endl ;
 }

  // shift

  if (diff_shift[0] == -1) {
    cout << "shift_x = 0" << endl ; 
  }
  else {
    cout << "max relative difference for shift_x(" << l_shift[0] 
	 << "," << k_shift[0] << "," << j_shift[0] << "," 
	 << i_shift[0] << ") = " << diff_shift[0] << endl ;
      cout << "average difference for shift_x : " << moy_shift[0]
	  << endl ;
 }
  
  if (diff_shift[1] == -1) {
    cout << "shift_y = 0" << endl ; 
  }
  else {
    cout << "max relative difference for shift_y(" << l_shift[1] 
	 << "," << k_shift[1] << "," << j_shift[1] << "," 
	 << i_shift[1] << ") = " << diff_shift[1] << endl ;
     cout << "average difference for shift_y : " << moy_shift[1]
	  << endl ; 
 }
  
  if (diff_shift[2] == -1) {
    cout << "shift_z = 0" << endl << endl ; 
  }
  else {
    cout << "max relative difference for shift_z(" << l_shift[2] 
	 << "," << k_shift[2] << "," << j_shift[2] << "," 
	 << i_shift[2] << ") = " << diff_shift[2] << endl ;
     cout << "average difference for shift_z : " << moy_shift[2]
	  << endl << endl ;
  }
  
 
  // a_car
  //--------

  if (diff_a_car == -1) {
    cout << "a_car = 0" << endl << endl ; 
  }
  else {
    cout << "max relative difference for a_car(" << l_a_car << "," 
	 << k_a_car << "," << j_a_car << "," << i_a_car << ") = "
	 << diff_a_car << endl ;
     cout << "average difference for a_car : " << moy_a_car << endl << endl ;
 }
  
  // beta
  //------

  if (diff_beta == -1) {
    cout << "beta = 0" << endl << endl ; 
  }
  else {
    cout << "max relative difference for beta(" << l_beta << "," 
	 << k_beta << "," << j_beta << "," << i_beta << ") = "
	 << diff_beta << endl  ;
    cout << "average difference for beta : " << moy_beta << endl << endl ;
  }

  // Hydro quantities
  //-----------------

  // ener

  if (diff_ener == -1) {
    cout << "ener = 0" << endl ; 
  }
  else {
    cout << "max relative difference for ener(" << l_ener << "," 
	 << k_ener << "," << j_ener << "," << i_ener << ") = "
	 << diff_ener << endl ;
    cout << "average difference for ener : " << moy_ener << endl ;
  }
 
  // press

  if (diff_press == -1) {
    cout << " press= 0" << endl ; 
  }
  else {
    cout << "max relative difference for press(" << l_press << "," 
	 << k_press << "," << j_press << "," << i_press << ") = "
	 << diff_press << endl ;
     cout << "average difference for press : " << moy_press << endl ;
 }
 
  // ener_euler

  if (diff_ener_euler == -1) {
    cout << "ener_euler = 0" << endl ; 
  }
  else {
    cout << "max relative difference for ener_euler(" << l_ener_euler<< "," 
	 << k_ener_euler << "," << j_ener_euler << "," << i_ener_euler << ") = "
	 << diff_ener_euler << endl ;
     cout << "average difference for ener_euler : " << moy_ener_euler << endl ;
 }
 
  // s_euler

  if (diff_s_euler == -1) {
    cout << "s_euler = 0" << endl ; 
  }
  else {
    cout << "max relative difference for s_euler(" << l_s_euler << "," 
	 << k_s_euler << "," << j_s_euler << "," << i_s_euler << ") = "
	 << diff_s_euler << endl ;
     cout << "average difference for s_euler : " << moy_s_euler << endl ;
 }
 
  // gam_euler

  if (diff_gam_euler == -1) {
    cout << "gam_euler = 0" << endl ; 
  }
  else {
    cout << "max relative difference for gam_euler(" << l_gam_euler << "," 
	 << k_gam_euler << "," << j_gam_euler << "," << i_gam_euler << ") = "
	 << diff_gam_euler << endl ;
     cout << "average difference for gam_euler : " << moy_gam_euler << endl ;
 }
 
  // u_euler

  if (diff_u_euler[0] == -1) {
    cout << "u_euler_x = 0" << endl ; 
  }
  else {
    cout << "max relative difference for u_euler_x(" << l_u_euler[0] 
	 << "," << k_u_euler[0] << "," << j_u_euler[0] << "," 
	 << i_u_euler[0] << ") = " << diff_u_euler[0] << endl ;
    cout << "average difference for u_euler_x : " << moy_u_euler[0] << endl ;
  }
 
  if (diff_u_euler[1] == -1) {
    cout << "u_euler_y = 0" << endl ; 
  }
  else {
    cout << "max relative difference for u_euler_y(" << l_u_euler[1] 
	 << "," << k_u_euler[1] << "," << j_u_euler[1] << "," 
	 << i_u_euler[1] << ") = " << diff_u_euler[1] << endl ;
     cout << "average difference for u_euler_y : " << moy_u_euler[1] << endl ;
 }

  if (diff_u_euler[2] == -1) {
    cout << "u_euler_z = 0" << endl << endl ; 
  }
  else {
    cout << "max relative difference for u_euler_z(" << l_u_euler[2] 
	 << "," << k_u_euler[2] << "," << j_u_euler[2] << "," 
	 << i_u_euler[2] << ") = " << diff_u_euler[2] << endl ;
     cout << "average difference for u_euler_z : " << moy_u_euler[2] << endl 
	  << endl ;
 }


  // Miscallaneous
  //---------------

  // tkij_auto

  if (diff_tkij_auto[0] == -1) {
    cout << "tkij_auto(0,0) = 0" << endl ; 
  }
  else {
    cout << "max relative difference for tkij_auto(0,0)(" << l_tkij_auto[0] 
	 << "," << k_tkij_auto[0] << "," << j_tkij_auto[0] << "," 
	 << i_tkij_auto[0] << ") = " << diff_tkij_auto[0] << endl ;
    cout << "average difference for tkij_auto(0,0) : " << moy_tkij_auto[0] 
	 << endl ;
  }
  
  if (diff_tkij_auto[1] == -1) {
    cout << "tkij_auto(0,1) = 0" << endl ; 
  }
  else {
    cout << "max relative difference for tkij_auto(0,1)(" << l_tkij_auto[1] 
	 << "," << k_tkij_auto[1] << "," << j_tkij_auto[1] << "," 
	 << i_tkij_auto[1] << ") = " << diff_tkij_auto[1] << endl ;
    cout << "average difference for tkij_auto(0,1) : " << moy_tkij_auto[1] 
	 << endl ;
 }
  
  if (diff_tkij_auto[2] == -1) {
    cout << "tkij_auto(0,2) = 0" << endl ; 
  }
  else {
    cout << "max relative difference for tkij_auto(0,2)(" << l_tkij_auto[2] 
	 << "," << k_tkij_auto[2] << "," << j_tkij_auto[2] << "," 
	 << i_tkij_auto[2] << ") = " << diff_tkij_auto[2] << endl ;
    cout << "average difference for tkij_auto(0,2) : " << moy_tkij_auto[2] 
	 << endl ;
  }

  if (diff_tkij_auto[3] == -1) {
    cout << "tkij_auto(1,0) = 0" << endl ; 
  }
  else {
    cout << "max relative difference for tkij_auto(1,0)(" << l_tkij_auto[3] 
	 << "," << k_tkij_auto[3] << "," << j_tkij_auto[3] << "," 
	 << i_tkij_auto[3] << ") = " << diff_tkij_auto[3] << endl ;
    cout << "average difference for tkij_auto(1,0) : " << moy_tkij_auto[3] 
	 << endl ;
  }
  
  if (diff_tkij_auto[4] == -1) {
    cout << "tkij_auto(1,1) = 0" << endl ; 
  }
  else {
    cout << "max relative difference for tkij_auto(1,1)(" << l_tkij_auto[4] 
	 << "," << k_tkij_auto[4] << "," << j_tkij_auto[4] << "," 
	 << i_tkij_auto[4] << ") = " << diff_tkij_auto[4] << endl ;
    cout << "average difference for tkij_auto(1,1) : " << moy_tkij_auto[4] 
	 << endl ;
  }
  
  if (diff_tkij_auto[5] == -1) {
    cout << "tkij_auto(1,2) = 0" << endl ; 
  }
  else {
    cout << "max relative difference for tkij_auto(1,2)(" << l_tkij_auto[5] 
	 << "," << k_tkij_auto[5] << "," << j_tkij_auto[5] << "," 
	 << i_tkij_auto[5] << ") = " << diff_tkij_auto[5] << endl ;
     cout << "average difference for tkij_auto(1,2) : " << moy_tkij_auto[5] 
	 << endl ;
 }

  if (diff_tkij_auto[6] == -1) {
    cout << "tkij_auto(2,0) = 0" << endl ; 
  }
  else {
    cout << "max relative difference for tkij_auto(2,0)(" << l_tkij_auto[6] 
	 << "," << k_tkij_auto[6] << "," << j_tkij_auto[6] << "," 
	 << i_tkij_auto[6] << ") = " << diff_tkij_auto[6] << endl ;
    cout << "average difference for tkij_auto(2,0) : " << moy_tkij_auto[6] 
	 << endl ;
  }
  
  if (diff_tkij_auto[7] == -1) {
    cout << "tkij_auto(2,1) = 0" << endl ; 
  }
  else {
    cout << "max relative difference for tkij_auto(2,1)(" << l_tkij_auto[7] 
	 << "," << k_tkij_auto[7] << "," << j_tkij_auto[7] << "," 
	 << i_tkij_auto[7] << ") = " << diff_tkij_auto[7] << endl ;
    cout << "average difference for tkij_auto(2,1) : " << moy_tkij_auto[7] 
	 << endl ;
  }

  if (diff_tkij_auto[8] == -1) {
    cout << "tkij_auto(2,2) = 0" << endl << endl ; 
  }
  else {
    cout << "max relative difference for tkij_auto(2,2)(" << l_tkij_auto[8] 
	 << "," << k_tkij_auto[8] << "," << j_tkij_auto[8] << "," 
	 << i_tkij_auto[8] << ") = " << diff_tkij_auto[8] << endl ;
    cout << "average difference for tkij_auto(2,2) : " << moy_tkij_auto[8] 
	 << endl << endl ;
  }


  // tkij_comp

  if (diff_tkij_comp[0] == -1) {
    cout << "tkij_comp(0,0) = 0" << endl ; 
  }
  else {
    cout << "max relative difference for tkij_comp(0,0)(" << l_tkij_comp[0] 
	 << "," << k_tkij_comp[0] << "," << j_tkij_comp[0] << "," 
	 << i_tkij_comp[0] << ") = " << diff_tkij_comp[0] << endl ;
   cout << "average difference for tkij_comp(0,0) : " << moy_tkij_comp[0] 
	 << endl ;
  }

  if (diff_tkij_comp[1] == -1) {
    cout << "tkij_comp(0,1) = 0" << endl ; 
  }
  else {
    cout << "max relative difference for tkij_comp(0,1)(" << l_tkij_comp[1] 
	 << "," << k_tkij_comp[1] << "," << j_tkij_comp[1] << "," 
	 << i_tkij_comp[1] << ") = " << diff_tkij_comp[1] << endl ;
   cout << "average difference for tkij_comp(0,1) : " << moy_tkij_comp[1] 
	 << endl ;
  }

  if (diff_tkij_comp[2] == -1) {
    cout << "tkij_comp(0,2) = 0" << endl ; 
  }
  else {
    cout << "max relative difference for tkij_comp(0,2)(" << l_tkij_comp[2] 
	 << "," << k_tkij_comp[2] << "," << j_tkij_comp[2] << "," 
	 << i_tkij_comp[2] << ") = " << diff_tkij_comp[2] << endl ;
   cout << "average difference for tkij_comp(0,2) : " << moy_tkij_comp[2] 
	 << endl ;
  }

  if (diff_tkij_comp[3] == -1) {
    cout << "tkij_comp(1,0) = 0" << endl ; 
  }
  else {
    cout << "max relative difference for tkij_comp(1,0)(" << l_tkij_comp[3] 
	 << "," << k_tkij_comp[3] << "," << j_tkij_comp[3] << "," 
	 << i_tkij_comp[3] << ") = " << diff_tkij_comp[3] << endl ;
    cout << "average difference for tkij_comp(1,0) : " << moy_tkij_comp[3] 
	 << endl ;
 }

  if (diff_tkij_comp[4] == -1) {
    cout << "tkij_comp(1,1) = 0" << endl ; 
  }
  else {
    cout << "max relative difference for tkij_comp(1,1)(" << l_tkij_comp[4] 
	 << "," << k_tkij_comp[4] << "," << j_tkij_comp[4] << "," 
	 << i_tkij_comp[4] << ") = " << diff_tkij_comp[4] << endl ;
   cout << "average difference for tkij_comp(1,1) : " << moy_tkij_comp[4] 
	 << endl ;
  }


  if (diff_tkij_comp[5] == -1) {
    cout << "tkij_comp(1,2) = 0" << endl ; 
  }
  else {
    cout << "max relative difference for tkij_comp(1,2)(" << l_tkij_comp[5] 
	 << "," << k_tkij_comp[5] << "," << j_tkij_comp[5] << "," 
	 << i_tkij_comp[5] << ") = " << diff_tkij_comp[5] << endl ;
    cout << "average difference for tkij_comp(1,2) : " << moy_tkij_comp[5] 
	 << endl ;
 }


  if (diff_tkij_comp[6] == -1) {
    cout << "tkij_comp(2,0) = 0" << endl ; 
  }
  else {
    cout << "max relative difference for tkij_comp(2,0)(" << l_tkij_comp[6] 
	 << "," << k_tkij_comp[6] << "," << j_tkij_comp[6] << "," 
	 << i_tkij_comp[6] << ") = " << diff_tkij_comp[6] << endl ;
    cout << "average difference for tkij_comp(2,0) : " << moy_tkij_comp[6] 
	 << endl ;
 }

  if (diff_tkij_comp[7] == -1) {
    cout << "tkij_comp(2,1) = 0" << endl ; 
  }
  else {
    cout << "max relative difference for tkij_comp(2,1)(" << l_tkij_comp[7] 
	 << "," << k_tkij_comp[7] << "," << j_tkij_comp[7] << "," 
	 << i_tkij_comp[7] << ") = " << diff_tkij_comp[7] << endl ;
   cout << "average difference for tkij_comp(2,1) : " << moy_tkij_comp[7] 
	 << endl ;
  }

  if (diff_tkij_comp[8] == -1) {
    cout << "tkij_comp(2,2) = 0" << endl << endl ; 
  }
  else {
    cout << "max relative difference for tkij_comp(2,2)(" << l_tkij_comp[8] 
	 << "," << k_tkij_comp[8] << "," << j_tkij_comp[8] << "," 
	 << i_tkij_comp[8] << ") = " << diff_tkij_comp[8] << endl ;
    cout << "average difference for tkij_comp(2,2) : " << moy_tkij_comp[8] 
	 << endl << endl ;
 }

  // kcar_auto

  if (diff_kcar_auto == -1) {
    cout << "kcar_auto = 0" << endl << endl ; 
  }
  else {
    cout << "max relative difference for kcar_auto(" << l_kcar_auto << "," 
	 << k_kcar_auto << "," << j_kcar_auto << "," << i_kcar_auto << ") = "
	 << diff_kcar_auto << endl ;
     cout << "average difference for kcar_auto : " << moy_kcar_auto 
	  << endl << endl ;
  }
 
  // kcar_comp

  if (diff_kcar_comp == -1) {
    cout << "kcar_comp = 0" << endl << endl ; 
  }
  else {
    cout << "max relative difference for kcar_comp(" << l_kcar_comp << "," 
	 << k_kcar_comp << "," << j_kcar_comp << "," << i_kcar_comp << ") = "
	 << diff_kcar_comp << endl ;
     cout << "average difference for kcar_comp : " << moy_kcar_comp 
	  << endl << endl ;
  }
  
  cout << "relative difference for mass_b : " << diff_mass_b << endl ;
  cout << "relative difference for mass_g : " << diff_mass_g << endl ;
  cout << "relative difference for ray_eq : " << diff_ray_eq << endl ;
  cout << "relative difference for ray_pole : " << diff_ray_pole << endl ;
  cout << "relative difference for omega : " << diff_omega << endl ;
  if (diff_x_axe == -1) {
    cout << "x_axe = 0" << endl ;
  }
  else {
  cout << "relative difference for x_axe : " << diff_x_axe << endl ;
  }
  cout << "---------------------------------------------------" << endl ;
  
  }

}

void Bin_ns_ncp::compare(FILE* fich) {

  int mer_ini ; 
  fread(&mer_ini, sizeof(int), 1, fich) ;	

  Mg3d mg1(fich) ;
  Map_et mp1(mg1, fich) ;
  Eos* peos1 = Eos::eos_from_file(fich) ; 

  Mg3d mg2(fich) ;
  Map_et mp2(mg2, fich) ;
  Eos* peos2 = Eos::eos_from_file(fich) ; 

  //Binaire star(mp1, *peos1, mp2, *peos2, fich) ; 
  Binaire star(et[0]->set_mp(), *peos1, et[1]->set_mp(), *peos2, fich) ; 

  for (int i=1; i<=2; i++) {
    (star.set(i)).update_metric( (star(3-i)) ) ;
  }
  for (int i=1; i<=2; i++) {
    (star.set(i)).update_metric_der_comp( (star(3-i)) ) ;
  }
  for (int i=1; i<=2; i++) {
  (star.set(i)).equation_of_state() ;
  (star.set(i)).kinematics(star.get_omega(), star.get_x_axe()) ; 
  (star.set(i)).fait_d_psi() ; 
  (star.set(i)).hydro_euler() ; 
  }

  compare( star ) ;

}
