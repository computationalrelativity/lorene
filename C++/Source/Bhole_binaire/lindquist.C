/*
 *   Copyright (c) 2000-2001 Philippe Grandclement
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


char lindquist_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/10/16 14:36:33  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  2000/12/13  15:42:57  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


//standard
#include <stdlib.h>
#include <math.h>

// LORENE
#include "type_parite.h"
#include "nbr_spx.h"
#include "proto.h"
#include "coord.h"
#include "tenseur.h"

double serie_lindquist_plus (double rayon, double distance, double xa, double ya, 
	double za, double precision, double itemax) {
	    
	    
    double result = 0.5 ;
    double c_n = rayon ;
    double d_n = distance/2. ;
    
    int indic = 1 ;
    int conte=0 ;
    while ((indic == 1) && (conte <= itemax)) {
	
	double norme_plus = sqrt ((xa+d_n)*(xa+d_n)+ya*ya+za*za) ;
		
	double terme = c_n * 1./norme_plus ;
	if (fabs(terme/result) < precision)
	    indic = -1 ;
    
        result = result + terme ;
	
	c_n *= rayon/(distance/2. + d_n) ;
	d_n = distance/2. - rayon*rayon/(distance/2.+d_n) ;
	conte ++ ;
    }

    if (conte > itemax)
	result = 1 ;
return result ;
}

double serie_lindquist_moins (double rayon, double distance, double xa, double ya, 
	double za, double precision, double itemax) {
	    
	    
    double result = 0.5 ;
    double c_n = rayon ;
    double d_n = distance/2. ;
    
    int indic = 1 ;
    int conte=0 ;
    while ((indic == 1) && (conte <= itemax)) {
	
	double norme_plus = sqrt ((xa-d_n)*(xa-d_n)+ya*ya+za*za) ;
		
	double terme = c_n * 1./norme_plus ;
	if (fabs(terme/result) < precision)
	    indic = -1 ;
    
        result = result + terme ;
	
	c_n *= rayon/(distance/2. + d_n) ;
	d_n = distance/2. - rayon*rayon/(distance/2.+d_n) ;
	conte ++ ;
    }

    if (conte > itemax)
	result = 1 ;
return result ;
}

double adm_serie (double rayon, double distance, double precision) {
    
    double result = 0 ;
    double c_n = rayon ;
    double d_n = distance/2. ;
    
    int indic =1 ;
    while (indic == 1) {
	
	result += 4*c_n ;
	if (fabs(c_n/result) < precision)
	    indic = -1 ;
	
	c_n *= rayon/(distance/2. + d_n) ;
	d_n = distance/2. - rayon*rayon/(distance/2.+d_n) ;
    }
return result ;
}

double bare_serie (double rayon, double distance, double precision) {
    
    double result = 0 ;
    double c_n = rayon ;
    double d_n = distance/2. ;
    
    int indic =1 ;
    int n = 1 ;
    while (indic == 1) {
	
	result += 4*c_n*n ;
	if (fabs(c_n/result) < precision)
	    indic = -1 ;
	
	c_n *= rayon/(distance/2. + d_n) ;
	d_n = distance/2. - rayon*rayon/(distance/2.+d_n) ;
	n++ ;
    }
return result ;
}

void set_lindquist (Cmp& psi_un, Cmp& psi_deux, double rayon, double precision) {
    
    // Pour les allocations !
    psi_un = 1 ;
    psi_deux = 1 ;
    
    double distance = psi_un.get_mp()->get_ori_x()-psi_deux.get_mp()->get_ori_x() ;
    
    // On regarde psi 1 :
    Mtbl xa_mtbl_un (psi_un.get_mp()->xa) ;
    Mtbl ya_mtbl_un (psi_un.get_mp()->ya) ;
    Mtbl za_mtbl_un (psi_un.get_mp()->za) ;
    
    int nz_un = psi_un.get_mp()->get_mg()->get_nzone() ;
    for (int l=1 ; l<nz_un ; l++) {
	int np = psi_un.get_mp()->get_mg()->get_np (l) ;
	int nt = psi_un.get_mp()->get_mg()->get_nt (l) ;
	int nr = psi_un.get_mp()->get_mg()->get_nr (l) ;
	double xa, ya, za ;
	for (int k=0 ; k<np ; k++) 
	    for (int j=0 ; j<nt ; j++)
		for (int i=0 ; i<nr ; i++) {
		    xa = xa_mtbl_un (l, k, j, i) ;
		    ya = ya_mtbl_un (l, k, j, i) ;
		    za = za_mtbl_un (l, k, j, i) ;
		    
		    psi_un.set(l, k, j, i) = 
serie_lindquist_moins (rayon, distance, xa, ya, za, precision, 30) ;
		}
    }
    
    psi_un.set_val_inf (0.5) ;
    psi_un.std_base_scal() ;
    
      // On regarde psi 2 :
    Mtbl xa_mtbl_deux (psi_deux.get_mp()->xa) ;
    Mtbl ya_mtbl_deux (psi_deux.get_mp()->ya) ;
    Mtbl za_mtbl_deux (psi_deux.get_mp()->za) ;
    
    int nz_deux = psi_deux.get_mp()->get_mg()->get_nzone() ;
    for (int l=1 ; l<nz_deux ; l++) {
	int np = psi_deux.get_mp()->get_mg()->get_np (l) ;
	int nt = psi_deux.get_mp()->get_mg()->get_nt (l) ;
	int nr = psi_deux.get_mp()->get_mg()->get_nr (l) ;
	double xa, ya, za ;
	for (int k=0 ; k<np ; k++) 
	    for (int j=0 ; j<nt ; j++)
		for (int i=0 ; i<nr ; i++) {
		    xa = xa_mtbl_deux (l, k, j, i) ;
		    ya = ya_mtbl_deux (l, k, j, i) ;
		    za = za_mtbl_deux (l, k, j, i) ;
		    
		    psi_deux.set(l, k, j, i) = 
serie_lindquist_plus (rayon, distance, xa, ya, za, precision, 30) ;
		}
    }
    psi_deux.set_val_inf (0.5) ;
    psi_deux.std_base_scal() ;
}
