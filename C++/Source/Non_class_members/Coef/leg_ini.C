/*
 *   Copyright (c) 2013 Jerome Novak
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


char leg_ini_C[] = "$Header$" ;


/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:53:13  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2013/06/06 15:31:33  j_novak
 * Functions to compute Legendre coefficients (not fully tested yet).
 *
 * Revision 1.1  2013/06/05 15:08:13  j_novak
 * Initial revision. Not ready yet...
 *
 *
 *
 * $Header$
 *
 */


// headers du C
#include <cstdlib>
#include <cassert>
#include <cmath>

//Lorene prototypes
#include "tbl.h"
#include "utilitaires.h"

namespace Lorene {

namespace {
  const int nmax = 50 ; //Maximal number of Legendre transforms sizes 
  int nwork_colloc = 0 ;
  int nwork_leg = 0 ;
  double* tab_colloc[nmax] ;
  int nb_colloc[nmax] ;
  Tbl* tab_pni[nmax] ;
  Tbl* tab_wn[nmax] ;
  int nb_leg[nmax] ;
}

void poly_leg (int n, double& poly, double& pder, double& polym1, double& pderm1, 
		double& polym2, double& pderm2, double x) {
		
	
	if (n==0) {
	     poly = 1 ;
	     pder = 0 ;
	     }
	else 
	    if (n==1) {
	         polym1 = 1 ;
		 pderm1 = 0 ;
		 poly = x ;
		 pder = 1 ;
		 }
	else {
	     polym1 = 1 ;
	     pderm1 = 0 ;
	     poly = x ;
	     pder = 1 ;
	     for (int i=1 ; i<n ; i++) {
	         polym2 = polym1 ;
		 pderm2 = pderm1 ;
		 polym1 = poly ;
		 pderm1 = pder ;
		 poly = ((2*i+1)*x*polym1 - i*polym2)/(i+1) ;
		 pder = ((2*i+1)*polym1+(2*i+1)*x*pderm1-i*pderm2)/(i+1) ;
		}
	}
}

/************************************************************************/
void legendre_collocation_points(int nr, double* colloc) {

  int index = -1 ;
  for (int i=0; ((i<nwork_colloc) && (index<0)); i++) 
    if (nb_colloc[i] == nr) index = i ; //Have the collocation points already been 
                                        // computed?

  if (index <0) { //New array needed
    index = nwork_colloc ;
    if (index >= nmax) {
      cout << "legendre_collocation_points: " << endl ;
      cout << "too many arrays!" << endl ;
      abort() ;
    }
    double*& t_colloc = tab_colloc[index] ;
    t_colloc = new double[nr] ;
    int nr0 = nr - 1 ;
	
    double x_plus = 1 ;
    double x_moins = -1 ;
    double p_plus, dp_plus, p_plus_m1, dp_plus_m1, p_plus_m2, dp_plus_m2 ;
    double p_moins, dp_moins, p_moins_m1, dp_moins_m1, p_moins_m2, dp_moins_m2 ;
    double p, dp, p_m1, dp_m1, p_m2, dp_m2 ;
    
    poly_leg (nr, p_plus, dp_plus, p_plus_m1, dp_plus_m1, p_plus_m2, 
	      dp_plus_m2, x_plus) ;
    poly_leg (nr, p_moins, dp_moins, p_moins_m1, dp_moins_m1, p_moins_m2, 
	      dp_moins_m2, x_moins) ;
    
    double det = p_plus_m1*p_moins_m2 - p_moins_m1*p_plus_m2 ;
    double r_plus = -p_plus ;
    double r_moins = -p_moins ;
    double a = (r_plus*p_moins_m2 - r_moins*p_plus_m2)/det ;
    double b = (r_moins*p_plus_m1 - r_plus*p_moins_m1)/det ;
    
    t_colloc[nr0] = 1 ;
    double dth = M_PI/(2*nr+1) ;
    double cd = cos (2*dth) ;
    double sd = sin (2*dth) ;
    double cs = cos(dth) ;
    double ss = sin(dth) ;
     
    int borne_sup = (nr%2==0) ? nr/2 : (nr+1)/2  ;
   
    for (int j=1 ; j<borne_sup ; j++) {
      double x = cs ;
      bool loop = true ;
      int ite = 0 ;
      while (loop) {
	poly_leg (nr, p, dp, p_m1, dp_m1, p_m2, dp_m2, x) ;
	double poly = p + a*p_m1 + b*p_m2 ;
	double pder = dp + a * dp_m1 + b*dp_m2 ;
	double sum = 0 ;
	for (int i=0 ; i<j ; i++)
	  sum += 1./(x-t_colloc[nr-i-1]) ;
		   
	double increm = -poly/(pder-sum*poly) ;
	       
	x += increm ;
	ite ++ ;
	if ((fabs(increm) < 1.e-14) || (ite >500))
	  loop = false ;
      }
      if (ite > 500) {
	cout << "leg_ini: too many iterations..." << endl ;
	abort() ;
      }
      t_colloc[nr-j-1] = x ;
      double auxi = cs*cd-ss*sd ;
      ss = cs*sd+ss*cd ;
      cs = auxi ;
    }
    if  (nr%2==1)
      t_colloc[(nr-1)/2] = 0 ; 
    // Copy of the symetric ones :
    for (int i=0 ; i<borne_sup ; i++)
      t_colloc[i] = - t_colloc[nr-i-1] ;
    nb_colloc[index] = nr ;
    nwork_colloc++ ;
  }
  assert((index>=0)&&(index<nmax)) ;
  for (int i=0; i<nr; i++)
    colloc[i] = (tab_colloc[index])[i] ;

  return ;

}



/************************************************************************/

void get_legendre_data(int np, Tbl*& p_Pni, Tbl*& p_wn) {

  int index = -1 ;
  for (int i=0; ((i<nwork_leg) && (index<0)); i++) 
    if (nb_leg[i] == np) index = i ; //Has the plan already been estimated?

  if (index <0) { //New plan needed
    index = nwork_leg ;
    if (index >= nmax) {
      cout << "get_legendre_data: " << endl ;
      cout << "too many transformation matrices!" << endl ;
      abort() ;
    }
    int np0 = np - 1 ;
    tab_pni[index] = new Tbl(np, np) ;
    Tbl& Pni = (*tab_pni[index]) ;
    Pni.set_etat_qcq() ;
    tab_wn[index] = new Tbl(np) ;
    Tbl& wn = (*tab_wn[index]) ;
    wn.set_etat_qcq() ;

    Tbl coloc(np) ;
    coloc.set_etat_qcq() ;
    legendre_collocation_points(np, coloc.t) ;
      
    for (int i=0; i<np; i++)  {
      Pni.set(0, i) = 1 ;
      Pni.set(1, i) = coloc(i) ;
    }
    for (int n=2; n<np; n++) {
      for (int i=0; i<np; i++) {
	Pni.set(n,i) = (2. - 1./double(n))*coloc(i)*Pni(n-1, i) 
	  - (1. - 1./double(n))*Pni(n-2, i) ;
      }
    }

    for (int j=0; j<np; j++)
      wn.set(j) = 2./( double(np0)*double(np) * Pni(np0,j) * Pni(np0, j) ) ;

    nb_leg[index] = np ;
    nwork_leg++ ;
  }
  assert((index>=0)&&(index<nmax)) ;
  p_Pni = tab_pni[index] ;
  p_wn = tab_wn[index] ;

  return ;
}

}
