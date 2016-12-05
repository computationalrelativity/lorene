/*
 *  Builds the operator for Ope_pois_tens_rr
 *
 *    (see file ope_elementary.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Jerome Novak
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.4  2016/12/05 16:18:12  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:53:34  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:16:13  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2004/12/23 16:30:15  j_novak
 * New files and class for the solution of the rr component of the tensor Poisson
 * equation.
 *
 *
 *
 * $Header$
 *
 */

//fichiers includes
#include <cstdlib>

#include "matrice.h"
#include "type_parite.h"
#include "proto.h"

/*
 * Routine caluclant l'operateur suivant :
 * 
 * -R_CHEB : r^2d2sdx2+6*rdsdx+6-l*(l+1)Id
 * 
 * -R_CHEBP et R_CHEBI : d2sdx2+6/r dsdx +(6-l(l+1))/r^2
 * 
 * -R_CHEBU : d2sdx2+(6-l(l+1))/x^2
 * 
 * Entree :
 *	-n nbre de points en r
 *      -l voire operateur. (doit etre >=2)
 *      -echelle utile uniquement pour R_CHEB : represente beta/alpha 
 *        (cf mapping)	
 *      - puis : exposant de multiplication dans la ZEC
 *      - base_r : base de developpement
 *  Sortie :
 *  La fonction renvoie la matrice.	
 */
 
namespace Lorene {
Matrice _ope_ptens_rr_mat_pas_prevu(int, int, double, int) ;
Matrice _ope_ptens_rr_mat_r_chebp(int, int, double, int) ;
Matrice _ope_ptens_rr_mat_r_chebi(int, int, double, int) ;
Matrice _ope_ptens_rr_mat_r_chebu(int, int, double, int) ;
Matrice _ope_ptens_rr_mat_r_cheb(int, int, double, int) ;

                       //-----------------------------------
                       // Routine pour les cas non prevus --
                       //-----------------------------------

Matrice _ope_ptens_rr_mat_pas_prevu(int n, int l, double echelle, int puis) {
  cout << "laplacien tens_rr pas prevu..." << endl ;
  cout << "n : " << n << endl ;
  cout << "l : " << l << endl ;
  cout << "puissance : " << puis << endl ;
  cout << "echelle : " << echelle << endl ;
  abort() ;
  exit(-1) ;
  Matrice res(1, 1) ;
  return res;
}


                             //-------------------------
                             //----  CAS R_CHEBP    ----
                             //-------------------------
		    

Matrice _ope_ptens_rr_mat_r_chebp (int n, int l, double, int) {
   
  const int nmax = 100 ;// Nombre de Matrices stockees
  static Matrice* tab[nmax] ;  // les matrices calculees
  static int nb_dejafait = 0 ; // nbre de matrices calculees
  static int l_dejafait[nmax] ;
  static int nr_dejafait[nmax] ;
    
  int indice = -1 ;
   
  // On determine si la matrice a deja ete calculee :
  for (int conte=0 ; conte<nb_dejafait ; conte ++)
    if ((l_dejafait[conte] == l) && (nr_dejafait[conte] == n))
      indice = conte ;
    
  // Calcul a faire : 
  if (indice  == -1) {
    if (nb_dejafait >= nmax) {
      cout << "_ope_ptens_rr_mat_r_chebp : trop de matrices" << endl ;
      abort() ;
      exit (-1) ;
    }
       

    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice dd(n, n) ;
    dd.set_etat_qcq() ;
    Matrice xd(n, n) ;
    xd.set_etat_qcq() ;
    Matrice xx(n, n) ;
    xx.set_etat_qcq() ;

    double* vect  = new double[n] ;
    
    for (int i=0 ; i<n ; i++) {
      for (int j=0 ; j<n ; j++)
	vect[j] = 0 ;
      vect[i] = 1 ;
      d2sdx2_1d (n, &vect, R_CHEBP) ;
	
      for (int j=0 ; j<n ; j++)
	dd.set(j, i) = vect[j] ; 
    }
    
    for (int i=0 ; i<n ; i++) {
      for (int j=0 ; j<n ; j++)
	vect[j] = 0 ;
      vect[i] = 1 ;
      sxdsdx_1d (n, &vect, R_CHEBP) ;
      for (int j=0 ; j<n ; j++)
	xd.set(j, i) = vect[j] ;
	
    }
    
    for (int i=0 ; i<n ; i++) {
      for (int j=0 ; j<n ; j++)
	vect[j] = 0 ;
      vect[i] = 1 ;
      sx2_1d (n, &vect, R_CHEBP) ;
      for (int j=0 ; j<n ; j++)
	xx.set(j, i) = vect[j] ;
    }
   
    delete [] vect ;
    
    Matrice res(n, n) ;
    res = dd+6*xd+(6-l*(l+1))*xx ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;	
    return res ;
  }
    
  // Cas ou le calcul a deja ete effectue :
  else
    return *tab[indice] ;
}



                             //------------------------
                             //--   CAS R_CHEBI    ----
                             //------------------------
		    

Matrice _ope_ptens_rr_mat_r_chebi (int n, int l, double, int) {
   
  const int nmax = 100 ;// Nombre de Matrices stockees
  static Matrice* tab[nmax] ;  // les matrices calculees
  static int nb_dejafait = 0 ; // nbre de matrices calculees
  static int l_dejafait[nmax] ;
  static int nr_dejafait[nmax] ;
    
  int indice = -1 ;
   
  // On determine si la matrice a deja ete calculee :
  for (int conte=0 ; conte<nb_dejafait ; conte ++)
    if ((l_dejafait[conte] == l) && (nr_dejafait[conte] == n))
      indice = conte ;
    
  // Calcul a faire : 
  if (indice  == -1) {
    if (nb_dejafait >= nmax) {
      cout << "_ope_ptens_rr_mat_r_chebi : trop de matrices" << endl ;
      abort() ;
      exit (-1) ;
    }
       
	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice dd(n, n) ;
    dd.set_etat_qcq() ;
    Matrice xd(n, n) ;
    xd.set_etat_qcq() ;
    Matrice xx(n, n) ;
    xx.set_etat_qcq() ;

    double* vect = new double[n] ;
    
    for (int i=0 ; i<n ; i++) {
      for (int j=0 ; j<n ; j++)
	vect[j] = 0 ;
      vect[i] = 1 ;
      d2sdx2_1d (n, &vect, R_CHEBI) ;  // appel dans le cas impair
      for (int j=0 ; j<n ; j++)
	dd.set(j, i) = vect[j] ;
    }
    
    for (int i=0 ; i<n ; i++) {
      for (int j=0 ; j<n ; j++)
	vect[j] = 0 ;
      vect[i] = 1 ;
      sxdsdx_1d (n, &vect, R_CHEBI) ;
      for (int j=0 ; j<n ; j++)
	xd.set(j, i) = vect[j] ;
    }
    
    for (int i=0 ; i<n ; i++) {
      for (int j=0 ; j<n ; j++)
	vect[j] = 0 ;
      vect[i] = 1 ;
      sx2_1d (n, &vect, R_CHEBI) ;
      for (int j=0 ; j<n ; j++)
	xx.set(j, i) = vect[j] ;
    }
    
    delete [] vect ;
    
    Matrice res(n, n) ;
    res = dd+6*xd+(6-l*(l+1))*xx ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;
    return res ;
  } 
    
  // Cas ou le calcul a deja ete effectue :
  else
    return *tab[indice] ;
}




                            //------------------------
                            //----  CAS R_CHEBU   ----
                            //------------------------

Matrice _ope_ptens_rr_mat_r_chebu( int n, int l, double, int puis) {

  if (puis != 4) {
    cout << "_ope_ptens_rr_mat_r_chebu : only the case dzpuis = 4 "
	 << '\n' << "is implemented! \n"
	 << "dzpuis = " << puis << endl ;
    abort() ;
  }
        
  const int nmax = 200 ;// Nombre de Matrices stockees
  static Matrice* tab[nmax] ;  // les matrices calculees
  static int nb_dejafait = 0 ; // nbre de matrices calculees
  static int l_dejafait[nmax] ;
  static int nr_dejafait[nmax] ;
    
  int indice = -1 ;
   
  // On determine si la matrice a deja ete calculee :
  for (int conte=0 ; conte<nb_dejafait ; conte ++)
    if ((l_dejafait[conte] == l) && (nr_dejafait[conte] == n))
      indice = conte ;
    
  // Calcul a faire : 
  if (indice  == -1) {
    if (nb_dejafait >= nmax) {
      cout << "_ope_ptens_rr_mat_r_chebu : trop de matrices" << endl ;
      abort() ;
      exit (-1) ;
    }
       
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
      
    Matrice dd(n, n) ;
    dd.set_etat_qcq() ;
    Matrice xd(n, n) ;
    xd.set_etat_qcq() ;
    Matrice xx(n, n) ;
    xx.set_etat_qcq() ;
      
    double* vect = new double[n] ;
      
    for (int i=0 ; i<n ; i++) {
      for (int j=0 ; j<n ; j++)
	vect[j] = 0 ;
      vect[i] = 1 ;
      d2sdx2_1d (n, &vect, R_CHEBU) ;  // appel dans le cas unsurr
      for (int j=0 ; j<n ; j++)
	dd.set(j, i) = vect[j] ;
    }
      
    for (int i=0 ; i<n ; i++) {
      for (int j=0 ; j<n ; j++)
	vect[j] = 0 ;
      vect[i] = 1 ;
      dsdx_1d (n, &vect, R_CHEBU) ;  // appel dans le cas unsurr
      sxm1_1d_cheb (n, vect) ;
      for (int j=0 ; j<n ; j++)
	xd.set(j, i) = vect[j] ;
    }
      
    for (int i=0 ; i<n ; i++) {
      for (int j=0 ; j<n ; j++)
	vect[j] = 0 ;
      vect[i] = 1 ;
      sx2_1d (n, &vect, R_CHEBU) ;
      for (int j=0 ; j<n ; j++)
	xx.set(j, i) = vect[j] ;
    }
    
    delete [] vect ;
    
    Matrice res(n, n) ;
    res = dd - 4*xd + (6 -l*(l+1))*xx ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;
    return res ;
  } 
    
  // Cas ou le calcul a deja ete effectue :
  else
    return *tab[indice] ;
}


                   //-----------------------
                   //----  CAS R_CHEB   ----
                   //-----------------------
		    

Matrice _ope_ptens_rr_mat_r_cheb (int n, int l, double echelle, int) {
            
  const int nmax = 200 ;// Nombre de Matrices stockees
  static Matrice* tab[nmax] ;  // les matrices calculees
  static int nb_dejafait = 0 ; // nbre de matrices calculees
  static int l_dejafait[nmax] ;
  static int nr_dejafait[nmax] ;
  static double vieux_echelle = 0;
   
  // Si on a change l'echelle : on detruit tout :
  if (vieux_echelle != echelle) {
    for (int i=0 ; i<nb_dejafait ; i++) {
      l_dejafait[i] = -1 ;
      nr_dejafait[i] = -1 ;
      delete tab[i] ;
    }
       
    nb_dejafait = 0 ;
    vieux_echelle = echelle ;
  }
      
  int indice = -1 ;
   
  // On determine si la matrice a deja ete calculee :
  for (int conte=0 ; conte<nb_dejafait ; conte ++)
    if ((l_dejafait[conte] == l) && (nr_dejafait[conte] == n))
      indice = conte ;
    
  // Calcul a faire : 
  if (indice  == -1) {
    if (nb_dejafait >= nmax) {
      cout << "_ope_ptens_rr_mat_r_cheb : trop de matrices" << endl ;
      abort() ;
      exit (-1) ;
    }
       
    	
    l_dejafait[nb_dejafait] = l ;
    nr_dejafait[nb_dejafait] = n ;
    
    Matrice dd(n, n) ;
    dd.set_etat_qcq() ;
    Matrice xd(n, n) ;
    xd.set_etat_qcq() ;
    Matrice xx(n, n) ;
    xx.set_etat_qcq() ;

    double* vect = new double[n] ;
    
    for (int i=0 ; i<n ; i++) {
      for (int j=0 ; j<n ; j++)
	vect[j] = 0 ;
      vect[i] = 1 ;
      d2sdx2_1d (n, &vect, R_CHEB) ;  // appel dans le cas fin
      for (int j=0 ; j<n ; j++)
	dd.set(j, i) = vect[j]*echelle*echelle ;
    }
    
    for (int i=0 ; i<n ; i++) {
      for (int j=0 ; j<n ; j++)
	vect[j] = 0 ;
      vect[i] = 1 ;
      d2sdx2_1d (n, &vect, R_CHEB) ;  // appel dans le cas fin
      multx_1d (n, &vect, R_CHEB) ;
      for (int j=0 ; j<(n>i+1 ? i+1 : n) ; j++)
	dd.set(j, i) += 2*echelle*vect[j] ;
    }
    
    for (int i=0 ; i<n ; i++) {
      for (int j=0 ; j<n ; j++)
	vect[j] = 0 ;
      vect[i] = 1 ;
      d2sdx2_1d (n, &vect, R_CHEB) ;  // appel dans le cas fin
      multx_1d (n, &vect, R_CHEB) ;
      multx_1d (n, &vect, R_CHEB) ;
      for (int j=0 ; j<(n>i+1 ? i+1 : n) ; j++)
	dd.set(j, i) += vect[j] ;
    }
    
    for (int i=0 ; i<n ; i++) {
      for (int j=0 ; j<n ; j++)
	vect[j] = 0 ;
      vect[i] = 1 ;
      sxdsdx_1d (n, &vect, R_CHEB) ;
      for (int j=0 ; j<n ; j++)
	xd.set(j, i) = vect[j]*echelle ;
    }
    
    for (int i=0 ; i<n ; i++) {
      for (int j=0 ; j<n ; j++)
	vect[j] = 0 ;
      vect[i] = 1 ;
      sxdsdx_1d (n, &vect, R_CHEB) ;
      multx_1d (n, &vect, R_CHEB) ;
      for (int j=0 ; j<(n>i+1 ? i+1 : n) ; j++)
	xd.set(j, i) += vect[j] ;
    }
	   
    for (int i=0 ; i<n ; i++) {
      for (int j=0 ; j<n ; j++)
	vect[j] = 0 ;
      vect[i] = 1 ;
      sx2_1d (n, &vect, R_CHEB) ;
      for (int j=0 ; j<n ; j++)
	xx.set(j, i) = vect[j] ;
    }
    
    delete [] vect ;
    
    Matrice res(n, n) ;
    res = dd + 6*xd + (6 - l*(l+1))*xx ;
    tab[nb_dejafait] = new Matrice(res) ;
    nb_dejafait ++ ;
    return res ;
  } 
    
  // Cas ou le calcul a deja ete effectue :
  else
    return *tab[indice] ;  
}


                        //----------------------------
                        //--- La routine a appeler ---
                        //----------------------------

Matrice ope_ptens_rr_mat(int n, int l, double echelle, int puis, int base_r)
{

  // Routines de derivation
  static Matrice (*ope_ptens_rr_mat[MAX_BASE])(int, int, double, int) ;
  static int nap = 0 ;

  // Premier appel
  if (nap==0) {
    nap = 1 ;
    for (int i=0 ; i<MAX_BASE ; i++) {
      ope_ptens_rr_mat[i] = _ope_ptens_rr_mat_pas_prevu ;
    }
    // Les routines existantes
    ope_ptens_rr_mat[R_CHEB >> TRA_R] = _ope_ptens_rr_mat_r_cheb ;
    ope_ptens_rr_mat[R_CHEBU >> TRA_R] = _ope_ptens_rr_mat_r_chebu ;
    ope_ptens_rr_mat[R_CHEBP >> TRA_R] = _ope_ptens_rr_mat_r_chebp ;
    ope_ptens_rr_mat[R_CHEBI >> TRA_R] = _ope_ptens_rr_mat_r_chebi ;
  }
  assert (l>1) ;
  Matrice res(ope_ptens_rr_mat[base_r](n, l, echelle, puis)) ;
  return res ;
}

}
