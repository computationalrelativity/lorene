 /*
  *   Copyright (c) 1999-2001 Jerome Novak
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
  * Revision 1.8  2023/05/24 09:53:46  g_servignat
  * Added multiplication by \xi for R_CHEB base
  *
  * Revision 1.7  2016/12/05 16:18:08  j_novak
  * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
  *
  * Revision 1.6  2015/03/05 08:49:32  j_novak
  * Implemented operators with Legendre bases.
  *
  * Revision 1.5  2014/10/13 08:53:25  j_novak
  * Lorene classes and functions now belong to the namespace Lorene.
  *
  * Revision 1.4  2008/08/27 11:22:25  j_novak
  * Minor modifications
  *
  * Revision 1.3  2008/08/27 08:50:10  jl_cornou
  * Added Jacobi(0,2) polynomials
  *
  * Revision 1.2  2004/11/23 15:16:01  m_forot
  *
  * Added the bases for the cases without any equatorial symmetry
  *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
  *
  * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
  * LORENE
  *
  * Revision 1.3  2000/09/07  12:49:53  phil
  * *** empty log message ***
  *
  * Revision 1.2  2000/02/24  16:42:18  eric
  * Initialisation a zero du tableau xo avant le calcul.
  *
  * Revision 1.1  1999/11/16  13:37:41  novak
  * Initial revision
  *
  *
  * $Header$
  *
  */
 
 /* 
  * Ensemble des routines de base de multiplication par x
  * (Utilisation interne)
  * 
  *  void _mult_x_XXXX(Tbl * t, int & b)
  *  t   pointeur sur le Tbl d'entree-sortie
  *  b   base spectrale
  * 
  */
  
  // Fichier includes
 #include "tbl.h"
 
 
         //-----------------------------------
         // Routine pour les cas non prevus --
         //-----------------------------------
 
 namespace Lorene {
 void _mult_x_pas_prevu(Tbl * tb, int& base) {
     cout << "mult_x pas prevu..." << endl ;
     cout << "Tbl: " << tb << " base: " << base << endl ;
     abort () ;
     exit(-1) ;
 }
 
             //-------------
             // Identite ---
             //-------------
 
 void _mult_x_identite(Tbl* , int& ) {
     return ;
 }
 
             //---------------
             // cas R_CHEBP --
             //--------------
 
 void _mult_x_r_chebp(Tbl* tb, int& base)
     {
     // Peut-etre rien a faire ?
     if (tb->get_etat() == ETATZERO) {
     int base_t = base & MSQ_T ;
     int base_p = base & MSQ_P ;
     base = base_p | base_t | R_CHEBI ;
     return ;
     }
     
     // Pour le confort
     int nr = (tb->dim).dim[0] ;     // Nombre
     int nt = (tb->dim).dim[1] ;     //   de points
     int np = (tb->dim).dim[2] ;     //      physiques REELS
     np = np - 2 ;           // Nombre de points physiques
     
     // pt. sur le tableau de double resultat
     double* xo = new double [tb->get_taille()];
     
     // Initialisation a zero :
     for (int i=0; i<tb->get_taille(); i++) {
     xo[i] = 0 ; 
     }
     
     // On y va...
     double* xi = tb->t ;
     double* xci = xi ;  // Pointeurs
     double* xco = xo ;  //  courants
 
     int borne_phi = np + 1 ; 
     if (np == 1) {
     borne_phi = 1 ; 
     }
     
     for (int k=0 ; k< borne_phi ; k++)
     if (k==1) {
         xci += nr*nt ;
         xco += nr*nt ;
     }
     else {
     for (int j=0 ; j<nt ; j++) {
 
         xco[0] = xci[0] + 0.5*xci[1] ;
         for (int i = 1 ; i < nr-1 ; i++ ) {
         xco[i] = 0.5*(xci[i]+xci[i+1]) ;
         }   // Fin de la boucle sur r
         xco[nr-1] = 0 ;
         
         xci += nr ;
         xco += nr ;
     }   // Fin de la boucle sur theta
     }   // Fin de la boucle sur phi
     
     // On remet les choses la ou il faut
     delete [] tb->t ;
     tb->t = xo ;
     
     // base de developpement
     // pair -> impair
     int base_t = base & MSQ_T ;
     int base_p = base & MSQ_P ;
     base = base_p | base_t | R_CHEBI ;
 
 }
 
             //----------------
             // cas R_CHEBI ---
             //----------------
 
 void _mult_x_r_chebi(Tbl* tb, int& base)
 {
 
     // Peut-etre rien a faire ?
     if (tb->get_etat() == ETATZERO) {
     int base_t = base & MSQ_T ;
     int base_p = base & MSQ_P ;
     base = base_p | base_t | R_CHEBP ;
     return ;
     }
     
     // Pour le confort
     int nr = (tb->dim).dim[0] ;     // Nombre
     int nt = (tb->dim).dim[1] ;     //   de points
     int np = (tb->dim).dim[2] ;     //      physiques REELS
     np = np - 2 ;           // Nombre de points physiques
     
     // pt. sur le tableau de double resultat
     double* xo = new double [tb->get_taille()];
     
     // Initialisation a zero :
     for (int i=0; i<tb->get_taille(); i++) {
     xo[i] = 0 ; 
     }
     
     // On y va...
     double* xi = tb->t ;
     double* xci = xi ;  // Pointeurs
     double* xco = xo ;  //  courants
     
     int borne_phi = np + 1 ; 
     if (np == 1) {
     borne_phi = 1 ; 
     }
     
     for (int k=0 ; k< borne_phi ; k++) 
     if (k == 1)  {
         xci += nr*nt ;
         xco += nr*nt ;
         }
     else {
     for (int j=0 ; j<nt ; j++) {
         
         xco[0] = 0.5*xci[0] ;
         for (int i = 1 ; i < nr-1 ; i++ ) {
         xco[i] = (xci[i]+xci[i-1])*0.5 ;
         }   // Fin de la premiere boucle sur r
         xco[nr-1] = 0.5*xci[nr-2] ;
         
         xci += nr ;
         xco += nr ;
     }   // Fin de la boucle sur theta
     }   // Fin de la boucle sur phi
     
     // On remet les choses la ou il faut
     delete [] tb->t ;
     tb->t = xo ;
     
     // base de developpement
     // impair -> pair
     int base_t = base & MSQ_T ;
     int base_p = base & MSQ_P ;    
     base = base_p | base_t | R_CHEBP ;
 }
 
             //--------------------
             // cas R_CHEBPIM_P --
             //------------------
 
 void _mult_x_r_chebpim_p(Tbl* tb, int& base)
 {
   
     // Peut-etre rien a faire ?
     if (tb->get_etat() == ETATZERO) {
     int base_t = base & MSQ_T ;
     int base_p = base & MSQ_P ;
     base = base_p | base_t | R_CHEBPIM_I ;
     return ;
     }
     
     // Pour le confort
     int nr = (tb->dim).dim[0] ;     // Nombre
     int nt = (tb->dim).dim[1] ;     //   de points
     int np = (tb->dim).dim[2] ;     //      physiques REELS
     np = np - 2 ;          
    
     // pt. sur le tableau de double resultat
     double* xo = new double [tb->get_taille()];
     
     // Initialisation a zero :
     for (int i=0; i<tb->get_taille(); i++) {
     xo[i] = 0 ; 
     }
        
     
     // On y va...
     double* xi = tb->t ;
     double* xci ;   // Pointeurs
     double* xco ;   //  courants
     
     
     // Partie paire
     xci = xi ;
     xco = xo ;
   
     int auxiliaire ;
     
     for (int k=0 ; k<np+1 ; k += 4)  {
     
     auxiliaire = (k==np) ? 1 : 2 ; // evite de prendre le dernier coef
     for (int kmod=0 ; kmod<auxiliaire ; kmod++) {
         
         
         // On evite le sin (0*phi)
     
     if ((k==0) && (kmod == 1)) { 
         xco += nr*nt ;
         xci += nr*nt ;
         }
         
     else 
         for (int j=0 ; j<nt ; j++) {
 
           xco[0] = xci[0] + 0.5*xci[1] ;
           for (int i = 1 ; i < nr-1 ; i++ ) {
             xco[i] = 0.5*(xci[i]+xci[i+1]) ;
           } // Fin de la boucle sur r
           xco[nr-1] = 0 ;
 
           xci += nr ; // au
           xco += nr ; //    suivant
         }   // Fin de la boucle sur theta
     }   // Fin de la boucle sur kmod
     xci += 2*nr*nt ;    // saute
     xco += 2*nr*nt ;    //  2 phis
     }   // Fin de la boucle sur phi
 
     // Partie impaire
     xci = xi + 2*nr*nt ;
     xco = xo + 2*nr*nt ;
     for (int k=2 ; k<np+1 ; k += 4) {
     
     auxiliaire = (k==np) ? 1 : 2 ; // evite de prendre le dernier coef
     for (int kmod=0 ; kmod<auxiliaire ; kmod++) {
         for (int j=0 ; j<nt ; j++) {
 
           xco[0] = 0.5*xci[0] ;
           for (int i = 1 ; i < nr-1 ; i++ ) {
         xco[i] = (xci[i]+xci[i-1])*0.5 ;
           } // Fin de la premiere boucle sur r
           xco[nr-1] = 0.5*xci[nr-2] ;
         
           xci += nr ;
           xco += nr ;
         }   // Fin de la boucle sur theta
     }   // Fin de la boucle sur kmod
     xci += 2*nr*nt ;
     xco += 2*nr*nt ;
     }   // Fin de la boucle sur phi
     
     // On remet les choses la ou il faut
     delete [] tb->t ;
     tb->t = xo ;
     
     // base de developpement
     // (pair,impair) -> (impair,pair)
     int base_t = base & MSQ_T ;
     int base_p = base & MSQ_P ;
     base = base_p | base_t | R_CHEBPIM_I ;
 }
 
             //-------------------
             // cas R_CHEBPIM_I --
             //-------------------
 
 void _mult_x_r_chebpim_i(Tbl* tb, int& base)
 {
 
    // Peut-etre rien a faire ?
     if (tb->get_etat() == ETATZERO) {
     int base_t = base & MSQ_T ;
     int base_p = base & MSQ_P ;
     base = base_p | base_t | R_CHEBPIM_P ;
     return ;
     }
     
     // Pour le confort
     int nr = (tb->dim).dim[0] ;     // Nombre
     int nt = (tb->dim).dim[1] ;     //   de points
     int np = (tb->dim).dim[2] ;     //      physiques REELS
     np = np - 2 ;          
    
     // pt. sur le tableau de double resultat
     double* xo = new double [tb->get_taille()];
     
     // Initialisation a zero :
     for (int i=0; i<tb->get_taille(); i++) {
     xo[i] = 0 ; 
     }
     
     // On y va...
     double* xi = tb->t ;
     double* xci ;   // Pointeurs
     double* xco ;   //  courants
     
     // Partie impaire
     xci = xi ;
     xco = xo ;
     
       
     int auxiliaire ;
     
     for (int k=0 ; k<np+1 ; k += 4)  {
     
     auxiliaire = (k==np) ? 1 : 2 ; // evite de prendre le dernier coef
     for (int kmod=0 ; kmod<auxiliaire ; kmod++) {
         
         
         // On evite le sin (0*phi)
     
     if ((k==0) && (kmod == 1)) { 
         xco += nr*nt ;
         xci += nr*nt ;
         }
         
     else 
       for (int j=0 ; j<nt ; j++) {
 
         xco[0] = 0.5*xci[0] ;
         for (int i = 1 ; i < nr-1 ; i++ ) {
           xco[i] = (xci[i]+xci[i-1])*0.5 ;
         }   // Fin de la premiere boucle sur r
         xco[nr-1] = 0.5*xci[nr-2] ;
         
         xci += nr ;
         xco += nr ;
       }   // Fin de la boucle sur theta
     }   // Fin de la boucle sur kmod
     xci += 2*nr*nt ;
     xco += 2*nr*nt ;
     }   // Fin de la boucle sur phi
 
     // Partie paire
     xci = xi + 2*nr*nt ;
     xco = xo + 2*nr*nt ;
     for (int k=2 ; k<np+1 ; k += 4) {
     
       auxiliaire = (k==np) ? 1 : 2 ; // evite de prendre le dernier coef
       for (int kmod=0 ; kmod<auxiliaire ; kmod++) {
     for (int j=0 ; j<nt ; j++) {
       
       xco[0] = xci[0] + 0.5*xci[1] ;
       for (int i = 1 ; i < nr-1 ; i++ ) {
         xco[i] = 0.5*(xci[i]+xci[i+1]) ;
       } // Fin de la boucle sur r
       xco[nr-1] = 0 ;
       
       xci += nr ;
       xco += nr ;
     }   // Fin de la boucle sur theta
       }   // Fin de la boucle sur kmod
       xci += 2*nr*nt ;
       xco += 2*nr*nt ;
     }   // Fin de la boucle sur phi
     
     // On remet les choses la ou il faut
     delete [] tb->t ;
     tb->t = xo ;
     
     // base de developpement
     // (impair,pair) -> (pair,impair)
     int base_t = base & MSQ_T ;
     int base_p = base & MSQ_P ;
     base = base_p | base_t | R_CHEBPIM_P ;
 }
 
             //---------------
             // cas R_CHEBPI_P --
             //--------------
 
 void _mult_x_r_chebpi_p(Tbl* tb, int& base)
     {
     // Peut-etre rien a faire ?
     if (tb->get_etat() == ETATZERO) {
     int base_t = base & MSQ_T ;
     int base_p = base & MSQ_P ;
     base = base_p | base_t | R_CHEBPI_I ;
     return ;
     }
     
     // Pour le confort
     int nr = (tb->dim).dim[0] ;     // Nombre
     int nt = (tb->dim).dim[1] ;     //   de points
     int np = (tb->dim).dim[2] ;     //      physiques REELS
     np = np - 2 ;           // Nombre de points physiques
     
     // pt. sur le tableau de double resultat
     double* xo = new double [tb->get_taille()];
     
     // Initialisation a zero :
     for (int i=0; i<tb->get_taille(); i++) {
     xo[i] = 0 ; 
     }
     
     // On y va...
     double* xi = tb->t ;
     double* xci = xi ;  // Pointeurs
     double* xco = xo ;  //  courants
 
     int borne_phi = np + 1 ; 
     if (np == 1) {
     borne_phi = 1 ; 
     }
     
     for (int k=0 ; k< borne_phi ; k++)
     if (k==1) {
         xci += nr*nt ;
         xco += nr*nt ;
     }
     else {
     for (int j=0 ; j<nt ; j++) {
         int  l = j%2 ;
         if(l==0){
           xco[0] = xci[0] + 0.5*xci[1] ;
           for (int i = 1 ; i < nr-1 ; i++ ) {
         xco[i] = 0.5*(xci[i]+xci[i+1]) ;
           } // Fin de la boucle sur r
           xco[nr-1] = 0 ;
         } else {
           xco[0] = 0.5*xci[0] ;
           for (int i = 1 ; i < nr-1 ; i++ ) {
         xco[i] = (xci[i]+xci[i-1])*0.5 ;
           } // Fin de la premiere boucle sur r
           xco[nr-1] = 0.5*xci[nr-2] ;
         }
         xci += nr ;
         xco += nr ;
     }   // Fin de la boucle sur theta
     }   // Fin de la boucle sur phi
     
     // On remet les choses la ou il faut
     delete [] tb->t ;
     tb->t = xo ;
     
     // base de developpement
     // pair -> impair
     int base_t = base & MSQ_T ;
     int base_p = base & MSQ_P ;
     base = base_p | base_t | R_CHEBPI_I ;
 
 }
 
             //----------------
             // cas R_CHEBPI_I ---
             //----------------
 
 void _mult_x_r_chebpi_i(Tbl* tb, int& base)
 {
 
     // Peut-etre rien a faire ?
     if (tb->get_etat() == ETATZERO) {
     int base_t = base & MSQ_T ;
     int base_p = base & MSQ_P ;
     base = base_p | base_t | R_CHEBPI_P ;
     return ;
     }
     
     // Pour le confort
     int nr = (tb->dim).dim[0] ;     // Nombre
     int nt = (tb->dim).dim[1] ;     //   de points
     int np = (tb->dim).dim[2] ;     //      physiques REELS
     np = np - 2 ;           // Nombre de points physiques
     
     // pt. sur le tableau de double resultat
     double* xo = new double [tb->get_taille()];
     
     // Initialisation a zero :
     for (int i=0; i<tb->get_taille(); i++) {
     xo[i] = 0 ; 
     }
     
     // On y va...
     double* xi = tb->t ;
     double* xci = xi ;  // Pointeurs
     double* xco = xo ;  //  courants
     
     int borne_phi = np + 1 ; 
     if (np == 1) {
     borne_phi = 1 ; 
     }
     
     for (int k=0 ; k< borne_phi ; k++) 
     if (k == 1)  {
         xci += nr*nt ;
         xco += nr*nt ;
         }
     else {
     for (int j=0 ; j<nt ; j++) {
         int  l = j%2 ;
         if(l==1){
           xco[0] = xci[0] + 0.5*xci[1] ;
           for (int i = 1 ; i < nr-1 ; i++ ) {
         xco[i] = 0.5*(xci[i]+xci[i+1]) ;
           } // Fin de la boucle sur r
           xco[nr-1] = 0 ;
         } else {
           xco[0] = 0.5*xci[0] ;
           for (int i = 1 ; i < nr-1 ; i++ ) {
         xco[i] = (xci[i]+xci[i-1])*0.5 ;
           } // Fin de la premiere boucle sur r
           xco[nr-1] = 0.5*xci[nr-2] ;
         }
         xci += nr ;
         xco += nr ;
     }   // Fin de la boucle sur theta
     }   // Fin de la boucle sur phi
     
     // On remet les choses la ou il faut
     delete [] tb->t ;
     tb->t = xo ;
     
     // base de developpement
     // impair -> pair
     int base_t = base & MSQ_T ;
     int base_p = base & MSQ_P ;    
     base = base_p | base_t | R_CHEBPI_P ;
 }
 
             //---------------
             // cas R_JACO02 -
             //---------------
 
 void _mult_x_r_jaco02(Tbl* tb, int&)
     {
     // Peut-etre rien a faire ?
     if (tb->get_etat() == ETATZERO) {
     return ;
     }
     
     // Pour le confort
     int nr = (tb->dim).dim[0] ;     // Nombre
     int nt = (tb->dim).dim[1] ;     //   de points
     int np = (tb->dim).dim[2] ;     //      physiques REELS
     np = np - 2 ;           // Nombre de points physiques
     
     // pt. sur le tableau de double resultat
     double* xo = new double [tb->get_taille()];
     
     // Initialisation a zero :
     for (int i=0; i<tb->get_taille(); i++) {
     xo[i] = 0 ; 
     }
     
     // On y va...
     double* xi = tb->t ;
     double* xci = xi ;  // Pointeurs
     double* xco = xo ;  //  courants
 
     int borne_phi = np + 1 ; 
     if (np == 1) {
     borne_phi = 1 ; 
     }
     
     for (int k=0 ; k< borne_phi ; k++)
     if (k==1) {
         xci += nr*nt ;
         xco += nr*nt ;
     }
     else {
     for (int j=0 ; j<nt ; j++) {
 
             xco[0] = 1.5*xci[0] + 0.3*xci[1] ;
             for (int i = 1 ; i < nr-1 ; i++) {
         xco[i] = i*(i+2)/double((i+1)*(2*i+1))*xci[i-1] + (i*i+3*i+3)/double((i+1)*(i+2))*xci[i] + (i+1)*(i+3)/double((i+2)*(2*i+5))*xci[i+1] ;
             }
             xco[nr-1] = (nr*nr-1)/double((nr)*(2*nr-1))*xci[nr-2] + (1+1/double((nr)*(nr+1)))*xci[nr-1] ;
         
         xci += nr ;
         xco += nr ;
     }   // Fin de la boucle sur theta
     }   // Fin de la boucle sur phi
     
     // On remet les choses la ou il faut
     delete [] tb->t ;
     tb->t = xo ;
     
     // base de developpement
     // inchangée
 
 }
 
             //--------------
             // cas R_LEGP --
             //--------------
 
 void _mult_x_r_legp(Tbl* tb, int& base)
     {
     // Peut-etre rien a faire ?
     if (tb->get_etat() == ETATZERO) {
     int base_t = base & MSQ_T ;
     int base_p = base & MSQ_P ;
     base = base_p | base_t | R_LEGI ;
     return ;
     }
     
     // Pour le confort
     int nr = (tb->dim).dim[0] ;     // Nombre
     int nt = (tb->dim).dim[1] ;     //   de points
     int np = (tb->dim).dim[2] ;     //      physiques REELS
     np = np - 2 ;           // Nombre de points physiques
     
     // pt. sur le tableau de double resultat
     double* xo = new double [tb->get_taille()];
     
     // Initialisation a zero :
     for (int i=0; i<tb->get_taille(); i++) {
     xo[i] = 0 ; 
     }
     
     // On y va...
     double* xi = tb->t ;
     double* xci = xi ;  // Pointeurs
     double* xco = xo ;  //  courants
 
     int borne_phi = np + 1 ; 
     if (np == 1) {
     borne_phi = 1 ; 
     }
     
     for (int k=0 ; k< borne_phi ; k++)
     if (k==1) {
         xci += nr*nt ;
         xco += nr*nt ;
     }
     else {
     for (int j=0 ; j<nt ; j++) {
 
         xco[0] = xci[0] + 0.4*xci[1] ;
         for (int i = 1 ; i < nr-1 ; i++ ) {
           xco[i] = double(2*i+1)*xci[i]/double(4*i+1)
         + double(2*i+2)*xci[i+1]/double(4*i+5) ;
         }   // Fin de la boucle sur r
         xco[nr-1] = 0 ;
         
         xci += nr ;
         xco += nr ;
     }   // Fin de la boucle sur theta
     }   // Fin de la boucle sur phi
     
     // On remet les choses la ou il faut
     delete [] tb->t ;
     tb->t = xo ;
     
     // base de developpement
     // pair -> impair
     int base_t = base & MSQ_T ;
     int base_p = base & MSQ_P ;
     base = base_p | base_t | R_LEGI ;
 
 }
 
             //----------------
             // cas R_LEGI ---
             //----------------
 
 void _mult_x_r_legi(Tbl* tb, int& base)
 {
 
     // Peut-etre rien a faire ?
     if (tb->get_etat() == ETATZERO) {
     int base_t = base & MSQ_T ;
     int base_p = base & MSQ_P ;
     base = base_p | base_t | R_LEGP ;
     return ;
     }
     
     // Pour le confort
     int nr = (tb->dim).dim[0] ;     // Nombre
     int nt = (tb->dim).dim[1] ;     //   de points
     int np = (tb->dim).dim[2] ;     //      physiques REELS
     np = np - 2 ;           // Nombre de points physiques
     
     // pt. sur le tableau de double resultat
     double* xo = new double [tb->get_taille()];
     
     // Initialisation a zero :
     for (int i=0; i<tb->get_taille(); i++) {
     xo[i] = 0 ; 
     }
     
     // On y va...
     double* xi = tb->t ;
     double* xci = xi ;  // Pointeurs
     double* xco = xo ;  //  courants
     
     int borne_phi = np + 1 ; 
     if (np == 1) {
     borne_phi = 1 ; 
     }
     
     for (int k=0 ; k< borne_phi ; k++) 
     if (k == 1)  {
         xci += nr*nt ;
         xco += nr*nt ;
         }
     else {
     for (int j=0 ; j<nt ; j++) {
         
         xco[0] = xci[0]/3. ;
         for (int i = 1 ; i < nr-1 ; i++ ) {
           xco[i] = double(2*i+1)*xci[i]/double(4*i+3)
         + double(2*i)*xci[i-1]/double(4*i-1) ;
         }   // Fin de la premiere boucle sur r
         xco[nr-1] = double(2*nr-2)*xci[nr-2]/double(4*nr-5) ;
         
         xci += nr ;
         xco += nr ;
     }   // Fin de la boucle sur theta
     }   // Fin de la boucle sur phi
     
     // On remet les choses la ou il faut
     delete [] tb->t ;
     tb->t = xo ;
     
     // base de developpement
     // impair -> pair
     int base_t = base & MSQ_T ;
     int base_p = base & MSQ_P ;    
     base = base_p | base_t | R_LEGP ;
 }
 
            //--------------
            // cas R_CHEB --
            //--------------

void _mult_x_r_cheb(Tbl* tb, int& base)
     {
     // Peut-etre rien a faire ?
     if (tb->get_etat() == ETATZERO) {
     int base_t = base & MSQ_T ;
     int base_p = base & MSQ_P ;
     base = base_p | base_t | R_CHEB ;
     return ;
     }
     
     // Pour le confort
     int nr = (tb->dim).dim[0] ;     // Nombre
     int nt = (tb->dim).dim[1] ;     //   de points
     int np = (tb->dim).dim[2] ;     //      physiques REELS
     np = np - 2 ;           // Nombre de points physiques
     
     // pt. sur le tableau de double resultat
     double* xo = new double [tb->get_taille()];
     
     // Initialisation a zero :
     for (int i=0; i<tb->get_taille(); i++) {
     xo[i] = 0 ; 
     }
     
     // On y va...
     double* xi = tb->t ;
     double* xci = xi ;  // Pointeurs
     double* xco = xo ;  //  courants
 
     int borne_phi = np + 1 ; 
     if (np == 1) {
     borne_phi = 1 ; 
     }
     
     for (int k=0 ; k< borne_phi ; k++)
     if (k==1) {
         xci += nr*nt ;
         xco += nr*nt ;
     }
     else {
     for (int j=0 ; j<nt ; j++) {
 
         xco[0] = .5*xci[1] ;
         xco[1] = xci[0] + 0.5*xci[2] ;
         for (int i = 2 ; i < nr-1 ; i++ ) {
         xco[i] = 0.5*(xci[i-1]+xci[i+1]) ;
         }   // Fin de la boucle sur r
         xco[nr-1] = .5*xci[nr-2] ;
         
         xci += nr ;
         xco += nr ;
     }   // Fin de la boucle sur theta
     }   // Fin de la boucle sur phi
     
     // On remet les choses la ou il faut
     delete [] tb->t ;
     tb->t = xo ;
     
     // base de developpement
     // pair -> impair
     int base_t = base & MSQ_T ;
     int base_p = base & MSQ_P ;
     base = base_p | base_t | R_CHEB ;
 
 }
 
 
 
 }