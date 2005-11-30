/*
 *  Methods of class Map_et
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


char map_et_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.11  2005/11/30 11:09:07  p_grandclement
 * Changes for the Bin_ns_bh project
 *
 * Revision 1.10  2004/03/25 10:29:23  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.9  2004/01/29 08:50:03  p_grandclement
 * Modification of Map::operator==(const Map&) and addition of the surface
 * integrales using Scalar.
 *
 * Revision 1.8  2003/10/15 10:36:52  e_gourgoulhon
 * In method fait_poly(): changed local variable name x to x1, not to shadow
 *  Coord's x.
 *
 * Revision 1.7  2003/07/07 20:01:43  p_grandclement
 * change assert in constructor for map_et from a surface
 *
 * Revision 1.6  2003/06/04 21:11:55  p_grandclement
 * Correction of separation in odd-even harmonics
 *
 * Revision 1.5  2002/10/16 14:36:41  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.4  2002/05/07 07:10:44  e_gourgoulhon
 * Compatibilty with xlC compiler on IBM SP2:
 * 	suppressed the parenthesis around argument of instruction new:
 * 	e.g.   aa = new (Tbl*[nzone])  --->  aa = new Tbl*[nzone]
 * 		result = new (Param)   --->  result = new Param
 *
 * Revision 1.3  2002/01/15 15:53:06  p_grandclement
 * I have had a constructor fot map_et using the equation of the surface
 * of the star.
 *
 * Revision 1.2  2001/12/04 21:27:53  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.11  2001/02/28  11:04:20  eric
 * 1ere version testee de resize.
 *
 * Revision 1.10  2001/02/26  17:29:42  eric
 * Ajout de la fonction resize.
 *
 * Revision 1.9  2000/08/18  11:10:48  eric
 * Ajout de l'operateur d'affectation a un autre Map_et.
 *
 * Revision 1.8  2000/01/24  16:42:36  eric
 * Ajout de la fonction virtuelle operator=(const Map_af& ).
 *
 * Revision 1.7  2000/01/24  11:03:28  eric
 * Correction d'une erreur dans le constructeur par lecture de fichier:
 *  ff et gg doivent etre construits sur mgi.get_angu() et non sur mgi.
 *
 * Revision 1.6  1999/12/20  10:24:49  eric
 * Ajout des fonctions de lecture des parametres de Map_et:
 *   get_alpha(), get_beta(), get_ff(), get_gg().
 *
 * Revision 1.5  1999/12/17  11:20:08  eric
 * Ajout de la fonction homothetie.
 *
 * Revision 1.4  1999/12/17  09:14:30  eric
 * Amelioration de l'affichage.
 *
 * Revision 1.3  1999/11/24  16:31:41  eric
 * Ajout des fonctions set_ff et set_gg.
 *
 * Revision 1.2  1999/11/24  11:22:44  eric
 * Map_et : fonctions de constructions amies.
 *
 * Revision 1.1  1999/11/22  10:37:36  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// headers C
#include <math.h>

// headers Lorene
#include "proto.h"
#include "map.h"
#include "utilitaires.h"
#include "unites.h"

			//--------------//
			// Constructors //
			//--------------//

// -----------------------
// Constructor from a grid
// -----------------------
Map_et::Map_et(const Mg3d& mgrille, const double* bornes) 
	      : Map_radial(mgrille),
	        aasx( mgrille.get_nr(0) ), 
	        aasx2( mgrille.get_nr(0) ), 
	        zaasx( mgrille.get_nr(mgrille.get_nzone()-1) ), 
	        zaasx2( mgrille.get_nr(mgrille.get_nzone()-1) ), 
	        bbsx( mgrille.get_nr(0) ), 
	        bbsx2( mgrille.get_nr(0) ), 
		ff(mgrille.get_angu()) , 
		gg(mgrille.get_angu())   
{
    // The Coord rsxdxdr and rsx2drdx are constructed by the default Coord 
    // constructor
    
    // Assignement of the building functions of the Coord's
    // ----------------------------------------------------
    set_coord() ; 
    
    
    // alpha and beta
    // --------------
    int nzone = mg->get_nzone() ;
    
    alpha = new double[nzone] ;
    beta = new double[nzone] ;

    for (int l=0 ; l<nzone ; l++) {
	switch (mg->get_type_r(l)) {
	    case RARE:	{
		alpha[l] = bornes[l+1] - bornes[l] ;
		beta[l] = bornes[l] ;
		break ; 
	    }
	    
	    case FIN:	{
		alpha[l] = (bornes[l+1] - bornes[l]) * .5 ;
		beta[l] = (bornes[l+1] + bornes[l]) * .5 ;
		break ;
	    }
	    
	    case UNSURR: {
		double umax = double(1) / bornes[l] ;
		double umin = double(1) /bornes[l+1] ;
		alpha[l] = (umin - umax) * double(0.5) ;  // u est une fonction decroissante
		beta[l] = (umin + umax) * double(0.5) ;   //  de l'indice i en r
		break ;
	    }
	    
	    default:	{
		cout << "Map_et::Map_et: unkown type_r ! " << endl ;
		abort () ;
		break ;
	    }
	    
	}
    }	    // End of the loop onto the domains


    // Radial polynomials A(x) and B(x)
    // --------------------------------
    
    fait_poly() ; 

    // Initialisation at zero of the functions F(theta',phi') and G(theta',phi')
    // -------------------------------------------------------------------------

    ff.set_etat_zero() ;
    gg.set_etat_zero() ;
    
    ff.std_base_scal() ;    // Standard spectral bases for F 
    gg.std_base_scal() ;    // Standard spectral bases for G

}

Map_et::Map_et(const Mg3d& grille, const double* r_lim, const Tbl& S_0) : 
                Map_radial(grille),
	        aasx(grille.get_nr(0) ), 
	        aasx2(grille.get_nr(0) ), 
	        zaasx(grille.get_nr(grille.get_nzone()-1) ), 
	        zaasx2(grille.get_nr(grille.get_nzone()-1) ), 
	        bbsx(grille.get_nr(0) ), 
	        bbsx2(grille.get_nr(0) ), 
		ff(grille.get_angu()) , 
		gg(grille.get_angu()) {
  
  assert (S_0.get_ndim() == 2) ;
  assert (S_0.get_dim(0) == grille.get_nt(0)) ;
  assert (S_0.get_dim(1) == grille.get_np(0)) ;

  Map_et mapping (grille, r_lim) ;

  int nz = grille.get_nzone() ;
  assert (nz >2) ;

  // Le noyau :
  int np = grille.get_np(0) ;
  int nt = grille.get_nt(0) ;
  
  double * cf = new double [nt*(np+2)] ;
  for (int k=0 ; k<np ; k++)
    for (int j=0 ; j<nt ; j++)
      cf[k*nt+j] = S_0(k,j) - S_0(0,0) ;

  int* deg = new int [3] ;
  deg[0] = np ;
  deg[1] = nt ;
  deg[2] = 1 ;

  int* dim = new int [3] ;
  dim[0] = np+2 ;
  dim[1] = nt ;
  dim[2] = 1 ;
  
  Tbl ff_nucleus (np,nt) ;
  ff_nucleus.set_etat_qcq() ;

  Tbl gg_nucleus (np,nt) ;
  gg_nucleus.set_etat_qcq() ;
  
  // On recupere la base en phi :
  int base_p = grille.std_base_scal().get_base_p(0) ;
  // Selon les cas (pas tres propre mais bon ...)
  double * odd ; 
  double * even ;
  double * coloc_odd ;
  double * coloc_even ;

  switch (base_p) {
  case P_COSSIN:
     cfpcossin (deg,dim,cf) ;

     // Separation des harmoniques paires et impaires :
     odd = new double [nt*(np+2)] ;
     even = new double [nt*(np+2)] ;


     for (int k=0 ; k<np+2 ; k++)
       if ((k%4 == 0) || (k%4==1))
	 for (int j=0 ; j<nt ; j++) {
	   odd[k*nt+j] = 0 ;
	   even[k*nt+j] = cf[k*nt+j] ;
	 }
       else
	 if ((k%4 == 2) || (k%4 == 3))
	 for (int j=0 ; j<nt ; j++) {
	   even[k*nt+j] = 0 ;
	   odd[k*nt+j] = cf[k*nt+j] ;
	 }

	 else {
	   cout << "Erreur bizzare..." << endl ;
	   abort() ;
	 }

     coloc_odd = new double [nt*np] ;
     coloc_even = new double [nt*np] ;

     cipcossin (deg,dim,deg,odd,coloc_odd) ;
     cipcossin (deg,dim,deg,even,coloc_even) ;
     for (int k=0 ; k<np ; k++)
       for (int j=0 ; j<nt ; j++) {
	 gg_nucleus.set(k,j) = coloc_even[k*nt+j] ;
	 ff_nucleus.set(k,j) = coloc_odd[k*nt+j] ;
       }

     delete [] even ;
     delete [] odd ;
     delete [] coloc_even ;
     delete [] coloc_odd ;
     delete[] dim ;
     delete [] deg ;
     delete [] cf ;
     
     break;
  default:
    cout << "Base_p != P_COSSIN not implemented in Map_et constructor" << 
      endl ;
    abort() ;
  }

  double mu_nucleus = - min(gg_nucleus) ;
  double alpha_nucleus = S_0(0,0)-mu_nucleus ;

  ff_nucleus /= alpha_nucleus ;
  gg_nucleus += mu_nucleus ;
  gg_nucleus /= alpha_nucleus ;
  
  // First shell : much simpler no ?
  Tbl ff_shell (np,nt) ;
  ff_shell.set_etat_qcq() ;
  ff_shell = S_0 - S_0(0,0) ;
  
  double lambda_shell = -max(ff_shell) ;
  
  double R_ext = r_lim[2] ;

  double beta_shell = (R_ext+S_0(0,0)-lambda_shell)/2. ;
  double alpha_shell = (R_ext-S_0(0,0)+lambda_shell)/2. ;
  
  ff_shell += lambda_shell ;
  ff_shell /= alpha_shell ;

  ff.annule_hard() ;
  gg.annule_hard() ;
  
  ff.set_etat_c_qcq() ;
  gg.set_etat_c_qcq() ;

  for (int k=0 ; k<np ; k++)
    for (int j=0 ; j<nt ; j++) {
      ff.set(0,k,j,0) = ff_nucleus(k,j) ;
      gg.set(0,k,j,0) = gg_nucleus(k,j) ;
      ff.set(1,k,j,0) = ff_shell(k,j) ;
    }

  gg.annule(1,nz-1) ;
  ff.annule(2,nz-1) ;
  
  ff.std_base_scal() ;
  gg.std_base_scal() ;
  
  alpha = new double[nz] ;
  alpha[0] = alpha_nucleus ;
  alpha[1] = alpha_shell ;
  
  beta = new double[nz] ;
  beta[0] = 0 ;
  beta[1] = beta_shell ;
  for (int i=2 ; i<nz ; i++) {
    alpha[i] = mapping.get_alpha()[i] ;
    beta[i] = mapping.get_beta()[i] ;
  }

  fait_poly() ;
  set_coord() ;
}
// ------------------
// Copy constructor 
// ------------------
Map_et::Map_et(const Map_et& mpi) : Map_radial(mpi) , 
	        aasx( mpi.aasx ), 
	        aasx2( mpi.aasx2 ), 
	        zaasx( mpi.zaasx ), 
	        zaasx2( mpi.zaasx2 ), 
	        bbsx( mpi.bbsx ), 
	        bbsx2( mpi.bbsx2 ), 
		ff(mpi.ff) , 
		gg(mpi.gg)   
{
    // Assignement of the building functions of the Coord's
    // ----------------------------------------------------
    set_coord() ; 

    // alpha and beta
    // --------------
    int nzone = mg->get_nzone() ;
    
    alpha = new double[nzone] ;
    beta = new double[nzone] ;

    for (int l=0 ; l<nzone ; l++) {
	alpha[l] = mpi.alpha[l] ;
	beta[l] = mpi.beta[l] ;
    }
    
    // Radial polynomials A(x) and B(x)
    // --------------------------------
    
    fait_poly() ; 

} 
			//------------------------------------------//
			//  Modification of the mapping parameters  //
			//------------------------------------------//

void Map_et::set_alpha(double alpha0, int l) {
    
    assert(l>=0) ; 
    assert(l<mg->get_nzone()) ; 
    
    alpha[l] = alpha0 ; 
    
    reset_coord() ; 
    
}

void Map_et::set_beta(double beta0, int l) {
    
    assert(l>=0) ; 
    assert(l<mg->get_nzone()) ; 
    
    beta[l] = beta0 ; 
    
    reset_coord() ; 
    
}

// ---------------------
// Constructor from file
// ---------------------
Map_et::Map_et(const Mg3d& mgi, FILE* fich) 
	      : Map_radial(mgi, fich),
	        aasx( mgi.get_nr(0) ), 
	        aasx2( mgi.get_nr(0) ), 
	        zaasx( mgi.get_nr(mgi.get_nzone()-1) ), 
	        zaasx2( mgi.get_nr(mgi.get_nzone()-1) ), 
	        bbsx( mgi.get_nr(0) ), 
	        bbsx2( mgi.get_nr(0) ), 
		ff(*(mgi.get_angu()), fich) , 
		gg(*(mgi.get_angu()), fich)   
{
    // The Coord rsxdxdr and rsx2drdx are constructed by the default Coord 
    // constructor
    
    // alpha and beta
    // --------------
    int nz = mg->get_nzone() ;
    alpha = new double[nz] ;
    beta = new double[nz] ;
    fread_be(alpha, sizeof(double), nz, fich) ;	
    fread_be(beta, sizeof(double), nz, fich) ;	

    // Assignement of the building functions of the Coord's
    // ----------------------------------------------------
    set_coord() ; 
    
    // Radial polynomials A(x) and B(x)
    // --------------------------------
    
    fait_poly() ; 

}

			//------------//
			// Destructor //
			//------------//

Map_et::~Map_et() {
    
    delete [] alpha ;
    delete [] beta ;
    
    for (int l=0 ; l<mg->get_nzone(); l++) {
	delete aa[l] ;
	delete daa[l] ;
	delete ddaa[l] ;
	delete bb[l] ;
	delete dbb[l] ;
	delete ddbb[l] ;
    }
    delete [] aa ; 
    delete [] daa ; 
    delete [] ddaa ; 
    delete [] bb ; 
    delete [] dbb ; 
    delete [] ddbb ; 
    
}


			    //------------//
			    // Assignment //
			    //------------//

void Map_et::operator=(const Map_et& mpi) {
    
    assert(mpi.get_mg() == mg) ; 
    
    set_ori( mpi.get_ori_x(), mpi.get_ori_y(), mpi.get_ori_z() ) ; 
    
    set_rot_phi( mpi.get_rot_phi() ) ; 
    
    // The members bvect_spher and bvect_cart are treated by the functions
    //  set_ori and set_rot_phi.
    
    for (int l=0; l<mg->get_nzone(); l++){
	alpha[l] = mpi.get_alpha()[l] ; 
	beta[l] = mpi.get_beta()[l] ; 
    }
    
    ff = mpi.ff ; 
    gg = mpi.gg ; 
    
    reset_coord() ;	// update of all the Coords
        
}




void Map_et::operator=(const Map_af& mpi) {
    
    assert(mpi.get_mg() == mg) ; 
    
    set_ori( mpi.get_ori_x(), mpi.get_ori_y(), mpi.get_ori_z() ) ; 
    
    set_rot_phi( mpi.get_rot_phi() ) ; 
    
    // The members bvect_spher and bvect_cart are treated by the functions
    //  set_ori and set_rot_phi.
    
    for (int l=0; l<mg->get_nzone(); l++){
	alpha[l] = mpi.get_alpha()[l] ; 
	beta[l] = mpi.get_beta()[l] ; 
    }
    
    ff = 0 ; 
    gg = 0 ; 
    
    reset_coord() ;	// update of all the Coords
        
}

void Map_et::set_ff(const Valeur& ffi) {
    
    ff = ffi ; 
    
    reset_coord() ;	// update of all the Coords
    
}

void Map_et::set_gg(const Valeur& ggi) {
    
    gg = ggi ; 
    
    reset_coord() ;	// update of all the Coords
    
}



	    //-------------------------------------------------//
	    //  Assignment of the Coord building functions    //
	    //-------------------------------------------------//
	    
void Map_et::set_coord(){

    // ... Coord's introduced by the base class Map : 
    r.set(this, map_et_fait_r) ;
    tet.set(this, map_et_fait_tet) ;
    phi.set(this, map_et_fait_phi) ;
    sint.set(this, map_et_fait_sint) ;
    cost.set(this, map_et_fait_cost) ;
    sinp.set(this, map_et_fait_sinp) ;
    cosp.set(this, map_et_fait_cosp) ;

    x.set(this, map_et_fait_x) ;
    y.set(this, map_et_fait_y) ;
    z.set(this, map_et_fait_z) ;

    xa.set(this, map_et_fait_xa) ;
    ya.set(this, map_et_fait_ya) ;
    za.set(this, map_et_fait_za) ;
    
    // ... Coord's introduced by the base class Map_radial : 
    xsr.set(this, map_et_fait_xsr) ;
    dxdr.set(this, map_et_fait_dxdr) ;
    drdt.set(this, map_et_fait_drdt) ;
    stdrdp.set(this, map_et_fait_stdrdp) ;
    srdrdt.set(this, map_et_fait_srdrdt) ;
    srstdrdp.set(this, map_et_fait_srstdrdp) ;
    sr2drdt.set(this, map_et_fait_sr2drdt) ;
    sr2stdrdp.set(this, map_et_fait_sr2stdrdp) ;
    d2rdx2.set(this, map_et_fait_d2rdx2) ;
    lapr_tp.set(this, map_et_fait_lapr_tp) ;
    d2rdtdx.set(this, map_et_fait_d2rdtdx) ;
    sstd2rdpdx.set(this, map_et_fait_sstd2rdpdx) ;
    sr2d2rdt2.set(this, map_et_fait_sr2d2rdt2) ;
    
    //... Coord's which belong to the class Map_et only : 
    rsxdxdr.set(this, map_et_fait_rsxdxdr) ; 
    rsx2drdx.set(this, map_et_fait_rsx2drdx) ; 

}

		    //--------------------------//
		    //  Reset of the Coord's	//
		    //--------------------------//

void Map_et::reset_coord() {

    // Coord's of all the class derived from Map_radial:
    
    Map_radial::reset_coord() ;
    
    // Coord's specific to Map_et
    
    rsxdxdr.del_t() ; 
    rsx2drdx.del_t() ; 
     
}
    
	    //------------------------------------------------------//
	    // Construction of the radial polynomials A(x) and B(x) //
	    //------------------------------------------------------//

void Map_et::fait_poly() {
    
    int nzone = mg->get_nzone() ;

    aa = new Tbl*[nzone] ;
    daa = new Tbl*[nzone] ;
    ddaa = new Tbl*[nzone] ;
    bb = new Tbl*[nzone] ;
    dbb = new Tbl*[nzone] ;
    ddbb = new Tbl*[nzone] ;

    for (int l=0; l<nzone; l++) {
	int nr = mg->get_nr(l) ;
	aa[l] = new Tbl(nr) ;
	daa[l] = new Tbl(nr) ;
	ddaa[l] = new Tbl(nr) ;
	bb[l] = new Tbl(nr) ;
	dbb[l] = new Tbl(nr) ;
	ddbb[l] = new Tbl(nr) ;
    }

    // Values in the nucleus
    // ---------------------
    assert( mg->get_type_r(0) == RARE ) ;
    
    aa[0]->set_etat_qcq() ;	    // Memory allocation for the Tbl 
    daa[0]->set_etat_qcq() ; 
    ddaa[0]->set_etat_qcq() ; 
    aasx.set_etat_qcq() ; 
    aasx2.set_etat_qcq() ; 

    bb[0]->set_etat_qcq() ;	    
    dbb[0]->set_etat_qcq() ; 
    ddbb[0]->set_etat_qcq() ; 
    bbsx.set_etat_qcq() ; 
    bbsx2.set_etat_qcq() ; 
    
    for (int i=0; i<mg->get_nr(0); i++) {

	double x1 =  (mg->get_grille3d(0))->x[i]  ;
	double x2 = x1 * x1 ;
	double x3 = x1 * x2 ; 

				//##...... A(x) = 2 x^2 - x^4 :	
				//	(aa[0])->t[i] = x2 * (2. - x2) ;
				//	(daa[0])->t[i] = 4. * x * (1. + x) * (1. - x)  ;
				//	(ddaa[0])->t[i] = 4. *(1. - 3.* x2)  ;
				//	aasx->t[i] = x * (2. - x2) ;
				//	aasx2->t[i] = 2. - x2 ;

	//...... A(x) = 3 x^4 - 2 x^6 :	

	aa[0]->set(i) = x2 * x2 * (3. - 2.*x2) ;
	daa[0]->set(i) = 12. * x3 * (1. + x1) * (1. - x1)  ;
	ddaa[0]->set(i) = 12. *x2 *(3. - 5.* x2)  ;
	aasx.set(i) = x3 * (3. - 2.*x2) ;
	aasx2.set(i) = x2 * (3. - 2.*x2) ;

	//... B(x) = 5/2 x^3 - 3/2 x^5   :  

	bb[0]->set(i) = 0.5 * x3 * (5. - 3.* x2) ; 
	dbb[0]->set(i) = 7.5 * x2 * (1. + x1) * (1. - x1) ;
	ddbb[0]->set(i) =  15. * x1 * ( 1. - 2.*x2 )   ;
	bbsx.set(i) = 0.5 * x2 * (5. - 3.* x2) ;
	bbsx2.set(i) = 0.5 * x1 * (5. - 3.* x2) ;
    }
    
    // Values in the shells and the outermost domain
    // ---------------------------------------------

    for (int l=1; l<nzone; l++) {

	assert( (mg->get_type_r(l) == FIN)|| (mg->get_type_r(l) == UNSURR) ) ;

	aa[l]->set_etat_qcq() ;	    // Memory allocation for the Tbl 
	daa[l]->set_etat_qcq() ; 
	ddaa[l]->set_etat_qcq() ; 

	bb[l]->set_etat_qcq() ;	    
	dbb[l]->set_etat_qcq() ; 
	ddbb[l]->set_etat_qcq() ; 
	
	for (int i=0; i<mg->get_nr(l); i++) {

	    double x1 = (mg->get_grille3d(l))->x[i] ;
	    double xm1 = x1 - 1. ; 
	    double xp1 = x1 + 1. ; 

	    //... A(x) = 1/4 (x-1)^2 (x+2) = 1/4(x^3 -3x +2)  : 

	    aa[l]->set(i) = 0.25* xm1 * xm1 * (x1 + 2.) ;
	    daa[l]->set(i) = 0.75* xm1 * xp1  ;
	    ddaa[l]->set(i) = 1.5* x1 ;

	    //... B(x) = 1/4 (x+1)^2 (-x+2) = 1/4(-x^3 +3x +2)   : 

	    bb[l]->set(i) = 0.25* xp1 * xp1 * (2. - x1) ;
	    dbb[l]->set(i) = - 0.75* xm1 * xp1  ;
	    ddbb[l]->set(i) = - 1.5* x1  ;
	}

    }	    // End of the loop onto the domains

    // Special case of a compactified outermost domain
    // -----------------------------------------------
    
    int nzm1 = nzone - 1 ; 
    if ( mg->get_type_r(nzm1) == UNSURR ) {

	zaasx.set_etat_qcq() ;	    // Memory allocation for the Tbl 
	zaasx2.set_etat_qcq() ;
	
	for (int i=0; i<mg->get_nr(nzm1); i++) {
	    
	    double x1 = (mg->get_grille3d(nzm1))->x[i] ; 
	    zaasx.set(i) = 0.25 * (x1 - 1.) * (2. + x1) ;	    // A(x)/(x-1)
	    zaasx2.set(i) = 0.25 * (2. + x1) ;		    // A(x)/(x-1)^2
	    
	}
	
	bb[nzm1]->set_etat_zero() ; 
	dbb[nzm1]->set_etat_zero() ; 
	ddbb[nzm1]->set_etat_zero() ; 
    
    }  
  
}



			//----------------//
			// Save in a file //
			//----------------//

void Map_et::sauve(FILE* fich) const {

    Map_radial::sauve(fich) ;	// Write of the elements common to all the
				// classes derived from Map_radial
    
    ff.sauve(fich) ;		// Write of F(theta',phi')
    gg.sauve(fich) ;		// Write of G(theta',phi')
    
    // Write of alpha and beta :
    int nz = mg->get_nzone() ;
    fwrite_be(alpha, sizeof(double), nz, fich) ;	
    fwrite_be(beta, sizeof(double), nz, fich) ;	
    
}

		//---------------------------//
		//        Printing	     //
		//---------------------------//

ostream & Map_et::operator>>(ostream & ost) const {

  using namespace Unites ;

    ost << 
    "Radial mapping of form r = xi + A(xi)F(t,p) + B(xi)G(t,p) (class Map_et)" 
    << endl ; 
    int nz = mg->get_nzone() ;
    for (int l=0; l<nz; l++) {
	ost << "     Domain #" << l << " : alpha_l = " << alpha[l] 
	  << " ,  beta_l = " << beta[l] << endl ;  
    }
    ost << endl << "Function F(theta', phi') : " << endl ; 
    ost << "-------------------------  " << endl ; 
    ff.affiche_seuil(ost) ; 
    ost << endl <<"Function G(theta', phi') : " << endl ; 
    ost << "-------------------------  " << endl ; 
    gg.affiche_seuil(ost) ; 

    int type_t = mg->get_type_t() ; 
    int type_p = mg->get_type_p() ; 

    ost << endl 
	<< "Values of r at the outer boundary of each domain [km] :" << endl ; 
    ost << "------------------------------------------------------" << endl ; 
    ost << "   1/ for theta = Pi/2 and phi = 0 : " << endl ; 
    ost << "       val_r :   " ;
    for (int l=0; l<nz; l++) {
	ost << " " << val_r(l, 1., M_PI/2, 0) / km ; 
    }
    ost << endl ; 

    if ( type_t == SYM ) {
	assert( (type_p == SYM) || (type_p == NONSYM) ) ; 
	ost << "       Coord r : " ; 
	for (int l=0; l<nz; l++) {
	    int nrm1 = mg->get_nr(l) - 1 ; 
	    int ntm1 = mg->get_nt(l) - 1 ; 
	    ost << " " << (+r)(l, 0, ntm1, nrm1) / km ; 
	}
	ost << endl ; 
    }

    ost << "   2/ for theta = Pi/2 and phi = Pi/2 : " << endl ; 
    ost << "       val_r :   " ;
    for (int l=0; l<nz; l++) {
	ost << " " << val_r(l, 1., M_PI/2, M_PI/2) / km ; 
    }
    ost << endl ; 

    if ( type_t == SYM ) {
	ost << "       Coord r : " ; 
	for (int l=0; l<nz; l++) {
	    int nrm1 = mg->get_nr(l) - 1 ; 
	    int ntm1 = mg->get_nt(l) - 1 ; 
	    int np = mg->get_np(l) ;
	    if ( (type_p == NONSYM) && (np % 4 == 0) ) {
		ost << " " << (+r)(l, np/4, ntm1, nrm1) / km ; 
	    }
	    if ( type_p == SYM ) {
		ost << " " << (+r)(l, np/2, ntm1, nrm1) / km ; 
	    }
	}
	ost << endl ; 
    }

    ost << "   3/ for theta = Pi/2 and phi = Pi : " << endl ; 
    ost << "       val_r :   " ;
    for (int l=0; l<nz; l++) {
	ost << " " << val_r(l, 1., M_PI/2, M_PI) / km ; 
    }
    ost << endl ; 

    if ( (type_t == SYM) && (type_p == NONSYM) ) {
	ost << "       Coord r : " ; 
	for (int l=0; l<nz; l++) {
	    int nrm1 = mg->get_nr(l) - 1 ; 
	    int ntm1 = mg->get_nt(l) - 1 ; 
	    int np = mg->get_np(l) ;
	    ost << " " << (+r)(l, np/2, ntm1, nrm1) / km ; 
	}
	ost << endl ; 
    }

    ost << "   4/ for theta = 0 : " << endl ; 
    ost << "       val_r :   " ;
    for (int l=0; l<nz; l++) {
	ost << " " << val_r(l, 1., 0., 0.) / km ; 
    }
    ost << endl ; 

    ost << "       Coord r : " ; 
    for (int l=0; l<nz; l++) {
	int nrm1 = mg->get_nr(l) - 1 ; 
	ost << " " << (+r)(l, 0, 0, nrm1) / km ; 
    }
    ost << endl ; 

    return ost ;
    
}    

			//------------------//
			//  Homothetie	    //
			//------------------//


void Map_et::homothetie(double fact) {

    int nz = mg->get_nzone() ; 
    
    for (int l=0; l<nz; l++) {
	if (mg->get_type_r(l) == UNSURR) {
	    alpha[l] /= fact ;
	    beta[l] /= fact ;
	}
	else {
	    alpha[l] *= fact ;
	    beta[l] *= fact ;
	}
    }
    
    reset_coord() ; 
    
}

		    //----------------------------//
		    //	Rescaling of one domain	  //
		    //----------------------------//

void Map_et::resize(int l, double lambda) {

    // Protections
    // -----------
    if (mg->get_type_r(l) != FIN) {
	cout << "Map_et::resize can be applied only to a shell !" << endl ; 
	abort() ;
    }

    // New values of alpha, beta, F and G in domain l :
    // ----------------------------------------------
    double n_alpha = 0.5 * ( (lambda + 1.) * alpha[l] +  
			     (lambda - 1.) * beta[l] ) ; 

    double n_beta = 0.5 * ( (lambda - 1.) * alpha[l] +  
			    (lambda + 1.) * beta[l] ) ; 
			    
    ff.set(l) = alpha[l] / n_alpha  * ff(l) ; 
    gg.set(l) = lambda * alpha[l] / n_alpha  * gg(l) ; 

    alpha[l] = n_alpha ; 
    beta[l] = n_beta ; 
    
    // New values of alpha, beta, F and G in  in domain l+1 :
    // ----------------------------------------------------
    assert(l<mg->get_nzone()-1) ; 
    int lp1 = l + 1 ; 
    
    if (mg->get_type_r(lp1) == UNSURR) {	    // compactified case
	
	assert(ff(lp1).get_etat() == ETATZERO ) ; 
	assert(gg(lp1).get_etat() == ETATZERO ) ; 
	
	alpha[lp1] = - 0.5 / ( alpha[l] + beta[l] ) ; 
	beta[lp1] = - alpha[lp1] ; 
	
    }
    else{	// non-compactified case
	
	assert( mg->get_type_r(lp1) == FIN ) ;
	n_alpha = 0.5 * ( alpha[lp1] - alpha[l] + beta[lp1] - beta[l] ) ; 
	n_beta =  0.5 * ( alpha[lp1] + alpha[l] + beta[lp1] + beta[l] ) ; 
	
	ff.set(lp1) = alpha[l] / n_alpha * gg(l) ; 
	gg.set(lp1) = alpha[lp1] / n_alpha * gg(lp1) ; 
	
	alpha[lp1] = n_alpha ; 
	beta[lp1] = n_beta ; 
    }
    
    // The coords are no longer up to date :
    reset_coord() ; 
} 


// Comparison operator :
bool Map_et::operator==(const Map& mpi) const {
  
  // Precision of the comparison
  double precis = 1e-10 ;
  bool resu = true ;

  // Dynamic cast pour etre sur meme Map...
  const Map_et* mp0 = dynamic_cast<const Map_et*>(&mpi) ;
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

    int nz = mg->get_nzone() ;
    for (int i=0 ; i<nz ; i++) {
      if (fabs(alpha[i]-mp0->alpha[i])/fabs(alpha[i]) > precis) 
	resu = false ;
      if ((i!=0) && (i!=nz-1))
	if (fabs(beta[i]-mp0->beta[i])/fabs(beta[i]) > precis) 
	resu = false ;
    }

    if (max(diffrelmax(ff, mp0->ff)) > precis)
      resu = false ;
    if (max(diffrelmax(gg, mp0->gg)) > precis)
      resu = false ;
  }

  return resu ;
}
		//--------------------------------------//
		// Extraction of the mapping parameters //
		//--------------------------------------//

const double* Map_et::get_alpha() const {
    return alpha ; 
}

const double* Map_et::get_beta() const {
    return beta ; 
}

const Valeur& Map_et::get_ff() const {
    return ff ; 
}

const Valeur& Map_et::get_gg() const {
    return gg ; 
}
