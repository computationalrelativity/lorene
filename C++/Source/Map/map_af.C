/*
 *  Methods of class Map_af
 *
 *   (see file map.h for documentation)
 *
 */

/*
 *   Copyright (c) 1999-2003 Eric Gourgoulhon
 *   Copyright (c) 1999-2001 Philippe Grandclement
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


char map_af_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.8  2005/08/29 15:10:18  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 * Revision 1.7  2004/12/02 09:33:06  p_grandclement
 * *** empty log message ***
 *
 * Revision 1.6  2004/03/25 10:29:23  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.5  2004/01/29 08:50:03  p_grandclement
 * Modification of Map::operator==(const Map&) and addition of the surface
 * integrales using Scalar.
 *
 * Revision 1.4  2003/10/15 10:33:11  e_gourgoulhon
 * Added new Coord's : drdt, srdrdp.
 *
 * Revision 1.3  2002/10/16 14:36:41  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
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
 * Revision 2.23  2001/02/28  11:04:05  eric
 * 1ere version testee de resize.
 *
 * Revision 2.22  2001/02/26  17:29:22  eric
 * Ajout de la fonction resize.
 *
 * Revision 2.21  2001/01/10  11:03:13  phil
 * ajout de homothetie interne
 *
 * Revision 2.20  2000/01/24  17:09:13  eric
 * suppression de la fonction convert.
 * suppression du constructeur par convertion d'un Map_et.
 * ajout du constructeur par conversion d'un Map.
 *
 * Revision 2.19  2000/01/24  16:42:48  eric
 * Dans operator=(const Map_af& ), appel de set_rot_phi.
 * Ajout de la fonction convert(const Map& ).
 *
 * Revision 2.18  1999/12/21  16:27:26  eric
 * Ajout du constructeur par conversion Map_af::Map_af(const Map_et&).
 * Ajout des fonctions Map_af::set_alpha et Map_af::set_beta.
 *
 * Revision 2.17  1999/12/17  09:14:43  eric
 * Amelioration de l'affichage.
 *
 * Revision 2.16  1999/12/06  13:12:23  eric
 * Introduction de la fonction homothetie.
 *
 * Revision 2.15  1999/11/25  16:28:53  eric
 * Le calcul des derivees partielles est transfere dans le fichier
 *   map_af_deriv.C.
 *
 * Revision 2.14  1999/11/24  14:32:54  eric
 * Les prototypes des fonctions de constructions des coords sont desormais
 *   dans map.h.
 * Introduction des fonctions get_alpha() et get_beta().
 * /
 *
 * Revision 2.13  1999/11/22  10:36:14  eric
 * Introduction de la fonction set_coord().
 * Fonction del_coord() rebaptisee reset_coord().
 *
 * Revision 2.12  1999/10/27  08:46:29  eric
 * Introduction de Cmp::va a la place de *(Cmp::c).
 *
 * Revision 2.11  1999/10/15  09:23:10  eric
 * *** empty log message ***
 *
 * Revision 2.10  1999/10/15  09:16:23  eric
 * Changement prototypes: const.
 *
 * Revision 2.9  1999/10/14  14:27:05  eric
 * Depoussierage.
 *
 * Revision 2.8  1999/04/12  12:54:05  phil
 * *** empty log message ***
 *
 * Revision 2.7  1999/04/12  12:09:03  phil
 * Mise a jour des bases dans dsdr()
 *
 * Revision 2.6  1999/03/04  13:11:48  eric
 * Ajout des Coord representant les derivees du changement de variable.
 *
 * Revision 2.5  1999/03/03  11:19:08  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// headers C
#include <math.h>

// headers Lorene
#include "coord.h"
#include "map.h"
#include "cmp.h"
#include "utilitaires.h"
#include "unites.h"

// Local prototype:
void c_est_pas_fait(char * ) ;
			
			//---------------//
			// Constructeurs //
			//---------------//

// Constructor from a grid
// -----------------------
Map_af::Map_af(const Mg3d& mgrille, const double* bornes) : Map_radial(mgrille)
{
    // Les coordonnees et les derivees du changement de variable
    set_coord() ; 
    
    // Les bornes
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
		double umax = 1./bornes[l] ;
		double umin = 1./bornes[l+1] ;
		alpha[l] = (umin - umax) * .5 ;  // u est une fonction decroissante
		beta[l] = (umin + umax) * .5 ;   //  de l'indice i en r
		break ;
	    }
	    
	    default:	{
		cout << "Map_af::Map_af: unkown type_r ! " << endl ;
		abort () ;
		break ;
	    }
	    
	}
    }	    // Fin de la boucle sur zone
}

// Constructor from a grid
// -----------------------
Map_af::Map_af(const Mg3d& mgrille, const Tbl& bornes) : Map_radial(mgrille)
{
    // Les coordonnees et les derivees du changement de variable
    set_coord() ; 
    
    // Les bornes
    int nzone = mg->get_nzone() ;
    
    alpha = new double[nzone] ;
    beta = new double[nzone] ;

    for (int l=0 ; l<nzone ; l++) {
	switch (mg->get_type_r(l)) {
	    case RARE:	{
		alpha[l] = bornes(l+1) - bornes(l) ;
		beta[l] = bornes(l) ;
		break ; 
	    }
	    
	    case FIN:	{
		alpha[l] = (bornes(l+1) - bornes(l)) * .5 ;
		beta[l] = (bornes(l+1) + bornes(l)) * .5 ;
		break ;
	    }
	    
	    case UNSURR: {
	    	assert (l==nzone-1) ;
		double umax = 1./bornes(l) ;
		double umin = 0 ;
		alpha[l] = (umin - umax) * .5 ;  // u est une fonction decroissante
		beta[l] = (umin + umax) * .5 ;   //  de l'indice i en r
		break ;
	    }
	    
	    default:	{
		cout << "Map_af::Map_af: unkown type_r ! " << endl ;
		abort () ;
		break ;
	    }
	    
	}
    }	    // Fin de la boucle sur zone
}

// Copy constructor 
// ----------------
Map_af::Map_af(const Map_af& mp) : Map_radial(mp)
{
    // Les coordonnees et les derivees du changement de variable
    set_coord() ; 
        
    // Les bornes
    int nzone = mg->get_nzone() ;
    
    alpha = new double[nzone] ;
    beta = new double[nzone] ;

    for (int l=0; l<nzone; l++){
	alpha[l] = mp.alpha[l] ; 
	beta[l] = mp.beta[l] ; 
    }
}
	
// Constructor from file
// ---------------------
Map_af::Map_af(const Mg3d& mgi, FILE* fd) : Map_radial(mgi, fd)
{
    int nz = mg->get_nzone() ;
    alpha = new double[nz] ;
    beta = new double[nz] ;
    fread_be(alpha, sizeof(double), nz, fd) ;	
    fread_be(beta, sizeof(double), nz, fd) ;	

    // Les coordonnees et les derivees du changement de variable
    set_coord() ; 
}



// Constructor from a Map
// -----------------------
Map_af::Map_af(const Map& mpi) : Map_radial(*(mpi.get_mg()))
{
    // Les coordonnees et les derivees du changement de variable
    set_coord() ; 
        
    // Les bornes
    int nz = mg->get_nzone() ; 
    
    alpha = new double[nz] ;
    beta = new double[nz] ;

    const Map_af* mp0 = dynamic_cast<const Map_af*>(&mpi) ;
    const Map_et* mp1 = dynamic_cast<const Map_et*>(&mpi) ;
     
    if( (mp0 == 0x0) && (mp1 == 0x0) ) {
	cout << "Map_af::Map_af(const Map& ) : unkown mapping type !" 
	     << endl ; 
	abort() ; 
    }

    if (mp0 != 0x0) {
	assert( mp1 == 0x0 ) ; 
	for (int l=0; l<nz; l++){
	    alpha[l] = mp0->get_alpha()[l] ; 
	    beta[l] = mp0->get_beta()[l] ; 
	}
    }

	
    if (mp1 != 0x0) {
	assert( mp0 == 0x0 ) ; 
	for (int l=0; l<nz; l++){
	    alpha[l] = mp1->get_alpha()[l] ; 
	    beta[l] = mp1->get_beta()[l] ; 
	}
    }

    
    set_ori( mpi.get_ori_x(), mpi.get_ori_y(), mpi.get_ori_z() ) ; 
    
    set_rot_phi( mpi.get_rot_phi() ) ; 
    
}




			//--------------//
			// Destructeurs //
			//--------------//

Map_af::~Map_af() {
    delete [] alpha ;
    delete [] beta ;
}

			//-------------//
			// Assignement //
			//-------------//
			
// From another Map_af
// -------------------

void Map_af::operator=(const Map_af & mpi) {
    
    assert(mpi.mg == mg) ;
    
    set_ori( mpi.ori_x, mpi.ori_y, mpi.ori_z ) ; 
    
    set_rot_phi( mpi.rot_phi ) ; 

    for (int l = 0; l<mg->get_nzone(); l++) {
	alpha[l] = mpi.alpha[l] ; 
	beta[l] = mpi.beta[l] ; 
    }

    reset_coord() ;      
}
    



	    //-------------------------------------------------//
	    //  Assignement of the Coord building functions    //
	    //-------------------------------------------------//
	    
void Map_af::set_coord(){

    // ... Coord's introduced by the base class Map : 
    r.set(this, map_af_fait_r) ;
    tet.set(this, map_af_fait_tet) ;
    phi.set(this, map_af_fait_phi) ;
    sint.set(this, map_af_fait_sint) ;
    cost.set(this, map_af_fait_cost) ;
    sinp.set(this, map_af_fait_sinp) ;
    cosp.set(this, map_af_fait_cosp) ;

    x.set(this, map_af_fait_x) ;
    y.set(this, map_af_fait_y) ;
    z.set(this, map_af_fait_z) ;

    xa.set(this, map_af_fait_xa) ;
    ya.set(this, map_af_fait_ya) ;
    za.set(this, map_af_fait_za) ;
    
    // ... Coord's introduced by the base class Map_radial : 
    xsr.set(this, map_af_fait_xsr) ;
    dxdr.set(this, map_af_fait_dxdr) ;
    drdt.set(this, map_af_fait_drdt) ;
    stdrdp.set(this, map_af_fait_stdrdp) ;
    srdrdt.set(this, map_af_fait_srdrdt) ;
    srstdrdp.set(this, map_af_fait_srstdrdp) ;
    sr2drdt.set(this, map_af_fait_sr2drdt) ;
    sr2stdrdp.set(this, map_af_fait_sr2stdrdp) ;
    d2rdx2.set(this, map_af_fait_d2rdx2) ;
    lapr_tp.set(this, map_af_fait_lapr_tp) ;
    d2rdtdx.set(this, map_af_fait_d2rdtdx) ;
    sstd2rdpdx.set(this, map_af_fait_sstd2rdpdx) ;
    sr2d2rdt2.set(this, map_af_fait_sr2d2rdt2) ;
    
}
// Comparison operator :
bool Map_af::operator==(const Map& mpi) const {
  
  // Precision of the comparison
  double precis = 1e-10 ;
  bool resu = true ;

  // Dynamic cast pour etre sur meme Map...
  const Map_af* mp0 = dynamic_cast<const Map_af*>(&mpi) ;
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
  }

  return resu ;
}

		//--------------------------------------//
		// Extraction of the mapping parameters //
		//--------------------------------------//

const double* Map_af::get_alpha() const {
    return alpha ; 
}

const double* Map_af::get_beta() const {
    return beta ; 
}

			//------------//
			// Sauvegarde //
			//------------//

void Map_af::sauve(FILE* fd) const {

    Map_radial::sauve(fd) ; 

    int nz = mg->get_nzone() ;
    fwrite_be(alpha, sizeof(double), nz, fd) ;	
    fwrite_be(beta, sizeof(double), nz, fd) ;	
    
}
    
			//------------//
			// Impression //
			//------------//

ostream & Map_af::operator>>(ostream & ost) const {

  using namespace Unites ;

    ost << "Affine mapping (class Map_af)" << endl ; 
    int nz = mg->get_nzone() ;
    for (int l=0; l<nz; l++) {
	ost << "     Domain #" << l << " : alpha_l = " << alpha[l] 
	  << " ,  beta_l = " << beta[l] << endl ;  
    }

    ost << endl << "     Values of r at the outer boundary of each domain [km] :" 
	<< endl ; 
    ost << "            val_r :   " ;
    for (int l=0; l<nz; l++) {
	ost << " " << val_r(l, 1., 0., 0.) / km ; 
    }
    ost << endl ; 

    ost << "            Coord r : " ; 
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


void Map_af::homothetie(double fact) {

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

void Map_af::resize(int l, double lambda) {

    // Protections
    // -----------
    if (mg->get_type_r(l) != FIN) {
	cout << "Map_af::resize can be applied only to a shell !" << endl ; 
	abort() ;
    }

    // New values of alpha and beta in domain l :
    // ----------------------------------------
    double n_alpha = 0.5 * ( (lambda + 1.) * alpha[l] +  
			     (lambda - 1.) * beta[l] ) ; 

    double n_beta = 0.5 * ( (lambda - 1.) * alpha[l] +  
			    (lambda + 1.) * beta[l] ) ; 
			    
    alpha[l] = n_alpha ; 
    beta[l] = n_beta ; 
    
    // New values of alpha and beta in domain l+1 :
    // ------------------------------------------
    assert(l<mg->get_nzone()-1) ; 
    int lp1 = l + 1 ; 
    
    if (mg->get_type_r(lp1) == UNSURR) {	    // compactified case
	
	alpha[lp1] = - 0.5 / ( alpha[l] + beta[l] ) ; 
	beta[lp1] = - alpha[lp1] ; 
	
    }
    else{	// non-compactified case
	
	assert( mg->get_type_r(lp1) == FIN ) ;
	n_alpha = 0.5 * ( alpha[lp1] - alpha[l] + beta[lp1] - beta[l] ) ; 
	n_beta =  0.5 * ( alpha[lp1] + alpha[l] + beta[lp1] + beta[l] ) ; 
	alpha[lp1] = n_alpha ; 
	beta[lp1] = n_beta ; 
    }
    
    // The coords are no longer up to date :
    reset_coord() ; 
    
} 

		    



			//---------------------------//
			//  Homothetie	partielle   //
			//-------------------------//


void Map_af::homothetie_interne(double fact) {
    
    // Dans le noyau
    alpha[0] *= fact ;

    // Dans la premiere coquille :
    double asauve  = alpha[1] ;
    alpha[1] = (1-fact)/2.*beta[1] + (1+fact)/2. * alpha[1] ;
    beta[1] = (1+fact)/2.*beta[1]+ (1-fact)/2. * asauve ;
    
    reset_coord() ; 
}
			//------------------------------------------//
			//  Modification of the mapping parameters  //
			//------------------------------------------//

void Map_af::set_alpha(double alpha0, int l) {
    
    assert(l>=0) ; 
    assert(l<mg->get_nzone()) ; 
    
    alpha[l] = alpha0 ; 
    
    reset_coord() ; 
    
}

void Map_af::set_beta(double beta0, int l) {
    
    assert(l>=0) ; 
    assert(l<mg->get_nzone()) ; 
    
    beta[l] = beta0 ; 
    
    reset_coord() ; 
    
}


// To be done
//-----------

void Map_af::adapt(const Cmp&, const Param&) {
    char* f = __FILE__ ;
    c_est_pas_fait(f) ;
}



