/*
 *  Methods of class Map_eps
 *
 *   (see file map.h for documentation)
 *
 */


/*
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


// headers C++
#include <stdexcept>

// headers C
#include <cmath>

// headers Lorene
#include "valeur.h"
#include "map.h"
#include "utilitaires.h"
#include "proto.h"
#include "unites.h"



			//---------------//
			// Constructeurs //
			//---------------//

// Constructor from a grid
// -----------------------
namespace Lorene {
Map_eps::Map_eps(const Mg3d& mgrille, const double* bornes) : Map_radial(mgrille),alpha(mgrille.get_angu()),beta(mgrille.get_angu())
{
    c_est_pas_fait("Map_eps::Map_eps(const Mg3d&,  const double*)") ;
}

// Constructor from a grid
// -----------------------
Map_eps::Map_eps(const Mg3d& mgrille, double a, double b, double c) : Map_radial(mgrille),aa(a),bb(b),cc(c),alpha(mgrille.get_angu()),beta(mgrille.get_angu())
{
    assert(aa > 0.) ;
    assert(bb > 0.) ;
    assert(cc > 0.) ;
    // Les coordonnees et les derivees du changement de variable
    set_coord() ; 

    // Allocate memory
    alpha.annule_hard() ;
    beta.annule_hard() ;
    
    // Les bornes
    int nzone = mg->get_nzone() ;
    assert(nzone==1) ; // Only one zone thank you
    const Grille3d* g = mg->get_grille3d(0) ;

    for (int l=0 ; l<nzone ; l++) {
        for (int k=0 ; k<mg->get_np(l); k++){
            for (int j=0; j<mg->get_nt(l); j++){
	switch (mg->get_type_r(l)) {
	    case RARE:	{
        double theta = (g->tet)[j] ;
        double pphi  = (g->phi)[k] ;
		alpha.set(l,k,j,0) = 1/sqrt(pow(cos(pphi)*sin(theta)/a, 2.) + pow(sin(pphi)*sin(theta)/b, 2.) + pow(cos(theta)/c, 2.)) ;
        beta.set(l) = 0. ;
		break ; 
	    }
	    
	    case FIN:	{
        cout << "Warning ! No shell allowed !" << endl;
        abort() ;
		break ;
	    }
	    
	    case UNSURR: {
	    cout << "Warning ! No compactified domain allowed !" << endl;
        abort() ;
        break ;
	    }
	    
	    default:	{
		cout << "Map_eps::Map_eps: unkown type_r ! " << endl ;
		abort () ;
		break ;
	    }
	}
    alpha.std_base_scal() ;
    beta.std_base_scal() ;
    }
    }
    }	    // Fin de la boucle sur zone
}

// Copy constructor 
// ----------------
Map_eps::Map_eps(const Map_eps& mp) : Map_radial(mp),alpha(mp.alpha),beta(mp.beta)
{
    // Les coordonnees et les derivees du changement de variable
    set_coord() ; 
}
      

// Constructor from file
// ---------------------
Map_eps::Map_eps(const Mg3d& mgi, FILE* fd) : Map_radial(mgi, fd),alpha(*mgi.get_angu(), fd),beta(*mgi.get_angu(), fd)
{
    // Les coordonnees et les derivees du changement de variable
    set_coord() ;
}

			//--------------//
			// Destructeurs //
			//--------------//

Map_eps::~Map_eps() {

}

			//-------------//
			// Assignement //
			//-------------//
			
// From another Map_eps
// -------------------

void Map_eps::operator=(const Map_eps & mpi) {
    
    assert(mpi.mg == mg) ;
    
    set_ori( mpi.ori_x, mpi.ori_y, mpi.ori_z ) ; 
    
    set_rot_phi( mpi.rot_phi ) ; 

	alpha = mpi.alpha ; 
    beta = mpi.beta ;

    reset_coord() ;      
}

//Assignment to a Map_af.

void Map_eps::operator=(const Map_af & mpi) {
    
    c_est_pas_fait("Map_eps::operator=(const Map_af)") ;    
}


	    //-------------------------------------------------//
	    //  Assignement of the Coord building functions    //
	    //-------------------------------------------------//
	    
void Map_eps::set_coord(){

    // ... Coord's introduced by the base class Map : 
    r.set(this, map_eps_fait_r) ;
    tet.set(this, map_eps_fait_tet) ;
    phi.set(this, map_eps_fait_phi) ;
    sint.set(this, map_eps_fait_sint) ;
    cost.set(this, map_eps_fait_cost) ;
    sinp.set(this, map_eps_fait_sinp) ;
    cosp.set(this, map_eps_fait_cosp) ;

    x.set(this, map_eps_fait_x) ;
    y.set(this, map_eps_fait_y) ;
    z.set(this, map_eps_fait_z) ;

    xa.set(this, map_eps_fait_xa) ;
    ya.set(this, map_eps_fait_ya) ;
    za.set(this, map_eps_fait_za) ;
    
    // ... Coord's introduced by the base class Map_radial : 
    xsr.set(this, map_eps_fait_xsr) ;
    dxdr.set(this, map_eps_fait_dxdr) ;
    drdt.set(this, map_eps_fait_drdt) ;
    stdrdp.set(this, map_eps_fait_stdrdp) ;
    srdrdt.set(this, map_eps_fait_srdrdt) ;
    srstdrdp.set(this, map_eps_fait_srstdrdp) ;
    sr2drdt.set(this, map_eps_fait_sr2drdt) ;
    sr2stdrdp.set(this, map_eps_fait_sr2stdrdp) ;
    d2rdx2.set(this, map_eps_fait_d2rdx2) ;
    lapr_tp.set(this, map_eps_fait_lapr_tp) ;
    d2rdtdx.set(this, map_eps_fait_d2rdtdx) ;
    sstd2rdpdx.set(this, map_eps_fait_sstd2rdpdx) ;
    sr2d2rdt2.set(this, map_eps_fait_sr2d2rdt2) ;
    
}
// Comparison operator :
bool Map_eps::operator==(const Map& mpi) const {
  
  // Precision of the comparison
  double precis = 1e-10 ;
  bool resu = true ;

  // Dynamic cast pour etre sur meme Map...
  const Map_eps* mp0 = dynamic_cast<const Map_eps*>(&mpi) ;
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
      if (max(abs(alpha(i)-mp0->alpha(i))/abs(alpha(i))) > precis) 
	resu = false ;
    }
  }

  return resu ;
}


		//--------------------------------------//
		// Extraction of the mapping parameters //
		//--------------------------------------//

const Valeur& Map_eps::get_alpha() const {
    return alpha ; 
}
const Valeur& Map_eps::get_beta() const {
    return beta ; 
}

			//------------//
			// Sauvegarde //
			//------------//

void Map_eps::sauve(FILE* fd) const {

    Map_radial::sauve(fd) ; 

    alpha.sauve(fd) ; //May not be enough to save everything...	
    beta.sauve(fd) ;

}

			//------------//
			// Impression //
			//------------//

ostream & Map_eps::operator>>(ostream & ost) const {

  using namespace Unites ;

    ost << "Ellipsoidal mapping (class Map_eps)" << endl ; 
    ost << "Parameters of the ellipsoid: (x-axis) a = " << aa << " (y-axis) b = " << bb << " (z-axis) c = " << cc << endl;
    int nz = mg->get_nzone() ;
    for (int l=0; l<nz; l++) {
        for (int k=0 ; k<mg->get_np(l); k++){
            for (int j=0; j<mg->get_nt(l); j++){
	            ost << "     Domain #" << l << " grid point index in theta : " << j << " phi : " << k << '\n' << " : alpha(l,k,j) = " << alpha(l,k,j,0)
                << '\n' << " : beta(l,k,j) = " << beta(l,k,j,0)
                << endl ;
            }
        }
    }

    ost << endl << "     Values of r at the outer boundary of each domain [km] :" 
	<< endl ; 
    ost << "            val_r :   " ;
    for (int l=0; l<nz; l++) {
        for (int k=0 ; k<mg->get_np(l); k++){
            for (int j=0; j<mg->get_nt(l); j++){
	            ost << " " << val_r(l, 1., (+tet)(l,k,j,0), (+phi)(l,k,j,0)) / km ; 
            }
        }
    }
    ost << endl ; 

    ost << "            Coord r : " ; 
    for (int l=0; l<nz; l++) {
        for (int k=0 ; k<mg->get_np(l); k++){
            for (int j=0; j<mg->get_nt(l); j++){
                int nrm1 = mg->get_nr(l) - 1 ; 
                ost << " " << (+r)(l, k, j, nrm1) / km ; 
            }
        }
    }
    ost << endl ; 

    return ost ;
}    

			//------------------------------------------//
			//  Modification of the mapping parameters  //
			//------------------------------------------//

void Map_eps::set_alpha(const Tbl& alpha0, int l) {
    
    assert(l>=0) ; 
    assert(l<mg->get_nzone()) ; 
    
    int Nt = mg->get_nt(l) ;
    int Np = mg->get_np(l) ;
    for (int k=0; k<Np; k++)
        for (int j=0; j<Nt; j++){
           alpha.set(l, k, j, 0) = alpha0(k,j) ; // alpha0 must represent the values of the radius. 
    }
    
    
    reset_coord() ; 
    
}

void Map_eps::set_alpha(const Valeur& alpha0) {
    
    alpha = alpha0 ;
    alpha.std_base_scal() ;
    
    reset_coord() ; 
    
}

void Map_eps::set_beta(const Valeur& beta0) {
    
    beta = beta0 ;
    beta.std_base_scal() ;

    reset_coord() ; 
    
}

void Map_eps::set_beta(const Tbl& beta0, int l) {
    
    assert(l>=0) ; 
    assert(l<mg->get_nzone()) ; 
    
    int Nt = mg->get_nt(l) ;
    int Np = mg->get_np(l) ;
    for (int k=0; k<Np; k++)
        for (int j=0; j<Nt; j++){
           beta.set(l, k, j, 0) = beta0(k,j) ;  // beta0 must represent the values of the radius. 
    }
    reset_coord() ; 
    
}



}