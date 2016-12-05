/*
 *  Methods of class Map
 *
 *   (see file map.h for documentation)
 *
 */

/*
 *   Copyright (c) 1999-2003 Eric Gourgoulhon
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
 * Revision 1.12  2016/12/05 16:17:56  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.11  2014/10/13 08:53:02  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.10  2014/10/06 15:13:12  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.9  2008/09/29 13:23:51  j_novak
 * Implementation of the angular mapping associated with an affine
 * mapping. Things must be improved to take into account the domain index.
 *
 * Revision 1.8  2004/01/29 08:50:03  p_grandclement
 * Modification of Map::operator==(const Map&) and addition of the surface
 * integrales using Scalar.
 *
 * Revision 1.7  2003/12/30 22:53:23  e_gourgoulhon
 * Added methods flat_met_spher() and flat_met_cart() to get
 * flat metric associated with the coordinates described by the mapping.
 *
 * Revision 1.6  2003/10/15 10:30:46  e_gourgoulhon
 * Map::set_ori: changed x,y,z to xo,yo,zo not to shadow Coord's x,y,z
 * Map::set_rot_phi: changed phi to newphi not to shadow Coord phi.
 *
 * Revision 1.5  2003/06/20 14:45:06  f_limousin
 * Add operator==
 *
 * Revision 1.4  2002/10/16 14:36:40  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.3  2002/05/22 12:44:04  f_limousin
 * Added print of ori_x, ori_y, ori_z and rot_phi in operator<<.
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
 * Revision 2.12  2000/02/11  14:26:54  eric
 * *** empty log message ***
 *
 * Revision 2.11  2000/02/11  13:38:23  eric
 * Ajout de la fonction convert_absolute.
 *
 * Revision 2.10  2000/02/09  13:25:29  eric
 * Mise en conformite avec le nouveau constructeur de Base_vect_spher
 * (passage des arguments ori_x, ori_y et ori_z).
 *
 * Revision 2.9  2000/01/12  12:54:53  eric
 * Ajout du Cmp null, *p_cmp_zero, et de la methode associee cmp_zero().
 *
 * Revision 2.8  2000/01/10  13:28:47  eric
 * Ajout des bases vectorielles associees aux coordonnees :
 *  membres bvect_spher et bvect_cart.
 *
 * Revision 2.7  1999/11/24  14:39:34  eric
 * Appel de reset_coord() dans les fonctions set_ori et set_rot_phi.
 *
 * Revision 2.6  1999/11/22  10:35:36  eric
 * Fonction del_coord() rebaptisee reset_coord().
 *
 * Revision 2.5  1999/10/15  14:12:29  eric
 * Suppression de l'appel a del_coord() dans le destructeur de Map.
 *
 * Revision 2.4  1999/10/15  09:23:04  eric
 * *** empty log message ***
 *
 * Revision 2.3  1999/10/14  14:26:51  eric
 * Depoussierage.
 * Documentation.
 *
 * Revision 2.2  1999/03/01  16:57:13  eric
 * Operateur <<
 *
 * Revision 2.1  1999/03/01  14:57:52  eric
 * *** empty log message ***
 *
 * Revision 2.0  1999/02/15  10:42:45  hyc
 * *** empty log message ***
 *
 * $Header$
 *
 */

// headers C
#include <cmath>
	
// headers Lorene
#include "map.h"
#include "cmp.h"
#include "metric.h"
#include "utilitaires.h"

			//---------------//
			// Constructeurs //
			//---------------//

// Constructor from a grid
// -----------------------
namespace Lorene {
Map::Map(const Mg3d& mgi) : mg(&mgi), 
			    ori_x(0), ori_y(0), ori_z(0), rot_phi(0), 
			    bvect_spher(ori_x, ori_y, ori_z, rot_phi, 
					"Mapping orthonormal spherical basis"), 
			    bvect_cart(rot_phi, "Mapping Cartesian basis"),
                            p_flat_met_spher(0x0), 
                            p_flat_met_cart(0x0),
			    p_mp_angu(0x0)
{
        // The Coord's are constructed by the default Coord constructor
	
	// The null Cmp : 
	p_cmp_zero = new Cmp(this) ; 
	p_cmp_zero->set_etat_zero() ; 
}

// Copy constructor 
// ----------------
Map::Map(const Map& mp) : mg(mp.mg), 
			  ori_x(mp.ori_x), ori_y(mp.ori_y), ori_z(mp.ori_z), 
			  rot_phi(mp.rot_phi),
			  bvect_spher(ori_x, ori_y, ori_z, rot_phi, 
				      "Mapping orthonormal spherical basis"), 
			  bvect_cart(rot_phi, "Mapping Cartesian basis"), 
                          p_flat_met_spher(0x0), 
                          p_flat_met_cart(0x0),
			  p_mp_angu(0x0) 
{
        // The Coord's are constructed by the default Coord constructor

	// The null Cmp : 
	p_cmp_zero = new Cmp(this) ; 
	p_cmp_zero->set_etat_zero() ; 
}
	
// Constructor from file
// ---------------------
Map::Map(const Mg3d& mgi, FILE* fd) : mg(&mgi), 
				      bvect_spher(0., 0., 0., 0., 
					"Mapping orthonormal spherical basis"), 
				      bvect_cart(0., "Mapping Cartesian basis"), 
                                      p_flat_met_spher(0x0), 
                                      p_flat_met_cart(0x0),
				      p_mp_angu(0x0) 
{
    Mg3d* mg_tmp = new Mg3d(fd) ;	// la multi-grille d'origine
    if (*mg != *mg_tmp) {
	cout << "Map::Map(const Mg3d&, FILE*): grid not consistent !" 
	     << endl ;
	abort() ;
    }
    delete mg_tmp ;
    
    fread_be(&ori_x, sizeof(double), 1, fd) ;		// ori_x
    fread_be(&ori_y, sizeof(double), 1, fd) ;		// ori_y
    fread_be(&ori_z, sizeof(double), 1, fd) ;		// ori_z
    fread_be(&rot_phi, sizeof(double), 1, fd) ;		// rot_phi

    bvect_spher.set_ori(ori_x, ori_y, ori_z) ; 
    bvect_spher.set_rot_phi(rot_phi) ; 
    bvect_cart.set_rot_phi(rot_phi) ; 

    // The Coord's are constructed by the default Coord constructor
    
    // The null Cmp : 
    p_cmp_zero = new Cmp(this) ; 
    p_cmp_zero->set_etat_zero() ; 
}

			//--------------//
			// Destructeurs //
			//--------------//

// Destructeur
Map::~Map() {
    if (p_flat_met_spher != 0x0) delete p_flat_met_spher ;
    if (p_flat_met_cart != 0x0) delete p_flat_met_cart ;
    if (p_mp_angu != 0x0) delete p_mp_angu ;
    delete p_cmp_zero ;
}

			//------------//
			// Sauvegarde //
			//------------//

void Map::sauve(FILE* fd) const {

    mg->sauve(fd) ;			    // la multi-grille

    fwrite_be(&ori_x, sizeof(double), 1, fd) ;		// ori_x
    fwrite_be(&ori_y, sizeof(double), 1, fd) ;		// ori_y
    fwrite_be(&ori_z, sizeof(double), 1, fd) ;		// ori_z
    fwrite_be(&rot_phi, sizeof(double), 1, fd) ;		// rot_phi
}
    
			//------------//
			// Impression //
			//------------//

// Operateurs <<
ostream& operator<<(ostream& o, const Map & cv)  {
  o << "Absolute coordinates of the mapping origin: " << endl ;
  o << "  X_0, Y_0, Z_0 : " << cv.get_ori_x() << "  "
    <<  cv.get_ori_y() << "  " << cv.get_ori_z() << endl ; 
  o << "Rotation angle between the x-axis and X-axis : " 
    << cv.get_rot_phi() << endl ;   
    cv >> o ;
    return o ;
}
    
			//-------------------//
			// Methodes diverses //
			//-------------------//

void Map::set_ori(double xo, double yo, double zo) {
    ori_x = xo ;
    ori_y = yo ;
    ori_z = zo ;

    bvect_spher.set_ori(ori_x, ori_y, ori_z) ; 
    
    reset_coord() ;  // Mise a jour des Coords 
}

void Map::set_rot_phi(double newphi) {

    rot_phi = newphi ;

    bvect_spher.set_rot_phi(rot_phi) ; 
    bvect_cart.set_rot_phi(rot_phi) ; 

    reset_coord() ;  // Mise a jour des Coords 

}

// Gestion de la memoire
// ---------------------
void Map::reset_coord() {

 	// Coordonnees spheriques centrees sur la grille : 
	 r.del_t() ;
	 tet.del_t() ;
	 phi.del_t() ;
	 sint.del_t() ;
	 cost.del_t() ;
	 sinp.del_t() ;
	 cosp.del_t() ;

	// Coordonnees cartesiennes centrees sur la grille : 
	 x.del_t() ;
	 y.del_t() ;
	 z.del_t() ;

	// Coordonnees cartesiennes "absolues" : 
	 xa.del_t() ;
	 ya.del_t() ;
	 za.del_t() ;
}


// Conversion coordonnees (X,Y,Z) --> (r, theta, phi)
// --------------------------------------------------

void Map::convert_absolute(double xx, double yy, double zz, 
			   double& rr, double& theta, double& pphi) const {

    // Cartesian coordinates aligned with the absolute ones but centered
    // on the mapping origin : 
    double x1 = xx - ori_x ;
    double y1 = yy - ori_y ;
    double z1 = zz - ori_z ;
 
    // Spherical coordinates : 		    
    double rho2 = x1*x1 + y1*y1 ; 
    double rho = sqrt( rho2 ) ;
    rr = sqrt(rho2 + z1*z1) ;
    theta = atan2(rho, z1) ;
    pphi = atan2(y1, x1) - rot_phi ;		// (rotation)
    if (pphi < 0) pphi += 2*M_PI ;
			 
}

const Metric_flat& Map::flat_met_spher() const {

    if (p_flat_met_spher == 0x0) {
        p_flat_met_spher = new Metric_flat(*this, bvect_spher) ; 
    }
    
    return *p_flat_met_spher ;

} 

const Metric_flat& Map::flat_met_cart() const {

    if (p_flat_met_cart == 0x0) {
        p_flat_met_cart = new Metric_flat(*this, bvect_cart) ; 
    }
    
    return *p_flat_met_cart ;

} 











}
