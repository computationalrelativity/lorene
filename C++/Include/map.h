/*
 *  Definition of Lorene classes Map
 *				 Map_radial
 *				 Map_af
 *				 Map_et
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *   Copyright (c) 2000-2001 Jerome Novak
 *   Copyright (c) 2000-2001 Keisuke Taniguchi
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


#ifndef __MAP_H_ 
#define __MAP_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.7  2003/06/20 09:27:09  j_novak
 * Modif commentaires.
 *
 * Revision 1.6  2002/10/16 14:36:29  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.5  2002/09/13 09:17:33  j_novak
 * Modif. commentaires
 *
 * Revision 1.4  2002/06/17 14:05:16  j_novak
 * friend functions are now also declared outside the class definition
 *
 * Revision 1.3  2002/05/07 07:06:37  e_gourgoulhon
 * Compatibily with xlC compiler on IBM SP2:
 *   added declaration of functions map_af_fait_* and map_et_fait_*
 *   outside the classes declarations.
 *
 * Revision 1.2  2002/01/15 15:53:06  p_grandclement
 * I have had a constructor fot map_et using the equation of the surface
 * of the star.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.110  2001/10/29  15:31:55  novak
 * Ajout de Map_radial::div_r
 *
 * Revision 2.109  2001/10/16 10:02:49  novak
 * *** empty log message ***
 *
 * Revision 2.108  2001/07/19 14:01:00  novak
 * new arguments for Map_af::dalembert
 *
 * Revision 2.107  2001/02/26 17:28:31  eric
 * Ajout de la fonction virtuelle resize.
 *
 * Revision 2.106  2001/01/10  11:03:00  phil
 * ajout de homothetie interne
 *
 * Revision 2.105  2001/01/02  10:51:55  phil
 * ajout integrale de surface a l'infini
 *
 * Revision 2.104  2000/10/23  13:59:48  eric
 * Map_et::adapt: changement des arguments (en autre, ajout de nz_search).
 *
 * Revision 2.103  2000/10/20  09:39:19  phil
 * changement commentaires
 *
 * Revision 2.102  2000/10/19  14:33:23  novak
 * corrige oubli pour Map_et?
 *
 * Revision 2.101  2000/10/19 14:11:20  novak
 * Ajout des fonctions membres Map::dalembert et Map_af::dalembert
 * (etat experimental)
 *
 * Revision 2.100  2000/10/09  13:46:39  eric
 * Ajout de la fonction virtuelle poisson2d.
 *
 * Revision 2.99  2000/09/19  15:28:55  phil
 * *** empty log message ***
 *
 * Revision 2.98  2000/09/19  15:24:19  phil
 * ajout du passage de cartesienne en spheriques
 *
 * Revision 2.97  2000/09/19  13:05:38  phil
 * ajout integrale_surface
 *
 * Revision 2.96  2000/09/11  15:54:03  eric
 * Suppression des methodes deriv_x, deriv_y et deriv_z.
 * Introduction des methodes comp_x_from_spherical, etc...
 *
 * Revision 2.95  2000/09/07  15:27:58  keisuke
 * Add a new argument Cmp& uu in Map_af::poisson_regular and Map_et::poisson_regular.
 *
 * Revision 2.94  2000/09/04  15:30:56  keisuke
 * Modify the arguments of Map_af::poisson_regular and Map_et::poisson_regular.
 *
 * Revision 2.93  2000/09/04  13:36:19  keisuke
 * Modify the explanation for "uu_div" in Map_et::poisson_regular.
 *
 * Revision 2.92  2000/08/31  15:50:12  keisuke
 * Modify Map_af::poisson_regular.
 * Add Map_et::poisson_regular and Map::poisson_regular.
 *
 * Revision 2.91  2000/08/31  13:03:22  eric
 * Ajout de la fonction virtuelle mult_rsint.
 *
 * Revision 2.90  2000/08/28  16:17:37  keisuke
 * Add "int nzet" in the argumant of Map_af::poisson_regular.
 *
 * Revision 2.89  2000/08/18  11:10:12  eric
 * Classe Map_et: ajout de l'operateur d'affectation a un autre Map_et.
 *
 * Revision 2.88  2000/08/11  08:50:18  keisuke
 * Modif Map_af::poisson_regular
 *
 * Revision 2.87  2000/08/10  12:54:00  keisuke
 * Ajout de Map_af::poisson_regular
 *
 * Revision 2.86  2000/07/20  14:21:07  eric
 * Ajout de la fonction div_rsint.
 *
 * Revision 2.85  2000/05/25  13:54:41  eric
 * Modif commentaires
 *
 * Revision 2.84  2000/05/22  14:38:51  phil
 * ajout de inc_dzpuis et dec_dzpuis
 *
 * Revision 2.83  2000/04/27  15:18:54  phil
 * *** empty log message ***
 *
 * Revision 2.82  2000/03/20  13:33:23  phil
 * commentaires
 *
 * Revision 2.81  2000/03/17  17:32:48  phil
 * *** empty log message ***
 *
 * Revision 2.80  2000/03/17  17:01:54  phil
 * *** empty log message ***
 *
 * Revision 2.79  2000/03/17  16:58:48  phil
 * ajout de poisson_frontiere
 *
 * Revision 2.78  2000/03/06  11:29:51  eric
 * Ajout du membre reeavaluate_symy.
 *
 * Revision 2.77  2000/02/15  15:08:21  eric
 * Changement du Param dans Map_et::adapt : fact_echelle est desormais
 * passe en double_mod.
 *
 * Revision 2.76  2000/02/15  10:26:25  phil
 * commentaire +
 * suppression de poisson_vect et poisson_vect_oohara
 *
 * Revision 2.75  2000/02/11  13:37:43  eric
 * Ajout de la fonction convert_absolute.
 *
 * Revision 2.74  2000/02/09  09:53:37  phil
 * ajout de poisson_vect_oohara
 *
 * Revision 2.73  2000/01/26  13:07:02  eric
 * Reprototypage complet des routines de derivation:
 *  le resultat est desormais suppose alloue a l'exterieur de la routine
 *  et est passe en argument (Cmp& resu), si bien que le prototypage
 *  complet devient:
 * 		void DERIV(const Cmp& ci, Cmp& resu)
 *
 * Revision 2.72  2000/01/24  17:08:21  eric
 * Class Map_af : suppression de la fonction convert.
 *                suppression du constructeur par convertion d'un Map_et.
 *                ajout du constructeur par conversion d'un Map.
 *
 * Revision 2.71  2000/01/24  16:41:43  eric
 * Ajout de la fonction virtuelle operator=(const Map_af& ).
 * Classe Map_af : ajout de la fonction convert(const Map& ).
 *
 * Revision 2.70  2000/01/21  12:48:34  phil
 * changement prototypage de Map::poisson_vect
 *
 * Revision 2.69  2000/01/20  16:35:05  phil
 * *** empty log message ***
 *
 * Revision 2.68  2000/01/20  15:44:42  phil
 * *** empty log message ***
 *
 * Revision 2.67  2000/01/20  15:31:56  phil
 * *** empty log message ***
 *
 * Revision 2.66  2000/01/20  14:18:06  phil
 * *** empty log message ***
 *
 * Revision 2.65  2000/01/20  13:16:34  phil
 * *** empty log message ***
 *
 * Revision 2.64  2000/01/20  12:51:24  phil
 * *** empty log message ***
 *
 * Revision 2.63  2000/01/20  12:45:28  phil
 * *** empty log message ***
 *
 * Revision 2.62  2000/01/20  12:40:27  phil
 * *** empty log message ***
 *
 * Revision 2.61  2000/01/20  11:27:54  phil
 * ajout de poisson_vect
 *
 * Revision 2.60  2000/01/13  15:31:55  eric
 * Modif commentaires/
 *
 * Revision 2.59  2000/01/12  16:02:57  eric
 * Modif commentaires poisson_compact.
 *
 * Revision 2.58  2000/01/12  12:54:23  eric
 * Ajout du Cmp null, *p_cmp_zero, et de la methode associee cmp_zero().
 *
 * Revision 2.57  2000/01/10  13:27:43  eric
 * Ajout des bases vectorielles associees aux coordonnees :
 *  membres bvect_spher et bvect_cart.
 *
 * Revision 2.56  2000/01/10  09:12:47  eric
 * Reprototypage de poisson_compact : Valeur -> Cmp, Tenseur.
 * Suppression de poisson_compact_boucle.
 * poisson_compact est desormais implementee au niveau Map_radial.
 *
 * Revision 2.55  2000/01/04  15:23:11  eric
 * Classe Map_radial : les data sont listees en premier
 * Introduction de la fonction reevalutate.
 *
 * Revision 2.54  2000/01/03  13:30:32  eric
 * Ajout de la fonction adapt.
 *
 * Revision 2.53  1999/12/22  17:09:52  eric
 * Modif commentaires.
 *
 * Revision 2.52  1999/12/21  16:26:25  eric
 * Ajout du constructeur par conversion Map_af::Map_af(const Map_et&).
 * Ajout des fonctions Map_af::set_alpha et Map_af::set_beta.
 *
 * Revision 2.51  1999/12/21  13:01:29  eric
 * Changement de prototype de la routine poisson : la solution est
 *  desormais passee en argument (et non plus en valeur de retour)
 *  pour permettre l'initialisation de methodes de resolution iterative.
 *
 * Revision 2.50  1999/12/21  10:12:09  eric
 * Modif commentaires.
 *
 * Revision 2.49  1999/12/21  10:06:05  eric
 * Ajout de l'argument Param& a poisson.
 *
 * Revision 2.48  1999/12/20  15:44:35  eric
 * Modif commentaires.
 *
 * Revision 2.47  1999/12/20  10:47:45  eric
 * Modif commentaires.
 *
 * Revision 2.46  1999/12/20  10:24:12  eric
 * Ajout des fonctions de lecture des parametres de Map_et:
 *   get_alpha(), get_beta(), get_ff(), get_gg().
 *
 * Revision 2.45  1999/12/16  14:50:08  eric
 * Modif commentaires.
 *
 * Revision 2.44  1999/12/16  14:17:54  eric
 * Introduction de l'argument const Param& par dans val_lx et val_lx_jk.
 * (en remplacement de l'argument Tbl& param).
 *
 * Revision 2.43  1999/12/09  10:45:24  eric
 * Ajout de la fonction virtuelle integrale.
 *
 * Revision 2.42  1999/12/07  14:50:47  eric
 * Changement ordre des arguments val_r, val_lx
 * val_r_kj --> val_r_jk
 * val_lx_kj -->val_lx_jk
 *
 * Revision 2.41  1999/12/06  16:45:20  eric
 * Surcharge de val_lx avec la version sans param.
 *
 * Revision 2.40  1999/12/06  15:33:44  eric
 * Ajout des fonctions val_r_kj et val_lx_kj.
 *
 * Revision 2.39  1999/12/06  13:11:54  eric
 * Introduction des fonctions val_r, val_lx et homothetie.
 *
 * Revision 2.38  1999/12/02  14:28:22  eric
 * Reprototypage de la fonction poisson: Valeur -> Cmp.
 *
 * Revision 2.37  1999/11/30  14:19:33  eric
 * Reprototypage complet des fonctions membres mult_r, mult_r_zec,
 *  dec2_dzpuis et inc2_dzpuis : Valeur --> Cmp
 *
 * Revision 2.36  1999/11/29  13:17:57  eric
 * Modif commentaires.
 *
 * Revision 2.35  1999/11/29  12:55:42  eric
 * Changement prototype de la fonction laplacien : Valeur --> Cmp.
 *
 * Revision 2.34  1999/11/25  16:27:27  eric
 * Reorganisation complete du calcul des derivees partielles.
 *
 * Revision 2.33  1999/11/24  16:31:17  eric
 * Map_et: ajout des fonctions set_ff et set_gg.
 *
 * Revision 2.32  1999/11/24  14:31:48  eric
 * Map_af: les membres alpha et beta deviennent prives.
 * Map_af: introduction des fonctions get_alpha() et get_beta().
 *
 * Revision 2.31  1999/11/24  11:22:09  eric
 * Map_et : Coords rendus publics
 * Map_et : fonctions de constructions amies.
 *
 * Revision 2.30  1999/11/22  10:32:39  eric
 * Introduction de la classe Map_et.
 * Constructeurs de Map rendus protected.
 * Fonction del_coord() rebaptisee reset_coord().
 *
 * Revision 2.29  1999/10/27  16:44:41  phil
 * ajout de mult_r_zec
 *
 * Revision 2.28  1999/10/19  14:40:37  phil
 * ajout de inc2_dzpuis()
 *
 * Revision 2.27  1999/10/15  14:12:20  eric
 * *** empty log message ***
 *
 * Revision 2.26  1999/10/14  14:26:06  eric
 * Depoussierage.
 * Documentation.
 *
 * Revision 2.25  1999/10/11  11:16:29  phil
 * changement prototypage de poisson_compact_boucle
 *
 * Revision 2.24  1999/10/11  10:48:51  phil
 * changement de nom pour poisson a support compact
 *
 * Revision 2.23  1999/10/04  09:20:58  phil
 * changement de prototypage de void Map_af::poisson_nul
 *
 * Revision 2.22  1999/09/30  18:38:32  phil
 * *** empty log message ***
 *
 * Revision 2.21  1999/09/30  18:33:10  phil
 * ajout de poisson_nul et poisson_nul_boucle
 *
 * Revision 2.20  1999/09/30  16:45:54  phil
 * ajout de Map_af::poisson_nul(const Valeur&, int, int)
 *
 * Revision 2.19  1999/09/16  13:15:40  phil
 * ajout de Valeur mult_r (const Valeur &)
 *
 * Revision 2.18  1999/09/15  10:42:11  phil
 * ajout de Valeur dec2_dzpuis(const Valeur&)
 *
 * Revision 2.17  1999/09/14  13:45:45  phil
 * suppression de la divergence
 *
 * Revision 2.16  1999/09/13  15:09:07  phil
 * ajout de Map_af::divergence
 *
 * Revision 2.15  1999/09/13  13:52:23  phil
 * ajout des derivations partielles par rapport a x,y et z.
 *
 * Revision 2.14  1999/09/07  14:35:20  phil
 * ajout de la fonction Valeur** gradient(const Valeur&)
 *
 * Revision 2.13  1999/04/26  16:37:43  phil
 * *** empty log message ***
 *
 * Revision 2.12  1999/04/26  16:33:28  phil
 * *** empty log message ***
 *
 * Revision 2.11  1999/04/26  13:53:04  phil
 * *** empty log message ***
 *
 * Revision 2.10  1999/04/26  13:51:19  phil
 * ajout de Map_af::laplacien  (2versions)
 *
 * Revision 2.9  1999/04/14  09:04:01  phil
 * *** empty log message ***
 *
 * Revision 2.8  1999/04/14  08:53:27  phil
 * *** empty log message ***
 *
 * Revision 2.7  1999/04/13  17:45:25  phil
 * *** empty log message ***
 *
 * Revision 2.6  1999/04/13  17:02:41  phil
 * ,
 *
 * Revision 2.5  1999/04/13  17:00:41  phil
 * ajout de la resolution de poisson affine
 *
 * Revision 2.4  1999/03/04  13:10:53  eric
 * Ajout des Coord representant les derivees du changement de variable
 * dans la classe Map_radial.
 *
 * Revision 2.3  1999/03/01  17:00:38  eric
 * *** empty log message ***
 *
 * Revision 2.2  1999/03/01  16:44:41  eric
 * Operateurs << et >> sur les ostream.
 * L'operateur >> est virtuel.
 *
 * Revision 2.1  1999/02/22  15:21:45  hyc
 * *** empty log message ***
 *
 *
 * Revision 2.0  1999/01/15  09:10:39  hyc
 * *** empty log message ***
 *
 * $Header$
 *
 */

#include <stdio.h>

#include "coord.h"
#include "base_vect.h"
#include "valeur.h"

class Cmp ;
class Param ; 
class Map_af ; 
class Map_et ; 
class Tenseur ; 

			//------------------------------------//
			//	    class Map		      //
			//------------------------------------//

/**
 * Base class for coordinate mappings.
 * 
 * This class is the basic class for describing the mapping between the
 * grid coordinates $(\xi, \theta', \phi')$ and the physical coordinates
 * $(r, \theta, \phi)$ [cf. Bonazzola, Gourgoulhon \& Marck, {\sl Phys. Rev. D}
 * {\bf 58}, 104020 (1998)]. 
 * The class {\tt Map} is an abstract one: it cannot be instanciated. 
 * Specific implementation of coordinate mappings will be performed by derived
 * classes of {\tt Map}. 
 *
 * The class {\tt Map} and its derived classes determine the methods for 
 * partial derivatives with respect to the physical coordinates, as well
 * as resolution of basic partial differential equations (e.g. Poisson
 * equations). 
 *
 * The mapping is defined with respect to some ``absolute'' reference frame, 
 * whose Cartesian coordinates are denoted by {\it (X, Y, Z)}. The coordinates
 * (X, Y, Z) of center of the mapping (i.e. the point {\it r}=0) are given by the data 
 * members {\tt  (ori\_x,  ori\_y,  ori\_z)}. 
 * The Cartesian coordinate relative to the mapping (i.e. defined from 
 * $(r, \theta, \phi)$ by the usual formul\ae $x=r\sin\theta\cos\phi, \ldots$)
 * are denoted by {\it (x, y, z)}. The planes {\it (x, y)} and {\it (X, Y)} are supposed to
 * coincide, so that the relative orientation of the mapping with respect to 
 * the absolute reference frame is described by only one angle (the data member
 * {\tt rot\_phi}).
 *
 * @version #$Id$#
 *
 */

class Map {

    // Data : 
    // ----
    protected:
	/// Pointer on the multi-grid {\tt Mgd3} on which {\tt this} is defined	
	const Mg3d* mg ; 
	 
	double ori_x ;		/// Absolute coordinate {\it X} of the origin
	double ori_y ;		/// Absolute coordinate {\it Y} of the origin
	double ori_z ;		/// Absolute coordinate {\it Z} of the origin
	double rot_phi ;	/// Angle between the {\it x}--axis and {\it X}--axis
    
	/** Orthonormal vectorial basis 
	 *  $(\partial/\partial r,1/r\partial/\partial \theta,
	 *	1/(r\sin\theta)\partial/\partial \phi)$
	 * associated with the coordinates $(r, \theta, \phi)$ of the 
	 * mapping.
	 */
	Base_vect_spher bvect_spher ; 

	/** Cartesian basis 
	 *  $(\partial/\partial x,\partial/\partial y,\partial/\partial z)$
	 * associated with the coordinates {\it (x, y, z)} of the
	 * mapping, i.e. the Cartesian coordinates related to 
	 * $(r, \theta, \phi)$ by means of usual formulae. 
	 */
	Base_vect_cart bvect_cart ; 
	
	/** The null Cmp. 
	 *  To be used by the {\tt Tenseur} class when necessary to 
	 *  return a null {\tt Cmp}. 
	 */
	Cmp* p_cmp_zero ; 

    public:
	Coord r ;	/// {\it r} coordinate centered on the grid
	Coord tet ;	/// $\theta$ coordinate centered on the grid
	Coord phi ;	/// $\phi$ coordinate centered on the grid
	Coord sint ;	/// $\sin\theta$
	Coord cost ;	/// $\cos\theta$
	Coord sinp ;	/// $\sin\phi$
	Coord cosp ;	/// $\cos\phi$

	Coord x ;	/// {\it x} coordinate centered on the grid
	Coord y ;	/// {\it y} coordinate centered on the grid
	Coord z ;	/// {\it z} coordinate centered on the grid

	Coord xa ;	/// Absolute {\it X} coordinate
	Coord ya ;	/// Absolute {\it Y} coordinate
	Coord za ;	/// Absolute {\it Z} coordinate
    

    // Constructors, destructor : 
    // ------------------------

    protected:
	explicit Map(const Mg3d& ) ;	/// Constructor from a multi-domain 3D grid
	Map(const Map &) ;		/// Copy constructor
	Map(const Mg3d&, FILE* ) ; /// Constructor from a file (see {\tt sauve(FILE* )})

    public:	
	virtual ~Map() ;		/// Destructor
	
    // Memory management
    // -----------------
    protected:
	virtual void reset_coord() ;  /// Resets all the member {\tt Coord}s	
	
    // Outputs
    // -------
    private:
	virtual ostream& operator>>(ostream &) const = 0 ;    /// Operator >>

    public:
	virtual void sauve(FILE* ) const ;	/// Save in a file
	

    // Extraction of information
    // -------------------------

    public:
	/// Gives the {\tt Mg3d} on which the mapping is defined
	const Mg3d* get_mg() const {return mg; };

	double get_ori_x() const {return ori_x;} ; /// Returns the {\it X} coordinate of the origin 
	double get_ori_y() const {return ori_y;} ; /// Returns the {\it Y} coordinate of the origin
	double get_ori_z() const {return ori_z;} ; /// Returns the {\it Z} coordinate of the origin

	/// Returns the angle between the {\it x}--axis and {\it X}--axis
	double get_rot_phi() const {return rot_phi;} ; 
	
	/** Returns the orthonormal vectorial basis 
	 *  $(\partial/\partial r,1/r\partial/\partial \theta,
	 *	1/(r\sin\theta)\partial/\partial \phi)$
	 * associated with the coordinates $(r, \theta, \phi)$ of the 
	 * mapping.
	 */
	const Base_vect_spher& get_bvect_spher() const {return bvect_spher;} ; 

	/** Returns the Cartesian basis 
	 *  $(\partial/\partial x,\partial/\partial y,\partial/\partial z)$
	 * associated with the coordinates {\it (x, y, z)} of the
	 * mapping, i.e. the Cartesian coordinates related to 
	 * $(r, \theta, \phi)$ by means of usual formulae. 
	 */
	const Base_vect_cart& get_bvect_cart() const {return bvect_cart;} ; 

	/** Returns the null {\tt Cmp} defined on {\tt *this}. 
	 *  To be used by the {\tt Tenseur} class when necessary to 
	 *  return a null {\tt Cmp}. 
	 */
	const Cmp& cmp_zero() const {return *p_cmp_zero;} ;	
	
	/** Determines the coordinates $(r,\theta,\phi)$
	 *  corresponding to given absolute Cartesian coordinates
	 *  {\it (X, Y, Z)}. 
	 *	@param xx [input] value of the coordinate {\it X} (absolute frame)
	 *	@param yy [input] value of the coordinate {\it Y} (absolute frame)
	 *	@param zz [input] value of the coordinate {\it Z} (absolute frame)
	 *	@param rr [output] value of {\it r}
	 *	@param theta [output] value of $\theta$
	 *	@param pphi [output] value of $\phi$
	 */
	void convert_absolute(double xx, double yy, double zz, 
			      double& rr, double& theta, double& pphi) const ; 

	/** Returns the value of the radial coordinate {\it r} for a given
	 *  $(\xi, \theta', \phi')$ in a given domain. 
	 *	@param l [input] index of the domain
	 *	@param xi [input] value of $\xi$
	 *	@param theta [input] value of $\theta'$
	 *	@param pphi [input] value of $\phi'$
	 *	@return value of $r=R_l(\xi, \theta', \phi')$
	 */
	virtual double val_r(int l, double xi, double theta, double pphi) 
			     const = 0 ; 
	
	/** Computes the domain index {\it l} and the value of $\xi$ corresponding
	 * to a point given by its physical coordinates $(r, \theta, \phi)$. 
	 *	@param rr [input] value of {\it r}
	 *	@param theta [input] value of $\theta$
	 *	@param pphi [input] value of $\phi$
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of $\xi$
	 */
	virtual void val_lx(double rr, double theta, double pphi, 
			    int& l, double& xi) const = 0 ; 
	
	/** Computes the domain index {\it l} and the value of $\xi$ corresponding
	 * to a point given by its physical coordinates $(r, \theta, \phi)$. 
	 * This version enables to pass some parameters to control the
	 * accuracy of the computation. 
	 *	@param rr [input] value of {\it r}
	 *	@param theta [input] value of $\theta$
	 *	@param pphi [input] value of $\phi$
	 *	@param par [input/output] parameters to control the
	 *	    accuracy of the computation
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of $\xi$
	 */
	virtual void val_lx(double rr, double theta, double pphi, 
			    const Param& par, int& l, double& xi) const = 0 ; 
	
    // Modification of the origin, the orientation and the radial scale:
    // ----------------------------------------------------------------
    public:
	void set_ori(double xa0, double ya0, double za0) ;  /// Sets a new origin
	void set_rot_phi(double phi0) ;		/// Sets a new rotation angle

	/** Sets a new radial scale.
	 *	@param lambda [input] factor by which the value of {\it r} is to
	 *	    be multiplied
	 */
	virtual void homothetie(double lambda) = 0 ;	

	/** Rescales the outer boundary of one domain.
	 *  The inner boundary is unchanged. The inner boundary 
	 *  of the next domain is changed to match the new outer
	 *  boundary. 
	 *	@param l [input] index of the domain
	 *	@param lambda [input] factor by which the value of 
	 *	    $R(\theta, \varphi)$ defining the outer boundary
	 *	    of the domain is to be multiplied. 
	 */
	virtual void resize(int l, double lambda) = 0 ; 

    // Modification of the mapping
    // ---------------------------
    public:
	/// Assignment to an affine mapping. 
	virtual void operator=(const Map_af& ) = 0 ;

	/** Adaptation of the mapping to a given scalar field.
	 *  This is a virtual function: see the actual implementations
	 *  in the derived classes for the meaning of the various 
	 *  parameters.
	 */
	virtual void adapt(const Cmp& ent, const Param& par) = 0 ; 

    // Values of a Cmp at the new grid points
    // --------------------------------------

	/** Recomputes the values of a {\tt Cmp} at the collocation points 
	 *  after a change in the mapping. 
	 *  @param mp_prev [input] Previous value of the mapping. 
	 *  @param nzet [input] Number of domains where the computation 
	 *    must be done: the computation is performed for the domains
	 *    of index $0\le {\tt l} \le {\tt nzet-1}$; {\tt uu} is set
	 *    to zero in the other domains. 
	 *  @param uu [input/output] input : {\tt Cmp} previously computed on 
	 *			     the mapping {\tt *mp\_prev}; ouput :
	 *			     values of (logically) the same {\tt Cmp}
	 *			     at the grid points defined by {\tt *this}. 
	 */
	virtual void reevaluate(const Map* mp_prev, int nzet, 
				Cmp& uu) const = 0 ; 

	/** Recomputes the values of a {\tt Cmp} at the collocation points 
	 *  after a change in the mapping. 
	 *  Case where the {\tt Cmp} is symmetric with respect the plane y=0.
	 *  @param mp_prev [input] Previous value of the mapping. 
	 *  @param nzet [input] Number of domains where the computation 
	 *    must be done: the computation is performed for the domains
	 *    of index $0\le {\tt l} \le {\tt nzet-1}$; {\tt uu} is set
	 *    to zero in the other domains. 
	 *  @param uu [input/output] input : {\tt Cmp} previously computed on 
	 *			     the mapping {\tt *mp\_prev}; ouput :
	 *			     values of (logically) the same {\tt Cmp}
	 *			     at the grid points defined by {\tt *this}. 
	 */
	virtual void reevaluate_symy(const Map* mp_prev, int nzet, 
				Cmp& uu) const = 0 ; 



    // Differential operators:
    // ----------------------
    public:
	/** Computes $\partial/ \partial r$ of a {\tt Cmp}.
	 *  Note that in the external compactified domain (ZEC), it computes
	 *  $-\partial/ \partial u = r^2 \partial/ \partial r$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of {\tt ci}
	 */
	virtual void dsdr(const Cmp& ci, Cmp& resu) const = 0 ;  	    

	/** Computes $1/r \partial/ \partial \theta$ of a {\tt Cmp}.
	 *  Note that in the external compactified domain (ZEC), it computes
	 *  $1/u \partial/ \partial \theta = r \partial/ \partial \theta$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of {\tt ci}
	 */
	virtual void srdsdt(const Cmp& ci, Cmp& resu) const = 0 ;  	    
	
	/** Computes $1/(r\sin\theta) \partial/ \partial \phi$ of a {\tt Cmp}.
	 *  Note that in the external compactified domain (ZEC), it computes
	 *  $1/(u\sin\theta) \partial/ \partial \phi = 
	 *	    r/\sin\theta \partial/ \partial \phi$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of {\tt ci}
	 */
	virtual void srstdsdp(const Cmp& ci, Cmp& resu) const = 0 ;  	    
    
	/** Computes the Laplacian of a scalar field.
	 *   @param uu	[input]  Scalar field {\it u} (represented as a {\tt Cmp})
	 *			 the Laplacian $\Delta u$ of which is to be computed
	 *   @param zec_mult_r [input] Determines the quantity computed in
	 *			 the external compactified domain (ZEC) :  \\
	 *		    zec\_mult\_r = 0 : $\Delta u$	\\
	 *		    zec\_mult\_r = 2 : $r^2 \,  \Delta u$	\\
	 *		    zec\_mult\_r = 4 (default) : $r^4 \, \Delta u$	
	 *   @param lap [output] Laplacian of {\it u}
	 */
	virtual void laplacien(const Cmp& uu, int zec_mult_r, 
			       Cmp& lap) const = 0 ; 
	
	 
    // Various linear operators
    // ------------------------
    public: 
	/** Multiplication by {\it r} of a {\tt Cmp}
	 */
	virtual void mult_r(Cmp& ) const = 0 ; 

	/** Multiplication by {\it r} (in the external compactified domain only)
	 * of a {\tt Cmp}
	 */
	virtual void mult_r_zec(Cmp& ) const = 0 ;

	/** Multiplication by $r\sin\theta$ of a {\tt Cmp}
	 */
	virtual void mult_rsint(Cmp& ) const = 0 ; 

	/** Division by $r\sin\theta$ of a {\tt Cmp}
	 */
	virtual void div_rsint(Cmp& ) const = 0 ; 

	/** Division by {\it r} of a {\tt Cmp}
	 */
	virtual void div_r(Cmp& ) const = 0 ; 

	/** Computes the Cartesian x component (with respect to 
	 *  {\tt bvect\_cart}) of a vector given
	 *  by its spherical components with respect to {\tt bvect\_spher}.
	 *  
	 *  @param v_r [input] {\it r}-component of the vector 
	 *  @param v_theta [input] $\theta$-component of the vector 
	 *  @param v_phi [input] $\phi$-component of the vector 
	 *  @param v_x [output] x-component of the vector 
	 */
	virtual void comp_x_from_spherical(const Cmp& v_r, const Cmp& v_theta, 
					   const Cmp& v_phi, Cmp& v_x) const = 0 ; 
	
	/** Computes the Cartesian y component (with respect to 
	 *  {\tt bvect\_cart}) of a vector given
	 *  by its spherical components with respect to {\tt bvect\_spher}.
	 *  
	 *  @param v_r [input] {\it r}-component of the vector 
	 *  @param v_theta [input] $\theta$-component of the vector 
	 *  @param v_phi [input] $\phi$-component of the vector 
	 *  @param v_y [output] y-component of the vector 
	 */
	virtual void comp_y_from_spherical(const Cmp& v_r, const Cmp& v_theta, 
					   const Cmp& v_phi, Cmp& v_y) const = 0 ; 
	
	/** Computes the Cartesian z component (with respect to 
	 *  {\tt bvect\_cart}) of a vector given
	 *  by its spherical components with respect to {\tt bvect\_spher}.
	 *  
	 *  @param v_r [input] {\it r}-component of the vector 
	 *  @param v_theta [input] $\theta$-component of the vector 
	 *  @param v_z [output] z-component of the vector 
	 */
	virtual void comp_z_from_spherical(const Cmp& v_r, const Cmp& v_theta, 
					   Cmp& v_z) const = 0 ; 
	
	 /** Computes the Spherical r component (with respect to 
	 *  {\tt bvect\_spher}) of a vector given
	 *  by its cartesian components with respect to {\tt bvect\_cart}.
	 *  
	 *  @param v_x [input] x-component of the vector 
	 *  @param v_y [input] y-component of the vector 
	 *  @param v_z [input] z-component of the vector 
	 *  @param v_r [output] {\it r}-component of the vector 
	 */
	virtual void comp_r_from_cartesian(const Cmp& v_x, const Cmp& v_y, 
					   const Cmp& v_z, Cmp& v_r) const = 0 ; 
	
	/** Computes the Spherical $\theta$ component (with respect to 
	 *  {\tt bvect\_spher}) of a vector given
	 *  by its cartesian components with respect to {\tt bvect\_cart}.
	 *  
	 *  @param v_x [input] x-component of the vector 
	 *  @param v_y [input] y-component of the vector 
	 *  @param v_z [input] z-component of the vector 
	 *  @param v_t [output] $\theta$-component of the vector 
	 */
	virtual void comp_t_from_cartesian(const Cmp& v_x, const Cmp& v_y, 
					   const Cmp& v_z, Cmp& v_t) const = 0 ; 
	
	/** Computes the Spherical $\phi$ component (with respect to 
	 *  {\tt bvect\_spher}) of a vector given
	 *  by its cartesian components with respect to {\tt bvect\_cart}.
	 *  
	 *  @param v_x [input] x-component of the vector 
	 *  @param v_y [input] y-component of the vector 
	 *  @param v_p [output] $\phi$-component of the vector 
	 */
	virtual void comp_p_from_cartesian(const Cmp& v_x, const Cmp& v_y, 
					   Cmp& v_p) const = 0 ; 
	
	/** Decreases by 1 the value of {\tt dzpuis} of a {\tt Cmp} 
	 *  and changes accordingly its values in the external 
	 *  compactified domain (ZEC).
	 */
	virtual void dec_dzpuis(Cmp& ) const = 0 ; 
	
	/** Decreases by 2 the value of {\tt dzpuis} of a {\tt Cmp} 
	 *  and changes accordingly its values in the external 
	 *  compactified domain (ZEC).
	 */
	virtual void dec2_dzpuis(Cmp& ) const = 0 ; 
	
	/** Increases by 1 the value of {\tt dzpuis} of a {\tt Cmp} 
	 *  and changes accordingly its values in the external 
	 *  compactified domain (ZEC).
	 */
	virtual void inc_dzpuis(Cmp& ) const = 0 ; 
	
	/** Increases by 2 the value of {\tt dzpuis} of a {\tt Cmp} 
	 *  and changes accordingly its values in the external 
	 *  compactified domain (ZEC).
	 */
	virtual void inc2_dzpuis(Cmp& ) const = 0 ; 
	
	/** Computes the integral over all space of a {\tt Cmp}.
	 *  The computed quantity is 
	 *    $\int u \, r^2 \sin\theta \,  dr\, d\theta \, d\phi$.
	 *  The routine allocates a {\tt Tbl} (size: {\tt mg->nzone}) to store 
	 *  the result (partial integral) in each domain and returns a pointer 
	 *  to it.
	 */
	virtual Tbl* integrale(const Cmp&) const = 0 ; 
	
    // PDE resolution :
    // --------------
    public:
	/** Computes the solution of a scalar Poisson equation.
	 * 
	 *   @param source [input] source $\sigma$ of the Poisson equation 
	 *	    $\Delta u = \sigma$.
	 *   @param par [input/output] possible parameters to control the
	 *   resolution of the Poisson equation. See the actual implementation 
	 *   in the derived class of {\tt Map} for documentation. 
	 *   @param uu [input/output] solution {\it u} with the boundary condition 
	 *	    {\it u}=0 at spatial infinity. 
	 */
	virtual void poisson(const Cmp& source, Param& par, Cmp& uu) const = 0 ;

	/** Computes the solution of a scalar Poisson equation.
	 *   The regularized source
	 *
	 *   @param source [input] source $\sigma$ of the Poisson equation 
	 *	    $\Delta u = \sigma$.
	 *   @param k_div [input] regularization degree of the procedure
	 *   @param nzet [input] number of domains covering the star
	 *   @param unsgam1 [input] parameter $1/(\gamma-1)$ where $\gamma$
         *          denotes the adiabatic index.
	 *   @param par [input/output] possible parameters to control the
	 *   resolution of the Poisson equation. See the actual implementation 
	 *   in the derived class of {\tt Map} for documentation.
	 *   @param uu [input/output] solution {\it u} with the boundary condition 
	 *	    {\it u}=0 at spatial infinity.
	 *   @param uu_regu [output] solution of the regular part of
	 *          the source.
         *   @param uu_div [output] solution of the diverging part of
         *          the source.
         *   @param duu_div [output] derivative of the diverging potential
         *   @param source_regu [output] regularized source
         *   @param source_div [output] diverging part of the source
	 */
	virtual void poisson_regular(const Cmp& source, int k_div, int nzet,
				     double unsgam1, Param& par, Cmp& uu,
				     Cmp& uu_regu, Cmp& uu_div,
				     Tenseur& duu_div, Cmp& source_regu,
				     Cmp& source_div) const = 0 ;

	/** Resolution of the elliptic equation 
	 * $ a \Delta\psi + {\bf b} \cdot \nabla \psi = \sigma$.
	 * 
	 * @param source [input] source $\sigma$ of the above equation
	 * @param aa [input] factor {\it a} in the above equation
	 * @param bb [input] vector {\it \bf b} in the above equation
	 * @param par [input/output] possible parameters to control the
	 *   resolution of the equation. See the actual implementation 
	 *   in the derived class of {\tt Map} for documentation. 
	 * @param psi [input/output] solution $\psi$ which satisfies
	 *	$\psi(0)=0$.  
	 */ 
	virtual void poisson_compact(const Cmp& source, const Cmp& aa, 
				     const Tenseur& bb, const Param& par, 
				     Cmp& psi) const = 0 ;
	
    public:
	/**
	 * Function intended to be used by {\tt Map::poisson\_vect}
	 * and {\tt Map::poisson\_vect\_oohara}. It constructs the sets of 
	 * parameters used for each scalar Poisson equation from the one for 
	 * the vectorial one.
	 * 
	 * @param para [input] : the {\tt Param} used for the resolution of 
	 * the vectorial Poisson equation : \\
	 * {\tt para.int()} maximum number of iteration.\\
	 * {\tt para.double(0)} relaxation parameter.\\
	 * {\tt para.double(1)} required precision. \\
	 * {\tt para.tenseur\_mod()} source of the vectorial part at the previous 
	 * step.\\
	 * {\tt para.cmp\_mod()} source of the scalar part at the previous 
	 * step.
	 * 
	 * @param i [input] number of the scalar Poisson equation that is being 
	 * solved (values from 0 to 2 for the componants of the vectorial part
	 * and 3 for the scalar one).
	 * 
	 * @return the pointer on the parameter set used for solving the
	 * scalar Poisson equation labelled by {\it i}.
	 */
	virtual Param* donne_para_poisson_vect (Param& para, int i) const = 0;
	
	/**
	 * Computes the solution of a Poisson equation from the domain 
	 * {\tt num\_front+1}.
	 * imposing a boundary condition at the boundary between the domains 
	 * {\tt num\_front} and {\tt num\_front+1}.
	 * 
	 * @param source [input] : source of the equation.
	 * @param limite [input] : {\tt limite[num\_front]} contains the angular 
	 * function being the boudary condition.
	 * @param raccord [input] : 1 for the Dirichlet problem and 2 for 
	 * the Neumann one.
	 * @param num_front [input] : index of the boudary at which the boundary 
	 * condition has to be imposed.
	 * @param pot [output] : result.
	 */
	virtual void poisson_frontiere (const Cmp& source, const Valeur& limite,
			 int raccord, int num_front, Cmp& pot) const = 0 ;
	
	virtual void poisson_frontiere_double (const Cmp& source, const Valeur& lim_func,
			const Valeur& lim_der, int num_zone, Cmp& pot) const = 0 ;
	
	/** Computes the solution of a 2-D Poisson equation.
	 *  The 2-D Poisson equation writes
	 *  \begin{equation}
	 *	{\partial^2 u\over\partial r^2} + 
	 *	    {1\over r} {\partial u \over \partial r} + 
	 *	    {1\over r^2} {\partial^2 u\over\partial \theta^2} = 
	 *		\sigma \ . 
	 *  \end{equation} 
	 *
	 *   @param source_mat [input] Compactly supported part of 
	 *	    the source $\sigma$ of the 2-D Poisson equation (typically
	 *	    matter terms)
	 *   @param source_quad [input] Non-compactly supported part of 
	 *	    the source $\sigma$ of the 2-D Poisson equation (typically
	 *	    quadratic terms)
	 *   @param par [input/output] possible parameters to control the
	 *   resolution of the Poisson equation. See the actual implementation 
	 *   in the derived class of {\tt Map} for documentation. 
	 *   @param uu [input/output] solution {\it u} with the boundary condition 
	 *	    {\it u}=0 at spatial infinity. 
	 */
	virtual void poisson2d(const Cmp& source_mat, const Cmp& source_quad, 
			       Param& par, Cmp& uu) const = 0 ;

	/** Performs one time-step integration of the d'Alembert scalar equation
	 *   @param par [input/output] possible parameters to control the
	 *   resolution of the Poisson equation. See the actual implementation 
	 *   in the derived class of {\tt Map} for documentation. Note that, 
	 *   at least, param must contain the time step as first {\tt double} 
	 *   parameter.
	 *   @param fJp1 [output] solution $f^{J+1}$ at time {\it J+1}
	 *   with boundary conditions of outgoing radiation (not exact!)
	 *   @param fJ [input] solution $f^J$ at time {\it J}
	 *   @param fJm1 [input] solution $f^{J-1}$ at time {\it J-1}
	 *   @param source [input] source $\sigma$ of the d'Alembert equation 
	 *	    $\diamond u = \sigma$.
	 */
	virtual void dalembert(Param& par, Cmp& fJp1, const Cmp& fJ, 
			       const Cmp& fJm1, const Cmp& source) const = 0 ;

    // Friend functions : 
    // ----------------
	friend ostream& operator<<(ostream& , const Map& ) ;	/// Operator <<
};
ostream& operator<<(ostream& , const Map& ) ; 



			//------------------------------------//
			//	    class Map_radial	      //
			//------------------------------------//



/**
 * Base class for pure radial mappings.
 * 
 * A pure radial mapping is a mapping of the type $r=R(\xi, \theta', \phi')$, 
 * $\theta=\theta'$, $\phi=\phi'$. 
 * The class {\tt Map\_radial} is an abstract class. Effective implementations 
 * of radial mapping are performed by the derived class {\tt Map\_af} and 
 * {\tt Map\_et}.
 * 
 * 
 * @version #$Id$#
 */

class Map_radial : public Map {

    // Data : 
    // ----

    // 0th order derivatives of the mapping
    // - - - - - - - - - - - - - - - - - - 
    public:
	/**
	 * $\xi/R$ in the nucleus; \\
	 * {\it 1/R} in the non-compactified shells; \\
	 * $(\xi-1)/U$ in the compactified outer domain.
	 */
	Coord xsr ;	    
	
    // 1st order derivatives of the mapping
    // - - - - - - - - - - - - - - - - - - 
    public:
	/**
	 * $1/(\partial R/\partial\xi) = \partial x /\partial r$ in the nucleus
	 *   and in the non-compactified shells; \\
	 * $-1/(\partial U/\partial\xi) = \partial x /\partial u$ in the 
	 *   compactified outer domain.
	 */
	Coord dxdr ;	    
	
	/**
	 * $1/R \times (\partial R/\partial\theta')$ in the nucleus
	 *   and in the non-compactified shells; \\
	 * $-1/U \times (\partial U/\partial\theta)$ in the 
	 *   compactified outer domain.
	 */
	Coord srdrdt ;	    

	/**
	 * $1/(R\sin\theta) \times (\partial R/\partial\varphi')$ in the nucleus
	 *   and in the non-compactified shells; \\
	 * $-1/(U\sin\theta) \times (\partial U/\partial\varphi')$ in the 
	 *   compactified outer domain.
	 */
	Coord srstdrdp ;	    

	/**
	 * $1/R^2 \times (\partial R/\partial\theta')$ in the nucleus
	 *   and in the non-compactified shells; \\
	 * $-1/U^2 \times (\partial U/\partial\theta')$ in the 
	 *   compactified outer domain.
	 */
	Coord sr2drdt ;	    
	
	/**
	 * $1/(R^2\sin\theta) \times (\partial R/\partial\varphi')$ in the nucleus
	 *   and in the non-compactified shells; \\
	 * $-1/(U^2\sin\theta) \times (\partial U/\partial\varphi')$ in the 
	 *   compactified outer domain.
	 */
	Coord sr2stdrdp ;   

    // 2nd order derivatives of the mapping
    // - - - - - - - - - - - - - - - - - - 
    public:
	/**
	 * $\partial^2 R/\partial\xi^2$ in the nucleus
	 *   and in the non-compactified shells; \\
	 * $-\partial^2 U/\partial\xi^2 $ in the 
	 *   compactified outer domain.
	 */
	Coord d2rdx2 ;	    
	
	/**
	 * $1/R^2 \times [ 1/\sin(\theta)\times \partial /\partial\theta'
	 *   (\sin\theta \partial R /\partial\theta') + 1/\sin^2\theta
	 *   \partial^2 R /\partial{\varphi'}^2] $ in the nucleus
	 *   and in the non-compactified shells; \\
	 * $- 1/U^2 \times [ 1/\sin(\theta)\times \partial /\partial\theta'
	 *   (\sin\theta \partial U /\partial\theta') + 1/\sin^2\theta
	 *   \partial^2 U /\partial{\varphi'}^2] $ in the 
	 *   compactified outer domain.
	 */
	Coord lapr_tp ;	    
	

	/**
	 * $\partial^2 R/\partial\xi\partial\theta'$ in the nucleus
	 *   and in the non-compactified shells; \\
	 * $-\partial^2 U/\partial\xi\partial\theta'$ in the 
	 *   compactified outer domain.
	 */
	Coord d2rdtdx ;	    
	
	/**
	 * $1/\sin\theta \times \partial^2 R/\partial\xi\partial\varphi'$ 
	 * in the nucleus and in the non-compactified shells; \\
	 * $-1/\sin\theta \times \partial^2 U/\partial\xi\partial\varphi' $ 
	 * in the compactified outer domain.
	 */
	Coord sstd2rdpdx ;  


	/**
	 * $1/R^2 \partial^2 R/\partial{\theta'}^2$ in the nucleus
	 *   and in the non-compactified shells; \\
	 * $-1/U^2 \partial^2 U/\partial{\theta'}^2$ in the 
	 *   compactified outer domain.
	 */
	Coord sr2d2rdt2 ;
	

    // Constructors, destructor : 
    // ------------------------

    protected:
	/// Constructor from a grid (protected to make {\tt Map\_radial} an abstract class) 
	Map_radial(const Mg3d& ) ;		    
	Map_radial(const Map_radial& ) ;   /// Copy constructor
	Map_radial (const Mg3d&, FILE* ) ; /// Constructor from a file (see {\tt sauve(FILE* )})

    public: 
	virtual ~Map_radial() ;		    /// Destructor

    // Memory management
    // -----------------
    protected:
	virtual void reset_coord() ;  /// Resets all the member {\tt Coord}s	
    // Modification of the mapping
    // ---------------------------
    public:
	/// Assignment to an affine mapping. 
	virtual void operator=(const Map_af& ) = 0 ;
		
    // Outputs
    // -------
    public:  
 	virtual void sauve(FILE* ) const ;	/// Save in a file
   
    // Extraction of information
    // -------------------------
	/** Returns the value of the radial coordinate {\it r} for a given
	 *  $\xi$ and a given collocation point in $(\theta', \phi')$ 
	 *   in a given domain. 
	 *	@param l [input] index of the domain
	 *	@param xi [input] value of $\xi$
	 *	@param j [input] index of the collocation point in $\theta'$
	 *	@param k [input] index of the collocation point in $\phi'$
	 *	@return value of $r=R_l(\xi, {\theta'}_j, {\phi'}_k)$
	 */
	virtual double val_r_jk(int l, double xi, int j, int k) const = 0 ; 
	
	/** Computes the domain index {\it l} and the value of $\xi$ corresponding
	 * to a point of arbitrary {\it r} but collocation values of $(\theta, \phi)$ 
	 *	@param rr [input] value of {\it r}
	 *	@param j [input] index of the collocation point in $\theta$
	 *	@param k [input] index of the collocation point in $\phi$
	 *	@param par [input/output] parameters to control the
	 *	    accuracy of the computation
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of $\xi$
	 */
	virtual void val_lx_jk(double rr, int j, int k, const Param& par, 
			       int& l, double& xi) const = 0 ; 

    // Values of a Cmp at the new grid points
    // --------------------------------------
	/** Recomputes the values of a {\tt Cmp} at the collocation points 
	 *  after a change in the mapping. 
	 *  @param mp_prev [input] Previous value of the mapping. 
	 *  @param nzet [input] Number of domains where the computation 
	 *    must be done: the computation is performed for the domains
	 *    of index $0\le {\tt l} \le {\tt nzet-1}$; {\tt uu} is set
	 *    to zero in the other domains. 
	 *  @param uu [input/output] input : {\tt Cmp} previously computed on 
	 *			     the mapping {\tt *mp\_prev}; ouput :
	 *			     values of (logically) the same {\tt Cmp}
	 *			     at the grid points defined by {\tt *this}. 
	 */
	virtual void reevaluate(const Map* mp_prev, int nzet, Cmp& uu) const ; 

	/** Recomputes the values of a {\tt Cmp} at the collocation points 
	 *  after a change in the mapping. 
	 *  Case where the {\tt Cmp} is symmetric with respect to the plane y=0.
	 *  @param mp_prev [input] Previous value of the mapping. 
	 *  @param nzet [input] Number of domains where the computation 
	 *    must be done: the computation is performed for the domains
	 *    of index $0\le {\tt l} \le {\tt nzet-1}$; {\tt uu} is set
	 *    to zero in the other domains. 
	 *  @param uu [input/output] input : {\tt Cmp} previously computed on 
	 *			     the mapping {\tt *mp\_prev}; ouput :
	 *			     values of (logically) the same {\tt Cmp}
	 *			     at the grid points defined by {\tt *this}. 
	 */
	virtual void reevaluate_symy(const Map* mp_prev, int nzet, Cmp& uu) 
				    const ; 

    
    // Various linear operators
    // ------------------------
    public: 
	/**
	 * Multiplication by {\it r} of a {\tt Cmp}
	 */
	virtual void mult_r(Cmp& ) const ; 

	/**
	 * Multiplication by {\it r} (in the external compactified domain only)
	 * of a {\tt Cmp}
	 */
	virtual void mult_r_zec(Cmp& ) const ;

	/** Multiplication by $r\sin\theta$ of a {\tt Cmp}
	 */
	virtual void mult_rsint(Cmp& ) const ; 

	/** Division by $r\sin\theta$ of a {\tt Cmp}
	 */
	virtual void div_rsint(Cmp& ) const ; 

	/** Division by {\it r} of a {\tt Cmp}
	 */
	virtual void div_r(Cmp& ) const ; 

	/** Computes the Cartesian x component (with respect to 
	 *  {\tt bvect\_cart}) of a vector given
	 *  by its spherical components with respect to {\tt bvect\_spher}.
	 *  
	 *  @param v_r [input] {\it r}-component of the vector 
	 *  @param v_theta [input] $\theta$-component of the vector 
	 *  @param v_phi [input] $\phi$-component of the vector 
	 *  @param v_x [output] x-component of the vector 
	 */
	virtual void comp_x_from_spherical(const Cmp& v_r, const Cmp& v_theta, 
					   const Cmp& v_phi, Cmp& v_x) const ; 
	
	/** Computes the Cartesian y component (with respect to 
	 *  {\tt bvect\_cart}) of a vector given
	 *  by its spherical components with respect to {\tt bvect\_spher}.
	 *  
	 *  @param v_r [input] {\it r}-component of the vector 
	 *  @param v_theta [input] $\theta$-component of the vector 
	 *  @param v_phi [input] $\phi$-component of the vector 
	 *  @param v_y [output] y-component of the vector 
	 */
	virtual void comp_y_from_spherical(const Cmp& v_r, const Cmp& v_theta, 
					   const Cmp& v_phi, Cmp& v_y) const ; 
	
	/** Computes the Cartesian z component (with respect to 
	 *  {\tt bvect\_cart}) of a vector given
	 *  by its spherical components with respect to {\tt bvect\_spher}.
	 *  
	 *  @param v_r [input] {\it r}-component of the vector 
	 *  @param v_theta [input] $\theta$-component of the vector 
	 *  @param v_z [output] z-component of the vector 
	 */
	virtual void comp_z_from_spherical(const Cmp& v_r, const Cmp& v_theta, 
					   Cmp& v_z) const ; 
	
	/** Computes the Spherical r component (with respect to 
	 *  {\tt bvect\_spher}) of a vector given
	 *  by its cartesian components with respect to {\tt bvect\_cart}.
	 *  
	 *  @param v_x [input] x-component of the vector 
	 *  @param v_y [input] y-component of the vector 
	 *  @param v_z [input] z-component of the vector 
	 *  @param v_r [output] {\it r}-component of the vector 
	 */
	virtual void comp_r_from_cartesian(const Cmp& v_x, const Cmp& v_y, 
					   const Cmp& v_z, Cmp& v_r) const ; 
	
	/** Computes the Spherical $\theta$ component (with respect to 
	 *  {\tt bvect\_spher}) of a vector given
	 *  by its cartesian components with respect to {\tt bvect\_cart}.
	 *  
	 *  @param v_x [input] x-component of the vector 
	 *  @param v_y [input] y-component of the vector 
	 *  @param v_z [input] z-component of the vector 
	 *  @param v_t [output] $\theta$-component of the vector 
	 */
	virtual void comp_t_from_cartesian(const Cmp& v_x, const Cmp& v_y, 
					   const Cmp& v_z, Cmp& v_t) const ; 
	
	/** Computes the Spherical $\phi$ component (with respect to 
	 *  {\tt bvect\_spher}) of a vector given
	 *  by its cartesian components with respect to {\tt bvect\_cart}.
	 *  
	 *  @param v_x [input] x-component of the vector 
	 *  @param v_y [input] y-component of the vector 
	 *  @param v_p [output] $\phi$-component of the vector 
	 */
	virtual void comp_p_from_cartesian(const Cmp& v_x, const Cmp& v_y, 
					   Cmp& v_p) const ; 
	
	/**
	 * Decreases by 1 the value of {\tt dzpuis} of a {\tt Cmp} 
	 *  and changes accordingly its values in the external 
	 *  compactified domain (ZEC).
	 */
	virtual void dec_dzpuis(Cmp& ) const ; 

	/**
	 * Decreases by 2 the value of {\tt dzpuis} of a {\tt Cmp} 
	 *  and changes accordingly its values in the external 
	 *  compactified domain (ZEC).
	 */
	virtual void dec2_dzpuis(Cmp& ) const ; 

	/**
	 * Increases by 1 the value of {\tt dzpuis} of a {\tt Cmp} 
	 *  and changes accordingly its values in the external 
	 *  compactified domain (ZEC).
	 */
	virtual void inc_dzpuis(Cmp& ) const ; 
	
	/**
	 * Increases by 2 the value of {\tt dzpuis} of a {\tt Cmp} 
	 *  and changes accordingly its values in the external 
	 *  compactified domain (ZEC).
	 */
	virtual void inc2_dzpuis(Cmp& ) const ; 
	

    // PDE resolution :
    // --------------
    public:
	/** Resolution of the elliptic equation 
	 * $ a \Delta\psi + {\bf b} \cdot \nabla \psi = \sigma$.
	 * 
	 * @param source [input] source $\sigma$ of the above equation
	 * @param aa [input] factor {\it a} in the above equation
	 * @param bb [input] vector {\it \bf b} in the above equation
	 * @param par [input/output] parameters of the iterative method of
	 *  resolution : \\
	 *  {\tt par.get\_int(0)} : [input] maximum number of iterations \\
	 *  {\tt par.get\_double(0)} : [input] required precision: the iterative
	 *	  method is stopped as soon as the relative difference between
	 *	 $\psi^J$ and $\psi^{J-1}$ is greater than 
	 *	 {\tt par.get\_double(0)}\\
	 *  {\tt par.get\_double(1)} : [input] relaxation parameter $\lambda$ \\
	 *  {\tt par.get\_int\_mod(0)} : [output] number of iterations 
	 *				    actually used to get the solution.
	 * @param psi [input/output]: input : previously computed value of $\psi$
	 *	to start the iteration (if nothing is known a priori, 
	 *	{\tt psi} must be set to zero); 
	 *	output: solution $\psi$ which satisfies $\psi(0)=0$.  
	 */ 
	virtual void poisson_compact(const Cmp& source, const Cmp& aa, 
				     const Tenseur& bb, const Param& par, 
				     Cmp& psi) const ;

};


			//------------------------------------//
			//	    class Map_af	      //
			//------------------------------------//



/**
 * Affine radial mapping.
 * 
 * The affine radial mapping is the simplest one between the grid coordinates
 * $(\xi, \theta', \phi')$ and the physical coordinates $(r, \theta, \phi)$. 
 * It is defined by $\theta=\theta'$, $\phi=\phi'$ and 
 * \begin{itemize}
 *  \item $r=\alpha \xi + \beta$, in non-compactified domains, 
 *  \item $ u={1\over r} = \alpha \xi + \beta$ in the (usually outermost) compactified
 *        domain, 
 * \end{itemize}
 * where $\alpha$ and $\beta$ are constant in each domain. 
 * 
 * 
 * @version #$Id$#
 */

class Map_af : public Map_radial {

    // Data :
    // ----
    private:
	/// Array (size: {\tt mg->nzone}) of the values of $\alpha$ in each domain
	double* alpha ;	 
	/// Array (size: {\tt mg->nzone}) of the values of $\beta$ in each domain
	double* beta ;	  
	
    // Constructors, destructor : 
    // ------------------------
    public:
	/**
	 * Standard Constructor
	 * @param mgrille  [input] Multi-domain grid on which the mapping is defined
	 * @param r_limits [input] Array (size: number of domains + 1) of the
	 *			   value of {\it r} at the boundaries of the various 
	 *			   domains : 
	 *			   \begin{itemize}
	 *			   \item {\tt r\_limits[l]}: inner boundary of the 
	 *				 domain no. {\tt l}
	 *			   \item {\tt r\_limits[l+1]}: outer boundary of the 
	 *				 domain no. {\tt l}
	 *			   \end{itemize}
	 */
	Map_af(const Mg3d& mgrille, const double* r_limits) ;	
	
	
	Map_af(const Map_af& ) ;      /// Copy constructor
	Map_af(const Mg3d&, FILE* ) ; /// Constructor from a file (see {\tt sauve(FILE* )})
	
	/** Constructor from a general mapping.
	 * 
	 *  If the input mapping belongs to the class {\tt Map\_af}, this
	 *  constructor does the same job as the copy constructor 
	 *  {\tt Map\_af(const Map\_af\& )}.
	 * 
	 *  If the input mapping belongs to the class {\tt Map\_et}, this
	 *  constructor sets in each domain, the values of 
	 *  $\alpha$ and $\beta$ to those of the {\tt Map\_et}.
	 *  
	 */
	explicit Map_af(const Map& ) ;      

	virtual ~Map_af() ;	      /// Destructor

    // Assignment
    // ----------
    public: 
	/// Assignment to another affine mapping. 
	virtual void operator=(const Map_af& ) ;
        
    // Memory management
    // -----------------
    private:
	/// Assignment of the building functions to the member {\tt Coord}s
	void set_coord() ;	

    // Extraction of information
    // -------------------------
    public:
	/// Returns the pointer on the array {\tt alpha}
	const double* get_alpha() const ; 

	/// Returns the pointer on the array {\tt beta}
	const double* get_beta() const ; 
	
	/**
	 *  Returns the value of the radial coordinate {\it r} for a given
	 *  $(\xi, \theta', \phi')$ in a given domain. 
	 *	@param l [input] index of the domain
	 *	@param xi [input] value of $\xi$
	 *	@param theta [input] value of $\theta'$
	 *	@param pphi [input] value of $\phi'$
	 *	@return value of $r=R_l(\xi, \theta', \phi')$
	 */
	virtual double val_r(int l, double xi, double theta, double pphi) const ; 

	/**
	 * Computes the domain index {\it l} and the value of $\xi$ corresponding
	 * to a point given by its physical coordinates $(r, \theta, \phi)$. 
	 *	@param rr [input] value of {\it r}
	 *	@param theta [input] value of $\theta$
	 *	@param pphi [input] value of $\phi$
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of $\xi$
	 */
	virtual void val_lx(double rr, double theta, double pphi,
			    int& l, double& xi) const ; 
	
	/** Computes the domain index {\it l} and the value of $\xi$ corresponding
	 * to a point given by its physical coordinates $(r, \theta, \phi)$. 
	 *	@param rr [input] value of {\it r}
	 *	@param theta [input] value of $\theta$
	 *	@param pphi [input] value of $\phi$
	 *	@param par [] unused by the {\tt Map\_af} version
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of $\xi$
	 */
	virtual void val_lx(double rr, double theta, double pphi, 
			    const Param& par, int& l, double& xi) const ; 
		
	/** Returns the value of the radial coordinate {\it r} for a given
	 *  $\xi$ and a given collocation point in $(\theta', \phi')$ 
	 *   in a given domain. 
	 *	@param l [input] index of the domain
	 *	@param xi [input] value of $\xi$
	 *	@param j [input] index of the collocation point in $\theta'$
	 *	@param k [input] index of the collocation point in $\phi'$
	 *	@return value of $r=R_l(\xi, {\theta'}_j, {\phi'}_k)$
	 */
	virtual double val_r_jk(int l, double xi, int j, int k) const ; 
	
	/** Computes the domain index {\it l} and the value of $\xi$ corresponding
	 * to a point of arbitrary {\it r} but collocation values of $(\theta, \phi)$ 
	 *	@param rr [input] value of {\it r}
	 *	@param j [input] index of the collocation point in $\theta$
	 *	@param k [input] index of the collocation point in $\phi$
	 *	@param par [] unused by the {\tt Map\_af} version
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of $\xi$
	 */
	virtual void val_lx_jk(double rr, int j, int k, const Param& par, 
			       int& l, double& xi) const ; 



    // Outputs
    // -------
    public:
	virtual void sauve(FILE* ) const ;	  /// Save in a file
    
    private:
	virtual ostream& operator>>(ostream &) const ;    /// Operator >>

    // Modification of the mapping
    // ---------------------------
    public:
	/** Sets a new radial scale.
	 *	@param lambda [input] factor by which the value of {\it r} is to
	 *	    be multiplied
	 */
	virtual void homothetie(double lambda) ;
	
	/** Rescales the outer boundary of one domain.
	 *  The inner boundary is unchanged. The inner boundary 
	 *  of the next domain is changed to match the new outer
	 *  boundary. 
	 *	@param l [input] index of the domain
	 *	@param lambda [input] factor by which the value of 
	 *	    $R(\theta, \varphi)$ defining the outer boundary
	 *	    of the domain is to be multiplied. 
	 */
	virtual void resize(int l, double lambda) ; 

	/** Sets a new radial scale at the bondary between the nucleus and the
	 * first shell.
	 *	@param lambda [input] factor by which the value of {\it r} is to
	 *	    be multiplied
	 */
	void homothetie_interne(double lambda) ;
	
	/** Adaptation of the mapping to a given scalar field.
	 */
	virtual void adapt(const Cmp& ent, const Param& par) ; 

	/// Modifies the value of $\alpha$ in domain no. {\it l}
	void set_alpha(double alpha0, int l) ;

	/// Modifies the value of $\beta$ in domain no. {\it l}
	void set_beta(double beta0, int l) ;  

    // Differential operators:
    // ----------------------
    public:
	/** Computes $\partial/ \partial r$ of a {\tt Cmp}.
	 *  Note that in the external compactified domain (ZEC), it computes
	 *  $-\partial/ \partial u = r^2 \partial/ \partial r$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of {\tt ci}
	 */
	virtual void dsdr(const Cmp& ci, Cmp& resu) const ;  	    

	/** Computes $1/r \partial/ \partial \theta$ of a {\tt Cmp}.
	 *  Note that in the external compactified domain (ZEC), it computes
	 *  $1/u \partial/ \partial \theta = r \partial/ \partial \theta$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of {\tt ci}
	 */
	virtual void srdsdt(const Cmp& ci, Cmp& resu) const ;  	    
	
	/** Computes $1/(r\sin\theta) \partial/ \partial \phi$ of a {\tt Cmp}.
	 *  Note that in the external compactified domain (ZEC), it computes
	 *  $1/(u\sin\theta) \partial/ \partial \phi = 
	 *	    r/\sin\theta \partial/ \partial \phi$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of {\tt ci}
	 */
	virtual void srstdsdp(const Cmp& ci, Cmp& resu) const ;  	    
    
	/** Computes the Laplacian of a scalar field.
	 *   @param uu	[input]  Scalar field {\it u} (represented as a {\tt Cmp})
	 *			 the Laplacian $\Delta u$ of which is to be computed
	 *   @param zec_mult_r [input] Determines the quantity computed in
	 *			 the external compactified domain (ZEC) :  \\
	 *		    zec\_mult\_r = 0 : $\Delta u$	\\
	 *		    zec\_mult\_r = 2 : $r^2 \,  \Delta u$	\\
	 *		    zec\_mult\_r = 4 (default) : $r^4 \, \Delta u$	
	 *  @param lap [output] Laplacian of {\it u}
	 */
	virtual void laplacien(const Cmp& uu, int zec_mult_r, Cmp& lap) const ; 
	
	
	/** Computes the integral over all space of a {\tt Cmp}.
	 *  The computed quantity is 
	 *    $\int u \, r^2 \sin\theta \,  dr\, d\theta \, d\phi$.
	 *  The routine allocates a {\tt Tbl} (size: {\tt mg->nzone}) to store 
	 *  the result (partial integral) in each domain and returns a pointer 
	 *  to it.
	 */
	virtual Tbl* integrale(const Cmp&) const ; 
	
	 
    // PDE resolution :
    // --------------
    public:
	/** Computes the solution of a scalar Poisson equation.
	 *   @param source [input] source $\sigma$ of the Poisson equation 
	 *	    $\Delta u = \sigma$.
	 *   @param par [] not used by this {\tt Map\_af} version. 
	 *   @param uu [output] solution {\it u} with the boundary condition 
	 *	    {\it u}=0 at spatial infinity. 
	 */
	virtual void poisson(const Cmp& source, Param& par, Cmp& uu) const ;

	/** Computes the solution of a scalar Poisson equation.
	 *   The regularized source
         *   $\sigma_{\rm regu} = \sigma - \sigma_{\rm div}$
         *   is constructed and solved.
	 *   @param source [input] source $\sigma$ of the Poisson equation 
	 *	    $\Delta u = \sigma$.
	 *   @param k_div [input] regularization degree of the procedure
	 *   @param nzet [input] number of domains covering the star
	 *   @param unsgam1 [input] parameter $1/(\gamma-1)$ where $\gamma$
         *          denotes the adiabatic index.
	 *   @param par [] not used by this {\tt Map\_af} version.
	 *   @param uu [output] solution {\it u} with the boundary condition 
	 *	    {\it u}=0 at spatial infinity.
	 *   @param uu_regu [output] solution of the regular part of
	 *          the source.
         *   @param uu_div [output] solution of the diverging part of
         *          the source.
         *   @param duu_div [output] derivative of the diverging potential
         *   @param source_regu [output] regularized source
         *   @param source_div [output] diverging part of the source
	 */
	virtual void poisson_regular(const Cmp& source, int k_div, int nzet,
				     double unsgam1, Param& par, Cmp& uu,
				     Cmp& uu_regu, Cmp& uu_div,
				     Tenseur& duu_div, Cmp& source_regu,
				     Cmp& source_div) const ;

	/**
	 * Internal function intended to be used by {\tt Map::poisson\_vect}
	 * and {\tt Map::poisson\_vect\_oohara}. It constructs the sets of 
	 * parameters used for each scalar Poisson equation from the one for 
	 * the vectorial one.
	 * 
	 * In the case of a {\tt Map\_af} the result is not used and the function 
	 * only returns {\tt \& par}.
	 */
	virtual Param* donne_para_poisson_vect (Param& par, int i) const ;
	
	/**
	 * Solver of the Poisson equation with boundary condition for the 
	 * affine mapping case.
	 */
	virtual void poisson_frontiere (const Cmp&, const Valeur&, int, int, Cmp&) const ;
	virtual void poisson_frontiere_double (const Cmp& source, const Valeur& lim_func,
			const Valeur& lim_der, int num_zone, Cmp& pot) const  ;
	/**
	 * Performs the surface integration of {\tt ci} on the sphere of 
	 * radius {\tt rayon}.
	 */
	double integrale_surface (const Cmp& ci, double rayon) const ;
	
	/**
	 * Performs the surface integration of {\tt ci} at infinity.
	 * {\tt ci} must have {\tt dzpuis} =2.
	 */
	double integrale_surface_infini (const Cmp& ci) const ;
	
	
	/** Computes the solution of a 2-D Poisson equation.
	 *  The 2-D Poisson equation writes
	 *  \begin{equation}
	 *	{\partial^2 u\over\partial r^2} + 
	 *	    {1\over r} {\partial u \over \partial r} + 
	 *	    {1\over r^2} {\partial^2 u\over\partial \theta^2} = 
	 *		\sigma \ . 
	 *  \end{equation} 
	 *
	 *   @param source_mat [input] Compactly supported part of 
	 *	    the source $\sigma$ of the 2-D Poisson equation (typically
	 *	    matter terms)
	 *   @param source_quad [input] Non-compactly supported part of 
	 *	    the source $\sigma$ of the 2-D Poisson equation (typically
	 *	    quadratic terms)
	 *   @param par [output] Parameter which contains the constant
	 *			    $\lambda$ used to fulfill the virial
	 *			 identity GRV2 : \\ 
	 *  {\tt par.get\_double\_mod(0)} : [output] constant {\tt lambda}
	 *	    such that the source of the equation effectively solved
	 *	    is {\tt source\_mat + lambda * source\_quad}. 
	 *   @param uu [input/output] solution {\it u} with the boundary condition 
	 *	    {\it u}=0 at spatial infinity. 
	 */
	virtual void poisson2d(const Cmp& source_mat, const Cmp& source_quad, 
			       Param& par, Cmp& uu) const ;

	/** Performs one time-step integration of the d'Alembert scalar equation
	 *   @param par [input/output] possible parameters to control the
	 *   resolution of the d'Alembert equation: \\
	 *   {\tt par.get\_double(0)} : [input] the time step {\it dt},\\
	 *   {\tt par.get\_int(0)} : [input] the type of boundary conditions
	 *   set at the outer boundary (0 : reflexion, 1 : Sommerfeld 
	 *   outgoing wave, valid only for {\it l=0} components, 2 : Bayliss 
	 *   \& Turkel outgoing wave, valid for {\it l=0, 1, 2} components)\\
	 *   {\tt par.get\_int\_mod(0)} : [input/output] set to 0 at first
	 *   call, is used as a working flag after (must not be modified after
	 *   first call)\\
	 *   {\tt par.get\_Cmp\_mod(0)} : [input] (optional) if the wave 
	 *   equation is on a curved space-time, this is the potential in front
	 *   of the Laplace operator. It has to be updated at every time-step
	 *   (for a potential depending on time).\\
	 *   Note: there are many other working objects attached to this
	 *   {\tt Param}, so one should not modify it.
	 *   @param fJp1 [output] solution $f^{J+1}$ at time {\it J+1}
	 *   with boundary conditions defined by {\tt par.get\_int(0)}
	 *   @param fJ [input] solution $f^J$ at time {\it J}
	 *   @param fJm1 [input] solution $f^{J-1}$ at time {\it J-1}
	 *   @param source [input] source $\sigma$ of the d'Alembert equation 
	 *	    $\diamond u = \sigma$.
	 */
	virtual void dalembert(Param& par, Cmp& fJp1, const Cmp& fJ, 
			       const Cmp& fJm1, const Cmp& source) const ;

    // Building functions for the Coord's
    // ----------------------------------
    friend Mtbl* map_af_fait_r(const Map* ) ;
    friend Mtbl* map_af_fait_tet(const Map* ) ;
    friend Mtbl* map_af_fait_phi(const Map* ) ;
    friend Mtbl* map_af_fait_sint(const Map* ) ;
    friend Mtbl* map_af_fait_cost(const Map* ) ;
    friend Mtbl* map_af_fait_sinp(const Map* ) ;
    friend Mtbl* map_af_fait_cosp(const Map* ) ;

    friend Mtbl* map_af_fait_x(const Map* ) ;
    friend Mtbl* map_af_fait_y(const Map* ) ;
    friend Mtbl* map_af_fait_z(const Map* ) ;

    friend Mtbl* map_af_fait_xa(const Map* ) ;
    friend Mtbl* map_af_fait_ya(const Map* ) ;
    friend Mtbl* map_af_fait_za(const Map* ) ;

    friend Mtbl* map_af_fait_xsr(const Map* ) ;
    friend Mtbl* map_af_fait_dxdr(const Map* ) ;
    friend Mtbl* map_af_fait_srdrdt(const Map* ) ;
    friend Mtbl* map_af_fait_srstdrdp(const Map* ) ;
    friend Mtbl* map_af_fait_sr2drdt(const Map* ) ;
    friend Mtbl* map_af_fait_sr2stdrdp(const Map* ) ;
    friend Mtbl* map_af_fait_d2rdx2(const Map* ) ;
    friend Mtbl* map_af_fait_lapr_tp(const Map* ) ;
    friend Mtbl* map_af_fait_d2rdtdx(const Map* ) ;
    friend Mtbl* map_af_fait_sstd2rdpdx(const Map* ) ;
    friend Mtbl* map_af_fait_sr2d2rdt2(const Map* ) ;

};

     Mtbl* map_af_fait_r(const Map* ) ;
     Mtbl* map_af_fait_tet(const Map* ) ;
     Mtbl* map_af_fait_phi(const Map* ) ;
     Mtbl* map_af_fait_sint(const Map* ) ;
     Mtbl* map_af_fait_cost(const Map* ) ;
     Mtbl* map_af_fait_sinp(const Map* ) ;
     Mtbl* map_af_fait_cosp(const Map* ) ;

     Mtbl* map_af_fait_x(const Map* ) ;
     Mtbl* map_af_fait_y(const Map* ) ;
     Mtbl* map_af_fait_z(const Map* ) ;

     Mtbl* map_af_fait_xa(const Map* ) ;
     Mtbl* map_af_fait_ya(const Map* ) ;
     Mtbl* map_af_fait_za(const Map* ) ;

     Mtbl* map_af_fait_xsr(const Map* ) ;
     Mtbl* map_af_fait_dxdr(const Map* ) ;
     Mtbl* map_af_fait_srdrdt(const Map* ) ;
     Mtbl* map_af_fait_srstdrdp(const Map* ) ;
     Mtbl* map_af_fait_sr2drdt(const Map* ) ;
     Mtbl* map_af_fait_sr2stdrdp(const Map* ) ;
     Mtbl* map_af_fait_d2rdx2(const Map* ) ;
     Mtbl* map_af_fait_lapr_tp(const Map* ) ;
     Mtbl* map_af_fait_d2rdtdx(const Map* ) ;
     Mtbl* map_af_fait_sstd2rdpdx(const Map* ) ;
     Mtbl* map_af_fait_sr2d2rdt2(const Map* ) ;




                        //------------------------------------//
			//	    class Map_et	      //
			//------------------------------------//



/**
 * Radial mapping of rather general form.
 *
 * This mapping relates the grid coordinates
 * $(\xi, \theta', \phi')$ and the physical coordinates $(r, \theta, \phi)$
 * in the following manner [see Bonazzola, Gourgoulhon \& Marck,
 * {\sl Phys. Rev. D} {\bf 58}, 104020 (1998) for details]:
 * $\theta=\theta'$, $\phi=\phi'$ and
 * \begin{itemize}
 *  \item $r = \alpha [\xi + A(\xi) F(\theta', \phi') + B(\xi) G(\theta', \phi')]
 *	       + \beta$ in non-compactified domains,
 *  \item $ u={1\over r} = \alpha [\xi + A(\xi) F(\theta', \phi')] + \beta$ in
 *	    the (usually outermost) compactified domain,
 * \end{itemize}
 * where $\alpha$ and $\beta$ are constant in each domain, $A(\xi)$ and
 * $B(\xi)$ are constant polynomials defined by
 * \begin{itemize}
 *	\item $A(\xi) = 3 \xi^4 - 2 \xi^6$ and 
 *	      $B(\xi) = (5\xi^3 - 3\xi^5)/2$ in the nucleus (innermost domain
 *		    which contains {\it r}=0);
 *	\item $A(\xi) = (\xi^3 - 3\xi + 2)/4$ and
 *	      $B(\xi) = (-\xi^3 + 3\xi +2)/4$ in the other domains.
 *  \end{itemize}
 * The functions $F(\theta', \phi')$ and $G(\theta', \phi')$ depend on the
 * domain under consideration and define the boundaries of this domain. 	
 *    
 * 
 * @version #$Id$#
 */

class Map_et : public Map_radial {

    // Data :
    // ----
    private:
	/// Array (size: {\tt mg->nzone}) of the values of $\alpha$ in each domain
	double* alpha ;	  
	/// Array (size: {\tt mg->nzone}) of the values of $\beta$ in each domain
	double* beta ;	  

	/** Array (size: {\tt mg->nzone}) of {\tt Tbl} which stores the 
	 *  values of $A(\xi)$ in each domain
	 */
	Tbl** aa ;	  

	/** Array (size: {\tt mg->nzone}) of {\tt Tbl} which stores the 
	 *  values of $A'(\xi)$ in each domain
	 */
	Tbl** daa ;	  
	
	/** Array (size: {\tt mg->nzone}) of {\tt Tbl} which stores the 
	 *  values of $A''(\xi)$ in each domain
	 */
	Tbl** ddaa ;	  
	
	/// Values at the {\tt nr} collocation points of $A(\xi)/\xi$ in the nucleus
	Tbl aasx ; 
	 
	/// Values at the {\tt nr} collocation points of $A(\xi)/\xi^2$ in the nucleus
	Tbl aasx2 ; 
	 
	/** Values at the {\tt nr} collocation points of $A(\xi)/(\xi-1)$
	 *  in the outermost compactified domain
	 */
	Tbl zaasx ; 
	  
	/** Values at the {\tt nr} collocation points of $A(\xi)/(\xi-1)^2$
	 *  in the outermost compactified domain
	 */
	Tbl zaasx2 ; 
	  
	/** Array (size: {\tt mg->nzone}) of {\tt Tbl} which stores the 
	 *  values of $B(\xi)$ in each domain
	 */
	Tbl** bb ;	  

	/** Array (size: {\tt mg->nzone}) of {\tt Tbl} which stores the 
	 *  values of $B'(\xi)$ in each domain
	 */
	Tbl** dbb ;	  
	
	/** Array (size: {\tt mg->nzone}) of {\tt Tbl} which stores the 
	 *  values of $B''(\xi)$ in each domain
	 */
	Tbl** ddbb ;	  
	
	/// Values at the {\tt nr} collocation points of $B(\xi)/\xi$ in the nucleus
	Tbl bbsx ; 
	 
	/// Values at the {\tt nr} collocation points of $B(\xi)/\xi^2$ in the nucleus
	Tbl bbsx2 ; 
	 
	/** Values of the function $F(\theta', \phi')$ at the {\tt nt*np}
	 * angular collocation points in each domain.
	 * The {\tt Valeur ff} is defined on the multi-grid {\tt mg->g\_angu}
	 * (cf. class {\tt Mg3d}). 
	 */
	Valeur ff ; 
	  
	/** Values of the function $G(\theta', \phi')$ at the {\tt nt*np}
	 *   angular collocation points in each domain.
	 * The {\tt Valeur gg} is defined on the multi-grid {\tt mg->g\_angu}
	 * (cf. class {\tt Mg3d}). 
	 */
	Valeur gg ; 
	  
    public:
	/** $1/(\partial R/\partial \xi) R/\xi$ in the nucleus; \\
	 *  $1/(\partial R/\partial \xi) R/(\xi + \beta/\alpha)$ in the shells; \\
	 *  $1/(\partial U/\partial \xi) U/(\xi-1)$ in the outermost 
	 *  compactified domain.
	 */
	Coord rsxdxdr ;
	 
	/** $[ R/ (\alpha \xi + \beta) ]^2 (\partial R/\partial \xi) / \alpha$
	 *  in the nucleus and the shells; \\
	 *  $\partial U/\partial \xi / \alpha$ in the outermost compactified
	 *  domain. 
	 */
	Coord rsx2drdx ;
	  
    // Constructors, destructor : 
    // ------------------------
    public:
	/**
	 * Standard Constructor
	 * @param mgrille  [input] Multi-domain grid on which the mapping is defined
	 * @param r_limits [input] Array (size: number of domains + 1) of the
	 *			   value of {\it r} at the boundaries of the various 
	 *			   domains : 
	 *			   \begin{itemize}
	 *			   \item {\tt r\_limits[l]}: inner boundary of the 
	 *				 domain no. {\tt l}
	 *			   \item {\tt r\_limits[l+1]}: outer boundary of the 
	 *				 domain no. {\tt l}
	 *			   \end{itemize}
	 */
	Map_et(const Mg3d& mgrille, const double* r_limits) ;	
	
	/**
	 * Constructor using the equation of the surface of the star.
	 * @param mgrille [input] Multi-domain grid on which the mapping is defined
	 * It must contains at least one shell.
	 * @param r_limits [input] Array (size: number of domains + 1) of the
	 *			   value of {\it r} at the boundaries of the various 
	 *			   domains : 
	 *			   \begin{itemize}
	 *			   \item {\tt r\_limits[l]}: inner boundary of the 
	 *				 domain no. {\tt l}
	 *			   \item {\tt r\_limits[l+1]}: outer boundary of the 
	 *				 domain no. {\tt l}
	 *			   \end{itemize}
	 * The first value is not used.
	 * @param tab [input] equation of the surface between the nucleus and the first
	 * shell in the form : ${\rm tab}(k,j) = r(\phi_k,\theta_j)$, where
	 * $\phi_k$ and $\theta_j$ are the values of the angular colocation points.
	 */

	Map_et(const Mg3d& mgrille, const double* r_limits,const Tbl& tab);	
	Map_et(const Map_et& ) ;      /// Copy constructor
	Map_et(const Mg3d&, FILE* ) ; /// Constructor from a file (see {\tt sauve(FILE* )})

	virtual ~Map_et() ;	      /// Destructor

    // Assignment
    // ----------
    public:
	/// Assignment to another {\tt Map\_et} 
	virtual void operator=(const Map_et& ) ;
	
	/// Assignment to an affine mapping. 
	virtual void operator=(const Map_af& ) ;
	
	/// Assigns a given value to the function $F(\theta',\phi')$ 
	void set_ff(const Valeur& ) ; 
	/// Assigns a given value to the function $G(\theta',\phi')$ 
	void set_gg(const Valeur& ) ; 

    // Memory management
    // -----------------
    private:
	/// Assignement of the building functions to the member {\tt Coord}s
	void set_coord() ;		
    protected:
	/// Resets all the member {\tt Coord}s	
	virtual void reset_coord() ;
	
    private: 
	/// Construction of the polynomials $A(\xi)$ and $B(\xi)$
	void fait_poly() ; 	
	
    // Extraction of information
    // -------------------------
    public:
	/** Returns a pointer on the array {\tt alpha} (values of $\alpha$
	 *   in each domain)
	 */
	const double* get_alpha() const ; 

	/** Returns a pointer on the array {\tt beta} (values of $\beta$
	 *   in each domain)
	 */
	const double* get_beta() const ; 

	/// Returns a (constant) reference to the function $F(\theta',\phi')$
	const Valeur& get_ff() const ; 

	/// Returns a (constant) reference to the function $G(\theta',\phi')$
	const Valeur& get_gg() const ; 
	
	/**
	 *  Returns the value of the radial coordinate {\it r} for a given
	 *  $(\xi, \theta', \phi')$ in a given domain. 
	 *	@param l [input] index of the domain
	 *	@param xi [input] value of $\xi$
	 *	@param theta [input] value of $\theta'$
	 *	@param pphi [input] value of $\phi'$
	 *	@return value of $r=R_l(\xi, \theta', \phi')$
	 */
	virtual double val_r(int l, double xi, double theta, double pphi) const ; 

	/**
	 * Computes the domain index {\it l} and the value of $\xi$ corresponding
	 * to a point given by its physical coordinates $(r, \theta, \phi)$. 
	 *	@param rr [input] value of {\it r}
	 *	@param theta [input] value of $\theta$
	 *	@param pphi [input] value of $\phi$
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of $\xi$
	 */
	virtual void val_lx(double rr, double theta, double pphi,
			    int& l, double& xi) const ; 
	
	/** Computes the domain index {\it l} and the value of $\xi$ corresponding
	 * to a point given by its physical coordinates $(r, \theta, \phi)$. 
	 * This version enables to pass some parameters to control the
	 * accuracy of the computation. 
	 *	@param rr [input] value of {\it r}
	 *	@param theta [input] value of $\theta$
	 *	@param pphi [input] value of $\phi$
	 *	@param par [input/output] parameters to control the
	 *	    accuracy of the computation:  \\
	 *  {\tt par.get\_int(0)} : [input] maximum number of iterations in the 
	 *		    secant method to locate $\xi$ \\
	 *  {\tt par.get\_int\_mod(0)} : [output] effective number of iterations 
	 *				    used \\
	 *  {\tt par.get\_double(0)} : [input] absolute precision in the secant 
	 *				method to locate $\xi$ 
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of $\xi$
	 */
	virtual void val_lx(double rr, double theta, double pphi, 
			    const Param& par, int& l, double& xi) const ; 
		
	/** Returns the value of the radial coordinate {\it r} for a given
	 *  $\xi$ and a given collocation point in $(\theta', \phi')$ 
	 *   in a given domain. 
	 *	@param l [input] index of the domain
	 *	@param xi [input] value of $\xi$
	 *	@param j [input] index of the collocation point in $\theta'$
	 *	@param k [input] index of the collocation point in $\phi'$
	 *	@return value of $r=R_l(\xi, {\theta'}_j, {\phi'}_k)$
	 */
	virtual double val_r_jk(int l, double xi, int j, int k) const ; 
	
	/** Computes the domain index {\it l} and the value of $\xi$ corresponding
	 * to a point of arbitrary {\it r} but collocation values of $(\theta, \phi)$ 
	 *	@param rr [input] value of {\it r}
	 *	@param j [input] index of the collocation point in $\theta$
	 *	@param k [input] index of the collocation point in $\phi$
	 *	@param par [input/output] parameters to control the
	 *	    accuracy of the computation:  \\
	 *  {\tt par.get\_int(0)} : [input] maximum number of iterations in the 
	 *		    secant method to locate $\xi$ \\
	 *  {\tt par.get\_int\_mod(0)} : [output] effective number of iterations 
	 *				    used \\
	 *  {\tt par.get\_double(0)} : [input] absolute precision in the secant 
	 *				method to locate $\xi$ 
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of $\xi$
	 */
	virtual void val_lx_jk(double rr, int j, int k, const Param& par, 
			       int& l, double& xi) const ; 



    // Outputs
    // -------
    public:
	virtual void sauve(FILE* ) const ;	  /// Save in a file
    
    private:
	virtual ostream& operator>>(ostream &) const ;    /// Operator >>

    // Modification of the radial scale
    // --------------------------------
    public:
	/** Sets a new radial scale.
	 *	@param lambda [input] factor by which the value of {\it r} is to
	 *	    be multiplied
	 */
	virtual void homothetie(double lambda) ;	

	/** Rescales the outer boundary of one domain.
	 *  The inner boundary is unchanged. The inner boundary 
	 *  of the next domain is changed to match the new outer
	 *  boundary. 
	 *	@param l [input] index of the domain
	 *	@param lambda [input] factor by which the value of 
	 *	    $R(\theta, \varphi)$ defining the outer boundary
	 *	    of the domain is to be multiplied. 
	 */
	virtual void resize(int l, double lambda) ; 

    // Modification of the mapping
    // ---------------------------
	/** Adaptation of the mapping to a given scalar field.
	 *  Computes the functions $F(\theta',\phi')$ and $G(\theta',\phi')$
	 *  as well as the factors $\alpha$ and $\beta$, so that the 
	 *  boundaries of some domains coincide with isosurfaces of the
	 *  scalar field {\tt ent}.
	 *  @param ent [input] scalar field, the isosurfaces of which are
	 *     used to determine the mapping
	 *  @param par [input/output] parameters of the computation: \\
	 *   {\tt par.get\_int(0)} : maximum number of iterations to locate
	 *     zeros by the secant method \\
	 *   {\tt par.get\_int(1)} : number of domains where the adjustment
	 *     to the isosurfaces of {\tt ent} is to be performed \\
	 *   {\tt par.get\_int(2)} : number of domains {\tt nz\_search} where 
	 *      the isosurfaces will be searched : the routine scans the 
	 *	{\tt nz\_search} innermost domains, starting from the domain
	 *	of index {\tt nz\_search-1}. NB: the field {\tt ent} must
	 *	be continuous over these domains \\
	 *   {\tt par.get\_int(3)} : 1 = performs the full computation, 
	 *			     0 = performs only the rescaling 
	 *				by the factor 
	 *				{\tt par.get\_double\_mod(0)} \\
	 *   {\tt par.get\_int(4)} : theta index of the collocation point
	 *			     $(\theta_*, \phi_*)$ [using the notations
	 *    of Bonazzola, Gourgoulhon \& Marck, {\sl Phys. Rev. D} {\bf 58}, 
	 *    104020 (1998)] defining an isosurface of {\tt ent} \\
	 *   {\tt par.get\_int(5)} : phi index of the collocation point
	 *			     $(\theta_*, \phi_*)$ [using the notations
	 *    of Bonazzola, Gourgoulhon \& Marck, {\sl Phys. Rev. D} {\bf 58}, 
	 *    104020 (1998)] defining an isosurface of {\tt ent} \\
	 *   {\tt par.get\_int\_mod(0)} [output] : number of iterations 
	 *	actually used in the secant method \\
	 *   {\tt par.get\_double(0)} : required absolute precision in the
	 *	 determination of zeros by the secant method \\
	 *   {\tt par.get\_double(1)} : factor by which the values of $\lambda$
	 *	    and $\mu$ [using the notations
	 *    of Bonazzola, Gourgoulhon \& Marck, {\sl Phys. Rev. D} {\bf 58}, 
	 *    104020 (1998)] will be multiplied : 1 = regular mapping, 
	 *	    0 = contracting mapping \\
	 *   {\tt par.get\_double(2)} : factor by which all the radial distances
	 *				will be multiplied \\
	 *   {\tt par.get\_tbl(0)} : array of values of the field {\tt ent} to
	 *			     define the isosurfaces. 
	 *  
	 */
	virtual void adapt(const Cmp& ent, const Param& par)  ; 

    // Differential operators:
    // ----------------------
    public:
	/** Computes $\partial/ \partial r$ of a {\tt Cmp}.
	 *  Note that in the external compactified domain (ZEC), it computes
	 *  $-\partial/ \partial u = r^2 \partial/ \partial r$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of {\tt ci}
	 */
	virtual void dsdr(const Cmp& ci, Cmp& resu) const ;  	    

	/** Computes $1/r \partial/ \partial \theta$ of a {\tt Cmp}.
	 *  Note that in the external compactified domain (ZEC), it computes
	 *  $1/u \partial/ \partial \theta = r \partial/ \partial \theta$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of {\tt ci}
	 */
	virtual void srdsdt(const Cmp& ci, Cmp& resu) const ;  	    
	
	/** Computes $1/(r\sin\theta) \partial/ \partial \phi$ of a {\tt Cmp}.
	 *  Note that in the external compactified domain (ZEC), it computes
	 *  $1/(u\sin\theta) \partial/ \partial \phi = 
	 *	    r/\sin\theta \partial/ \partial \phi$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of {\tt ci}
	 */
	virtual void srstdsdp(const Cmp& ci, Cmp& resu) const ;  	    
    
	/** Computes the Laplacian of a scalar field.
	 *   @param uu	[input]  Scalar field {\it u} (represented as a {\tt Cmp})
	 *			 the Laplacian $\Delta u$ of which is to be computed
	 *   @param zec_mult_r [input] Determines the quantity computed in
	 *			 the external compactified domain (ZEC) :  \\
	 *		    zec\_mult\_r = 0 : $\Delta u$	\\
	 *		    zec\_mult\_r = 2 : $r^2 \,  \Delta u$	\\
	 *		    zec\_mult\_r = 4 (default) : $r^4 \, \Delta u$	
	 *  @param lap [output] Laplacian of {\it u}
	 */
	virtual void laplacien(const Cmp& uu, int zec_mult_r, Cmp& lap) const ; 
	

	/** Computes the integral over all space of a {\tt Cmp}.
	 *  The computed quantity is 
	 *    $\int u \, r^2 \sin\theta \,  dr\, d\theta \, d\phi$.
	 *  The routine allocates a {\tt Tbl} (size: {\tt mg->nzone}) to store 
	 *  the result (partial integral) in each domain and returns a pointer 
	 *  to it.
	 */
	virtual Tbl* integrale(const Cmp&) const ; 
	
	 
    // PDE resolution :
    // --------------
    public:
	/** Computes the solution of a scalar Poisson equation.
	 * 
	 * Following the method explained in Sect. III.C of Bonazzola, 
	 * Gourgoulhon \& Marck, {\sl Phys. Rev. D} {\bf 58}, 104020 (1998),  
	 * the Poisson equation $\Delta u = \sigma$ is re-written
	 * as $a \tilde\Delta u = \sigma + R(u)$,  where $\tilde\Delta$
	 * is the Laplacian in an affine mapping and {\it R(u)} contains the
	 * terms generated by the deviation of the mapping {\tt *this}
	 * from spherical symmetry. This equation is solved by iterations.
	 * At each step {\it J} the equation effectively solved is 
	 *  $\tilde\Delta u^{J+1} = s^J$ where
	 * \begin{equation}
	 *   s^J = 1/a_l^{\rm max} \{ {\tt source} + R(u^J) + (a_l^{\rm max}-a)
	 *          [ \lambda s^{J-1} + (1-\lambda) s^{J-2} ] \} \ ,  
	 * \end{equation}
	 * with $a_l^{\rm max} := \max(a)$ in domain no. {\it l} and $\lambda$
	 * is a relaxation parameter. 
	 *  @param source [input] source $\sigma$ of the Poisson equation
	 *  @param par [input/output] parameters for the iterative method: \\
	 *  {\tt par.get\_int(0)} : [input] maximum number of iterations \\
	 *  {\tt par.get\_double(0)} : [input] relaxation parameter $\lambda$ \\
	 *  {\tt par.get\_double(1)} : [input] required precision: the iterative
	 *	  method is stopped as soon as the relative difference between
	 *	 $u^J$ and $u^{J-1}$ is greater than {\tt par.get\_double(1)}\\
	 *  {\tt par.get\_cmp\_mod(0)} : [input/output] input : {\tt Cmp} 
	 *		containing $s^{J-1}$ (cf. the above equation) to 
	 *		start the iteration (if nothing is known a priori, 
	 *		this {\tt Cmp} must be set to zero); output: value
	 *		of $s^{J-1}$, to used in a next call to the routine \\
	 *  {\tt par.get\_int\_mod(0)} : [output] number of iterations 
	 *				    actually used to get the solution.
	 *
	 *  @param uu [input/output] input : previously computed value of {\it u}
	 *	to start the iteration (term {\it R(u)}) (if nothing is known a 
	 *	priori, {\tt uu} must be set to zero); output: solution {\it u} 
	 *	with the boundary condition {\it u}=0 at spatial infinity. 
	 */
	virtual void poisson(const Cmp& source, Param& par, Cmp& uu) const ;

	/** Computes the solution of a scalar Poisson equation.
	 *   The regularized source
	 *   @param source [input] source $\sigma$ of the Poisson equation 
	 *	    $\Delta u = \sigma$.
	 *   @param k_div [input] regularization degree of the procedure
	 *   @param nzet [input] number of domains covering the star
	 *   @param unsgam1 [input] parameter $1/(\gamma-1)$ where $\gamma$
         *          denotes the adiabatic index.
	 *   @param par [input/output] parameters for the iterative method: \\
	 *  {\tt par.get\_int(0)} : [input] maximum number of iterations \\
	 *  {\tt par.get\_double(0)} : [input] relaxation parameter
	 *                             $\lambda$ \\
	 *  {\tt par.get\_double(1)} : [input] required precision: the
	 *        iterative method is stopped as soon as the relative
	 *        difference between $u^J$ and $u^{J-1}$ is greater than
	 *        {\tt par.get\_double(1)}\\
	 *  {\tt par.get\_cmp\_mod(0)} : [input/output] input : {\tt Cmp} 
	 *		containing $s^{J-1}$ (cf. the above equation) to 
	 *		start the iteration (if nothing is known a priori, 
	 *		this {\tt Cmp} must be set to zero); output: value
	 *		of $s^{J-1}$, to used in a next call to the routine \\
	 *  {\tt par.get\_int\_mod(0)} : [output] number of iterations 
	 *				    actually used to get the solution.
	 *   @param uu [input/output] input : previously computed value of {\it u}
	 *	to start the iteration (term {\it R(u)}) (if nothing is known a 
	 *	priori, {\tt uu} must be set to zero); output: solution {\it u} 
	 *	with the boundary condition {\it u}=0 at spatial infinity.
	 *   @param uu_regu [output] solution of the regular part of
	 *          the  source.
         *   @param uu_div [output] solution of the diverging part of
         *          the source.
         *   @param duu_div [output] derivative of the diverging potential
         *   @param source_regu [output] regularized source
         *   @param source_div [output] diverging part of the source
	 */
	virtual void poisson_regular(const Cmp& source, int k_div, int nzet,
				     double unsgam1, Param& par, Cmp& uu,
				     Cmp& uu_regu, Cmp& uu_div,
				     Tenseur& duu_div, Cmp& source_regu,
				     Cmp& source_div) const ;

	/**
	 * Internal function intended to be used by {\tt Map::poisson\_vect}
	 * and {\tt Map::poisson\_vect\_oohara}. It constructs the sets of 
	 * parameters used for each scalar Poisson equation from the one for 
	 * the vectorial one.
	 * 
	 * @param para [input] : the {\tt Param} used for the resolution of 
	 * the vectorial Poisson equation : \\
	 * {\tt para.int()} maximum number of iteration.\\
	 * {\tt para.double(0)} relaxation parameter.\\
	 * {\tt para.double(1)} required precision. \\
	 * {\tt para.tenseur\_mod()} source of the vectorial part at the previous 
	 * step.\\
	 * {\tt para.cmp\_mod()} source of the scalar part at the previous 
	 * step.
	 * 
	 * @param i [input] number of the scalar Poisson equation that is being 
	 * solved (values from 0 to 2 for the componants of the vectorial part
	 * and 3 for the scalar one).
	 * 
	 * @return the pointer on the parameter set used for solving the scalar 
	 * Poisson equation labelled by {\it i}.
	 */
	virtual Param* donne_para_poisson_vect (Param& para, int i) const ;
	
	/**
	 * Not yet implemented.
	 */
	virtual void poisson_frontiere (const Cmp&, const Valeur&, int, int, 
					Cmp&) const ;
	virtual void poisson_frontiere_double (const Cmp& source, 
			const Valeur& lim_func, const Valeur& lim_der, 
			int num_zone, Cmp& pot) const  ;

	/** Computes the solution of a 2-D Poisson equation.
	 *  The 2-D Poisson equation writes
	 *  \begin{equation}
	 *	{\partial^2 u\over\partial r^2} + 
	 *	    {1\over r} {\partial u \over \partial r} + 
	 *	    {1\over r^2} {\partial^2 u\over\partial \theta^2} = 
	 *		\sigma \ . 
	 *  \end{equation} 
	 *
	 *   @param source_mat [input] Compactly supported part of 
	 *	    the source $\sigma$ of the 2-D Poisson equation (typically
	 *	    matter terms)
	 *   @param source_quad [input] Non-compactly supported part of 
	 *	    the source $\sigma$ of the 2-D Poisson equation (typically
	 *	    quadratic terms)
	 *   @param par [input/output] Parameters to control the resolution : \\ 
	 *  {\tt par.get\_double\_mod(0)} : [output] constant {\tt lambda}
	 *	    such that the source of the equation effectively solved
	 *	    is {\tt source\_mat + lambda * source\_quad}, in order to
	 *	    fulfill the virial identity GRV2.  \\
	 *  If the theta basis is {\tt T\_SIN\_I}, the following arguments
	 *  are required: \\
	 *  {\tt par.get\_int(0)} : [input] maximum number of iterations \\
	 *  {\tt par.get\_double(0)} : [input] relaxation parameter \\
	 *  {\tt par.get\_double(1)} : [input] required precision: the iterative
	 *	  method is stopped as soon as the relative difference between
	 *	 $u^J$ and $u^{J-1}$ is greater than {\tt par.get\_double(1)}\\
	 *  {\tt par.get\_cmp\_mod(0)} : [input/output] input : {\tt Cmp} 
	 *		containing $s^{J-1}$ to 
	 *		start the iteration (if nothing is known a priori, 
	 *		this {\tt Cmp} must be set to zero); output: value
	 *		of $s^{J-1}$, to used in a next call to the routine \\
	 *  {\tt par.get\_int\_mod(0)} : [output] number of iterations 
	 *				    actually used to get the solution.
	 *	 
	 *   @param uu [input/output] solution {\it u} with the boundary condition 
	 *	    {\it u}=0 at spatial infinity. 
	 */
	virtual void poisson2d(const Cmp& source_mat, const Cmp& source_quad, 
			       Param& par, Cmp& uu) const ;

	/**
	 * Not yet implemented.
	 */
	virtual void dalembert(Param& par, Cmp& fJp1, const Cmp& fJ, 
			       const Cmp& fJm1, const Cmp& source) const ;




    // Building functions for the Coord's
    // ----------------------------------
    friend Mtbl* map_et_fait_r(const Map* ) ;
    friend Mtbl* map_et_fait_tet(const Map* ) ;
    friend Mtbl* map_et_fait_phi(const Map* ) ;
    friend Mtbl* map_et_fait_sint(const Map* ) ;
    friend Mtbl* map_et_fait_cost(const Map* ) ;
    friend Mtbl* map_et_fait_sinp(const Map* ) ;
    friend Mtbl* map_et_fait_cosp(const Map* ) ;

    friend Mtbl* map_et_fait_x(const Map* ) ;
    friend Mtbl* map_et_fait_y(const Map* ) ;
    friend Mtbl* map_et_fait_z(const Map* ) ;

    friend Mtbl* map_et_fait_xa(const Map* ) ;
    friend Mtbl* map_et_fait_ya(const Map* ) ;
    friend Mtbl* map_et_fait_za(const Map* ) ;

    friend Mtbl* map_et_fait_xsr(const Map* ) ;
    friend Mtbl* map_et_fait_dxdr(const Map* ) ;
    friend Mtbl* map_et_fait_srdrdt(const Map* ) ;
    friend Mtbl* map_et_fait_srstdrdp(const Map* ) ;
    friend Mtbl* map_et_fait_sr2drdt(const Map* ) ;
    friend Mtbl* map_et_fait_sr2stdrdp(const Map* ) ;
    friend Mtbl* map_et_fait_d2rdx2(const Map* ) ;
    friend Mtbl* map_et_fait_lapr_tp(const Map* ) ;
    friend Mtbl* map_et_fait_d2rdtdx(const Map* ) ;
    friend Mtbl* map_et_fait_sstd2rdpdx(const Map* ) ;
    friend Mtbl* map_et_fait_sr2d2rdt2(const Map* ) ;

    friend Mtbl* map_et_fait_rsxdxdr(const Map* ) ;
    friend Mtbl* map_et_fait_rsx2drdx(const Map* ) ;

};

     Mtbl* map_et_fait_r(const Map* ) ;
     Mtbl* map_et_fait_tet(const Map* ) ;
     Mtbl* map_et_fait_phi(const Map* ) ;
     Mtbl* map_et_fait_sint(const Map* ) ;
     Mtbl* map_et_fait_cost(const Map* ) ;
     Mtbl* map_et_fait_sinp(const Map* ) ;
     Mtbl* map_et_fait_cosp(const Map* ) ;

     Mtbl* map_et_fait_x(const Map* ) ;
     Mtbl* map_et_fait_y(const Map* ) ;
     Mtbl* map_et_fait_z(const Map* ) ;

     Mtbl* map_et_fait_xa(const Map* ) ;
     Mtbl* map_et_fait_ya(const Map* ) ;
     Mtbl* map_et_fait_za(const Map* ) ;

     Mtbl* map_et_fait_xsr(const Map* ) ;
     Mtbl* map_et_fait_dxdr(const Map* ) ;
     Mtbl* map_et_fait_srdrdt(const Map* ) ;
     Mtbl* map_et_fait_srstdrdp(const Map* ) ;
     Mtbl* map_et_fait_sr2drdt(const Map* ) ;
     Mtbl* map_et_fait_sr2stdrdp(const Map* ) ;
     Mtbl* map_et_fait_d2rdx2(const Map* ) ;
     Mtbl* map_et_fait_lapr_tp(const Map* ) ;
     Mtbl* map_et_fait_d2rdtdx(const Map* ) ;
     Mtbl* map_et_fait_sstd2rdpdx(const Map* ) ;
     Mtbl* map_et_fait_sr2d2rdt2(const Map* ) ;

     Mtbl* map_et_fait_rsxdxdr(const Map* ) ;
     Mtbl* map_et_fait_rsx2drdx(const Map* ) ;



#endif
