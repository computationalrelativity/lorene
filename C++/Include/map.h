/*
 *  Definition of Lorene classes Map
 *				 Map_radial
 *				 Map_af
 *				 Map_et
 *                               Map_log
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2003 Eric Gourgoulhon
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
 * Revision 1.37  2004/12/02 09:33:04  p_grandclement
 * *** empty log message ***
 *
 * Revision 1.36  2004/11/30 20:42:05  k_taniguchi
 * Addition of some functions with the falloff condition and a method
 * to resize the external shell.
 *
 * Revision 1.35  2004/11/23 12:39:12  f_limousin
 * Intoduce function poisson_dir_neu(...) to solve a scalar poisson
 * equation with a mixed boundary condition (Dirichlet + Neumann).
 *
 * Revision 1.34  2004/10/11 15:08:59  j_novak
 * The radial manipulation functions take Scalar as arguments, instead of Cmp.
 * Added a conversion operator from Scalar to Cmp.
 * The Cmp radial manipulation function make conversion to Scalar, call to the
 * Map_radial version with a Scalar argument and back.
 *
 * Revision 1.33  2004/10/08 13:34:35  j_novak
 * Scalar::div_r() does not need to pass through Cmp version anymore.
 *
 * Revision 1.32  2004/08/24 09:14:40  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.31  2004/07/27 08:24:26  j_novak
 * Modif. comments
 *
 * Revision 1.30  2004/07/26 16:02:21  j_novak
 * Added a flag to specify whether the primitive should be zero either at r=0
 * or at r going to infinity.
 *
 * Revision 1.29  2004/06/22 08:49:57  p_grandclement
 * Addition of everything needed for using the logarithmic mapping
 *
 * Revision 1.28  2004/06/14 15:23:53  e_gourgoulhon
 * Added virtual function primr for computation of radial primitives.
 *
 * Revision 1.27  2004/03/31 11:22:23  f_limousin
 * Methods Map_et::poisson_interne and Map_af::poisson_interne have been
 * implemented to solve the continuity equation for strange stars.
 *
 * Revision 1.26  2004/03/22 13:12:41  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.24  2004/03/01 09:57:02  j_novak
 * the wave equation is solved with Scalars. It now accepts a grid with a
 * compactified external domain, which the solver ignores and where it copies
 * the values of the field from one time-step to the next.
 *
 * Revision 1.23  2004/02/11 09:47:44  p_grandclement
 * Addition of a new elliptic solver, matching with the homogeneous solution
 * at the outer shell and not solving in the external domain (more details
 * coming soon ; check your local Lorene dealer...)
 *
 * Revision 1.22  2004/01/29 08:50:01  p_grandclement
 * Modification of Map::operator==(const Map&) and addition of the surface
 * integrales using Scalar.
 *
 * Revision 1.21  2004/01/28 16:46:22  p_grandclement
 * Addition of the sol_elliptic_fixe_der_zero stuff
 *
 * Revision 1.20  2004/01/28 10:35:52  j_novak
 * Added new methods mult_r() for Scalars. These do not change the dzpuis flag.
 *
 * Revision 1.19  2004/01/27 09:33:46  j_novak
 * New method Map_radial::div_r_zec
 *
 * Revision 1.18  2004/01/26 16:16:15  j_novak
 * Methods of gradient for Scalar s. The input can have any dzpuis.
 *
 * Revision 1.17  2004/01/19 21:38:21  e_gourgoulhon
 * Corrected sign error in comments of Map_radial::dxdr.
 *
 * Revision 1.16  2003/12/30 22:52:47  e_gourgoulhon
 * Class Map: added methods flat_met_spher() and flat_met_cart() to get
 * flat metric associated with the coordinates described by the mapping.
 *
 * Revision 1.15  2003/12/11 14:48:47  p_grandclement
 * Addition of ALL (and that is a lot !) the files needed for the general elliptic solver ... UNDER DEVELOPEMENT...
 *
 * Revision 1.14  2003/11/06 14:43:37  e_gourgoulhon
 * Gave a name to const arguments in certain method prototypes (e.g.
 * constructors) to correct a bug of DOC++.
 *
 * Revision 1.13  2003/11/04 22:54:49  e_gourgoulhon
 * Added new virtual methods mult_cost, mult_sint and div_sint.
 *
 * Revision 1.12  2003/10/16 08:49:21  j_novak
 * Added a flag to decide wether the output is in the Ylm or in the standard base.
 *
 * Revision 1.11  2003/10/15 21:08:22  e_gourgoulhon
 * Added method poisson_angu.
 *
 * Revision 1.10  2003/10/15 16:03:35  j_novak
 * Added the angular Laplace operator for Scalar.
 *
 * Revision 1.9  2003/10/15 10:27:33  e_gourgoulhon
 * Classes Map, Map_af and Map_et: added new methods dsdt, stdsdp and div_tant.
 * Class Map_radial: added new Coord's : drdt and stdrdp.
 *
 * Revision 1.8  2003/06/20 14:14:53  f_limousin
 * Add the operator== to compare two Cmp.
 *
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
#include "tbl.h"
#include "itbl.h"

class Scalar ;
class Cmp ;
class Param ; 
class Map_af ; 
class Map_et ; 
class Tenseur ;
class Param_elliptic ;
class Metric_flat ; 
class Tbl ;
class Itbl ;

			//------------------------------------//
			//            class Map               //
			//------------------------------------//

/**
 * Base class for coordinate mappings. \ingroup (map)
 * 
 * This class is the basic class for describing the mapping between the
 * grid coordinates \f$(\xi, \theta', \phi')\f$ and the physical coordinates
 * \f$(r, \theta, \phi)\f$ [cf. Bonazzola, Gourgoulhon & Marck, \e Phys. 
 * \e Rev. \e D \b 58, 104020 (1998)]. 
 * The class \c Map is an abstract one: it cannot be instanciated. 
 * Specific implementation of coordinate mappings will be performed by derived
 * classes of \c Map. 
 *
 * The class \c Map and its derived classes determine the methods for 
 * partial derivatives with respect to the physical coordinates, as well
 * as resolution of basic partial differential equations (e.g. Poisson
 * equations). 
 *
 * The mapping is defined with respect to some ``absolute'' reference frame, 
 * whose Cartesian coordinates are denoted by \e (X,Y,Z). The coordinates
 * (X, Y, Z) of center of the mapping (i.e. the point \e r =0) are given by the data 
 * members \c  (ori_x,ori_y,ori_z). 
 * The Cartesian coordinate relative to the mapping (i.e. defined from 
 * \f$(r, \theta, \phi)\f$ by the usual formul\ae \f$x=r\sin\theta\cos\phi, \ldots\f$)
 * are denoted by \e (x,y,z). The planes \e (x,y) and \e (X,Y) are supposed to
 * coincide, so that the relative orientation of the mapping with respect to 
 * the absolute reference frame is described by only one angle (the data member
 * \c rot_phi).
 *
 *
 */
class Map {

    // Data : 
    // ----
    protected:
	/// Pointer on the multi-grid \c Mgd3  on which \c this is defined	
	const Mg3d* mg ; 
	 
	double ori_x ;		///< Absolute coordinate \e x  of the origin
	double ori_y ;		///< Absolute coordinate \e y of the origin
	double ori_z ;		///< Absolute coordinate \e z  of the origin
	double rot_phi ;	///< Angle between the \e x --axis and \e X --axis
    
	/** Orthonormal vectorial basis 
	 *  \f$(\partial/\partial r,1/r\partial/\partial \theta,
	 *	1/(r\sin\theta)\partial/\partial \phi)\f$
	 * associated with the coordinates \f$(r, \theta, \phi)\f$ of the 
	 * mapping.
	 */
	Base_vect_spher bvect_spher ; 

	/** Cartesian basis 
	 *  \f$(\partial/\partial x,\partial/\partial y,\partial/\partial z)\f$
	 * associated with the coordinates \e (x,y,z) of the
	 * mapping, i.e. the Cartesian coordinates related to 
	 * \f$(r, \theta, \phi)\f$ by means of usual formulae. 
	 */
	Base_vect_cart bvect_cart ; 
	
        /** Pointer onto the flat metric associated with the spherical coordinates
         *  and with components expressed in the triad \c bvect_spher 
         */
        mutable Metric_flat* p_flat_met_spher ; 

        /** Pointer onto the flat metric associated with the Cartesian coordinates
         *  and with components expressed in the triad \c bvect_cart 
         */
        mutable Metric_flat* p_flat_met_cart ; 

	/** The null Cmp. 
	 *  To be used by the \c Tenseur  class when necessary to 
	 *  return a null \c Cmp . 
	 */
	Cmp* p_cmp_zero ; 

    public:
	Coord r ;	///< \e r  coordinate centered on the grid
	Coord tet ;	///< \f$\theta\f$ coordinate centered on the grid
	Coord phi ;	///< \f$\phi\f$ coordinate centered on the grid
	Coord sint ;	///< \f$\sin\theta\f$
	Coord cost ;	///< \f$\cos\theta\f$
	Coord sinp ;	///< \f$\sin\phi\f$
	Coord cosp ;	///< \f$\cos\phi\f$

	Coord x ;	///< \e x  coordinate centered on the grid
	Coord y ;	///< \e y coordinate centered on the grid
	Coord z ;	///< \e z  coordinate centered on the grid

	Coord xa ;	///< Absolute \e x  coordinate
	Coord ya ;	///< Absolute \e y coordinate
	Coord za ;	///< Absolute \e z  coordinate
    

    // Constructors, destructor : 
    // ------------------------

    protected:
	explicit Map(const Mg3d& ) ;	///< Constructor from a multi-domain 3D grid
	Map(const Map &) ;		///< Copy constructor
	Map(const Mg3d&, FILE* ) ; ///< Constructor from a file (see \c sauve(FILE* ) )

    public:	
	virtual ~Map() ;		///< Destructor
	
    // Memory management
    // -----------------
    protected:
	virtual void reset_coord() ;  ///< Resets all the member \c Coords	
	
    // Outputs
    // -------
    private:
	virtual ostream& operator>>(ostream &) const = 0 ;    ///< Operator >>

    public:
	virtual void sauve(FILE* ) const ;	///< Save in a file
	

    // Extraction of information
    // -------------------------

    public:
	/// Gives the \c Mg3d  on which the mapping is defined
	const Mg3d* get_mg() const {return mg; };

	/// Returns the \e x  coordinate of the origin
	double get_ori_x() const {return ori_x;} ; 
	/// Returns the \e y coordinate of the origin 
	double get_ori_y() const {return ori_y;} ; 
	/// Returns the \e z  coordinate of the origin
	double get_ori_z() const {return ori_z;} ; 

	/// Returns the angle between the \e x --axis and \e X --axis
	double get_rot_phi() const {return rot_phi;} ; 
	
	/** Returns the orthonormal vectorial basis 
	 *  \f$(\partial/\partial r,1/r\partial/\partial \theta,
	 *	1/(r\sin\theta)\partial/\partial \phi)\f$
	 * associated with the coordinates \f$(r, \theta, \phi)\f$ of the 
	 * mapping.
	 */
	const Base_vect_spher& get_bvect_spher() const {return bvect_spher;} ; 

	/** Returns the Cartesian basis 
	 *  \f$(\partial/\partial x,\partial/\partial y,\partial/\partial z)\f$
	 * associated with the coordinates \e (x,y,z) of the
	 * mapping, i.e. the Cartesian coordinates related to 
	 * \f$(r, \theta, \phi)\f$ by means of usual formulae. 
	 */
	const Base_vect_cart& get_bvect_cart() const {return bvect_cart;} ; 

        /** Returns the flat metric associated with the spherical coordinates
         *  and with components expressed in the triad \c bvect_spher 
         */
        const Metric_flat& flat_met_spher() const ; 

        /** Returns the flat metric associated with the Cartesian coordinates
         *  and with components expressed in the triad \c bvect_cart 
         */
        const Metric_flat& flat_met_cart() const ; 

	/** Returns the null \c Cmp  defined on \c *this. 
	 *  To be used by the \c Tenseur  class when necessary to 
	 *  return a null \c Cmp . 
	 */
	const Cmp& cmp_zero() const {return *p_cmp_zero;} ;	
	
	/** Determines the coordinates \f$(r,\theta,\phi)\f$
	 *  corresponding to given absolute Cartesian coordinates
	 *  \e (X,Y,Z). 
	 *	@param xx [input] value of the coordinate \e x  (absolute frame)
	 *	@param yy [input] value of the coordinate \e y (absolute frame)
	 *	@param zz [input] value of the coordinate \e z  (absolute frame)
	 *	@param rr [output] value of \e r 
	 *	@param theta [output] value of \f$\theta\f$
	 *	@param pphi [output] value of \f$\phi\f$
	 */
	void convert_absolute(double xx, double yy, double zz, 
			      double& rr, double& theta, double& pphi) const ; 

	/** Returns the value of the radial coordinate \e r  for a given
	 *  \f$(\xi, \theta', \phi')\f$ in a given domain. 
	 *	@param l [input] index of the domain
	 *	@param xi [input] value of \f$\xi\f$
	 *	@param theta [input] value of \f$\theta'\f$
	 *	@param pphi [input] value of \f$\phi'\f$
	 *	@return value of \f$r=R_l(\xi, \theta', \phi')\f$
	 */
	virtual double val_r(int l, double xi, double theta, double pphi) 
			     const = 0 ; 
	
	/** Computes the domain index \e l  and the value of \f$\xi\f$ corresponding
	 * to a point given by its physical coordinates \f$(r, \theta, \phi)\f$. 
	 *	@param rr [input] value of \e r 
	 *	@param theta [input] value of \f$\theta\f$
	 *	@param pphi [input] value of \f$\phi\f$
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of \f$\xi\f$
	 */
	virtual void val_lx(double rr, double theta, double pphi, 
			    int& l, double& xi) const = 0 ; 
	
	/** Computes the domain index \e l  and the value of \f$\xi\f$ corresponding
	 * to a point given by its physical coordinates \f$(r, \theta, \phi)\f$. 
	 * This version enables to pass some parameters to control the
	 * accuracy of the computation. 
	 *	@param rr [input] value of \e r 
	 *	@param theta [input] value of \f$\theta\f$
	 *	@param pphi [input] value of \f$\phi\f$
	 *	@param par [input/output] parameters to control the
	 *	    accuracy of the computation
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of \f$\xi\f$
	 */
	virtual void val_lx(double rr, double theta, double pphi, 
			    const Param& par, int& l, double& xi) const = 0 ;


	/// Comparison operator (egality)
	virtual bool operator==(const Map& ) const = 0;  

 
	
    // Modification of the origin, the orientation and the radial scale:
    // ----------------------------------------------------------------
    public:
	void set_ori(double xa0, double ya0, double za0) ;  ///< Sets a new origin
	void set_rot_phi(double phi0) ;		///< Sets a new rotation angle

	/** Sets a new radial scale.
	 *	@param lambda [input] factor by which the value of \e r  is to
	 *	    be multiplied
	 */
	virtual void homothetie(double lambda) = 0 ;	

	/** Rescales the outer boundary of one domain.
	 *  The inner boundary is unchanged. The inner boundary 
	 *  of the next domain is changed to match the new outer
	 *  boundary. 
	 *	@param l [input] index of the domain
	 *	@param lambda [input] factor by which the value of 
	 *	    \f$R(\theta, \varphi)\f$ defining the outer boundary
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

	/** Recomputes the values of a \c Cmp  at the collocation points 
	 *  after a change in the mapping. 
	 *  @param mp_prev [input] Previous value of the mapping. 
	 *  @param nzet [input] Number of domains where the computation 
	 *    must be done: the computation is performed for the domains
	 *    of index \f$0\le {\tt l} \le {\tt nzet-1}\f$; \c uu  is set
	 *    to zero in the other domains. 
	 *  @param uu [input/output] input : \c Cmp  previously computed on 
	 *			     the mapping \c *mp_prev ; ouput :
	 *			     values of (logically) the same \c Cmp
	 *			     at the grid points defined by \c *this. 
	 */
	virtual void reevaluate(const Map* mp_prev, int nzet, 
				Cmp& uu) const = 0 ; 

	/** Recomputes the values of a \c Cmp  at the collocation points 
	 *  after a change in the mapping. 
	 *  Case where the \c Cmp  is symmetric with respect the plane y=0.
	 *  @param mp_prev [input] Previous value of the mapping. 
	 *  @param nzet [input] Number of domains where the computation 
	 *    must be done: the computation is performed for the domains
	 *    of index \f$0\le {\tt l} \le {\tt nzet-1}\f$; \c uu is set
	 *    to zero in the other domains. 
	 *  @param uu [input/output] input : \c Cmp  previously computed on 
	 *			     the mapping \c *mp_prev ; ouput :
	 *			     values of (logically) the same \c Cmp 
	 *			     at the grid points defined by \c *this. 
	 */
	virtual void reevaluate_symy(const Map* mp_prev, int nzet, 
				Cmp& uu) const = 0 ; 



    // Differential operators:
    // ----------------------
    public:
	/** Computes \f$\partial/ \partial r\f$ of a \c Cmp .
	 *  Note that in the  compactified external domain (CED), it computes
	 *  \f$-\partial/ \partial u = r^2 \partial/ \partial r\f$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of \c ci 
	 */
	virtual void dsdr(const Cmp& ci, Cmp& resu) const = 0 ;  	    

	/** Computes \f$1/r \partial/ \partial \theta\f$ of a \c Cmp .
	 *  Note that in the  compactified external domain (CED), it computes
	 *  \f$1/u \partial/ \partial \theta = r \partial/ \partial \theta\f$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of \c ci 
	 */
	virtual void srdsdt(const Cmp& ci, Cmp& resu) const = 0 ;  	    
	
	/** Computes \f$1/(r\sin\theta) \partial/ \partial \phi\f$ of a \c Cmp .
	 *  Note that in the  compactified external domain (CED), it computes
	 *  \f$1/(u\sin\theta) \partial/ \partial \phi = 
	 *	    r/\sin\theta \partial/ \partial \phi\f$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of \c ci 
	 */
	virtual void srstdsdp(const Cmp& ci, Cmp& resu) const = 0 ;  	    
    
	/** Computes \f$\partial/ \partial r\f$ of a \c Scalar .
	 *  Note that in the  compactified external domain (CED), the \c dzpuis 
	 *  flag of the output is 2 if the input has \c dzpuis  = 0, and 
	 *  is increased by 1 in other cases.
	 *  @param uu [input] field to consider
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void dsdr(const Scalar& uu, Scalar& resu) const = 0 ;
	
	/** Computes \f$\partial/ \partial r\f$ of a \c Scalar if the description is affine and 
	 * \f$\partial/ \partial \ln r\f$ if it is logarithmic.
	 *  Note that in the  compactified external domain (CED), the \c dzpuis 
	 *  flag of the output is 2 if the input has \c dzpuis  = 0, and 
	 *  is increased by 1 in other cases.
	 *  @param uu [input] field to consider
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void dsdradial(const Scalar& uu, Scalar& resu) const = 0 ;  

	/** Computes \f$1/r \partial/ \partial \theta\f$ of a \c Scalar .
	 *  Note that in the  compactified external domain (CED), the \c dzpuis 
	 *  flag of the output is 2 if the input has \c dzpuis  = 0, and 
	 *  is increased by 1 in other cases.
	 *  @param uu [input] field to consider
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void srdsdt(const Scalar& uu, Scalar& resu) const = 0 ;  	    
	
	/** Computes \f$1/(r\sin\theta) \partial/ \partial \phi\f$ of a \c Scalar .
	 *  Note that in the  compactified external domain (CED), the \c dzpuis 
	 *  flag of the output is 2 if the input has \c dzpuis  = 0, and 
	 *  is increased by 1 in other cases.
	 *  @param uu [input] field to consider
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void srstdsdp(const Scalar& uu, Scalar& resu) const = 0 ;  	    
    
	/** Computes \f$\partial/ \partial \theta\f$ of a \c Scalar .
	 *  @param uu [input] scalar field 
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void dsdt(const Scalar& uu, Scalar& resu) const = 0 ;  	    

	/** Computes \f$1/\sin\theta \partial/ \partial \varphi\f$ of a \c Scalar .
	 *  @param uu [input] scalar field 
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void stdsdp(const Scalar& uu, Scalar& resu) const = 0 ;  	    

	/** Computes the Laplacian of a scalar field.
	 *   @param uu	[input]  Scalar field \e u  (represented as a \c Cmp )
	 *			 the Laplacian \f$\Delta u\f$ of which is to be computed
	 *   @param zec_mult_r [input] Determines the quantity computed in
	 *			 the  compactified external domain (CED) :  \\
	 *		    zec_mult_r = 0 : \f$\Delta u\f$	\\
	 *		    zec_mult_r = 2 : \f$r^2 \,  \Delta u\f$	\\
	 *		    zec_mult_r = 4 (default) : \f$r^4 \, \Delta u\f$	
	 *   @param lap [output] Laplacian of \e u 
	 */
	virtual void laplacien(const Cmp& uu, int zec_mult_r, 
			       Cmp& lap) const = 0 ; 
	
	/** Computes the angular Laplacian of a scalar field.
	 *   @param uu	[input]  Scalar field \e u  (represented as a \c Scalar )
	 *			 the Laplacian \f$\Delta u\f$ of which is to be computed
	 *   @param lap [output] Angular Laplacian of \e u  (see documentation 
	 *                       of \c Scalar 
	 */
	virtual void lapang(const Scalar& uu, Scalar& lap) const = 0 ; 


	/** Computes the radial primitive which vanishes for \f$r\to \infty\f$.
	 *  i.e. the function 
	 *   \f$ F(r,\theta,\varphi) = \int_r^\infty f(r',\theta,\varphi) \, dr' \f$
	 *
	 *      @param uu [input] function \e f (must have a \c dzpuis = 2)
	 *      @param resu [input] function \e F
	 *      @param null_infty if true (default), the primitive is null
	 *      at infinity (or on the grid boundary). \e F is null at the
	 *      center otherwise
	 */ 
	virtual void primr(const Scalar& uu, Scalar& resu, 
			   bool null_infty) const = 0 ;  	    

	 
    // Various linear operators
    // ------------------------
    public: 
	/** Multiplication by \e r  of a \c Scalar , the \c dzpuis  
	 *  of \c uu  is not changed.
	 */
	virtual void mult_r(Scalar& uu) const = 0 ; 

	/** Multiplication by \e r  of a \c Cmp . In the CED, 
	 *  there is only a decrement of \c dzpuis 
	 */
	virtual void mult_r(Cmp& ci) const = 0 ; 

	/** Multiplication by \e r  (in the  compactified external domain only)
	 * of a \c Scalar 
	 */
	virtual void mult_r_zec(Scalar& ) const = 0 ;

	/** Multiplication by \f$r\sin\theta\f$ of a \c Scalar 
	 */
	virtual void mult_rsint(Scalar& ) const = 0 ; 

	/** Division by \f$r\sin\theta\f$ of a \c Scalar
	 */
	virtual void div_rsint(Scalar& ) const = 0 ; 

	/** Division by \e r  of a \c Scalar 
	 */
	virtual void div_r(Scalar& ) const = 0 ; 

	/** Division by \e r  (in the  compactified external domain only)
	 * of a \c Scalar 
	 */
	virtual void div_r_zec(Scalar& ) const = 0 ;

	/** Multiplication by \f$\cos\theta\f$ of a \c Scalar 
	 */
	virtual void mult_cost(Scalar& ) const = 0 ; 

	/** Multiplication by \f$\sin\theta\f$ of a \c Scalar 
	 */
	virtual void mult_sint(Scalar& ) const = 0 ; 

	/** Division by \f$\sin\theta\f$ of a \c Scalar 
	 */
	virtual void div_sint(Scalar& ) const = 0 ; 

	/** Division by \f$\tan\theta\f$ of a \c Scalar 
	 */
	virtual void div_tant(Scalar& ) const = 0 ; 

	/** Computes the Cartesian x component (with respect to 
	 *  \c bvect_cart ) of a vector given
	 *  by its spherical components with respect to \c bvect_spher .
	 *  
	 *  @param v_r [input] \e r -component of the vector 
	 *  @param v_theta [input] \f$\theta\f$-component of the vector 
	 *  @param v_phi [input] \f$\phi\f$-component of the vector 
	 *  @param v_x [output] x-component of the vector 
	 */
	virtual void comp_x_from_spherical(const Cmp& v_r, const Cmp& v_theta, 
					   const Cmp& v_phi, Cmp& v_x) const = 0 ; 
	
	/** Computes the Cartesian y component (with respect to 
	 *  \c bvect_cart ) of a vector given
	 *  by its spherical components with respect to \c bvect_spher .
	 *  
	 *  @param v_r [input] \e r -component of the vector 
	 *  @param v_theta [input] \f$\theta\f$-component of the vector 
	 *  @param v_phi [input] \f$\phi\f$-component of the vector 
	 *  @param v_y [output] y-component of the vector 
	 */
	virtual void comp_y_from_spherical(const Cmp& v_r, const Cmp& v_theta, 
					   const Cmp& v_phi, Cmp& v_y) const = 0 ; 
	
	/** Computes the Cartesian z component (with respect to 
	 *  \c bvect_cart ) of a vector given
	 *  by its spherical components with respect to \c bvect_spher .
	 *  
	 *  @param v_r [input] \e r -component of the vector 
	 *  @param v_theta [input] \f$\theta\f$-component of the vector 
	 *  @param v_z [output] z-component of the vector 
	 */
	virtual void comp_z_from_spherical(const Cmp& v_r, const Cmp& v_theta, 
					   Cmp& v_z) const = 0 ; 
	
	 /** Computes the Spherical r component (with respect to 
	 *  \c bvect_spher ) of a vector given
	 *  by its cartesian components with respect to \c bvect_cart .
	 *  
	 *  @param v_x [input] x-component of the vector 
	 *  @param v_y [input] y-component of the vector 
	 *  @param v_z [input] z-component of the vector 
	 *  @param v_r [output] \e r -component of the vector 
	 */
	virtual void comp_r_from_cartesian(const Cmp& v_x, const Cmp& v_y, 
					   const Cmp& v_z, Cmp& v_r) const = 0 ; 
	
	/** Computes the Spherical \f$\theta\f$ component (with respect to 
	 *  \c bvect_spher ) of a vector given
	 *  by its cartesian components with respect to \c bvect_cart .
	 *  
	 *  @param v_x [input] x-component of the vector 
	 *  @param v_y [input] y-component of the vector 
	 *  @param v_z [input] z-component of the vector 
	 *  @param v_t [output] \f$\theta\f$-component of the vector 
	 */
	virtual void comp_t_from_cartesian(const Cmp& v_x, const Cmp& v_y, 
					   const Cmp& v_z, Cmp& v_t) const = 0 ; 
	
	/** Computes the Spherical \f$\phi\f$ component (with respect to 
	 *  \c bvect_spher ) of a vector given
	 *  by its cartesian components with respect to \c bvect_cart .
	 *  
	 *  @param v_x [input] x-component of the vector 
	 *  @param v_y [input] y-component of the vector 
	 *  @param v_p [output] \f$\phi\f$-component of the vector 
	 */
	virtual void comp_p_from_cartesian(const Cmp& v_x, const Cmp& v_y, 
					   Cmp& v_p) const = 0 ; 
	
	/** Decreases by 1 the value of \c dzpuis  of a \c Scalar  
	 *  and changes accordingly its values in the  
	 *  compactified external domain (CED).
	 */
	virtual void dec_dzpuis(Scalar& ) const = 0 ; 
	
	/** Decreases by 2 the value of \c dzpuis  of a \c Scalar  
	 *  and changes accordingly its values in the  
	 *  compactified external domain (CED).
	 */
	virtual void dec2_dzpuis(Scalar& ) const = 0 ; 
	
	/** Increases by 1 the value of \c dzpuis  of a \c Scalar  
	 *  and changes accordingly its values in the  
	 *  compactified external domain (CED).
	 */
	virtual void inc_dzpuis(Scalar& ) const = 0 ; 
	
	/** Increases by 2 the value of \c dzpuis  of a \c Scalar  
	 *  and changes accordingly its values in the  
	 *  compactified external domain (CED).
	 */
	virtual void inc2_dzpuis(Scalar& ) const = 0 ; 
	
	/** Computes the integral over all space of a \c Cmp .
	 *  The computed quantity is 
	 *    \f$\int u \, r^2 \sin\theta \,  dr\, d\theta \, d\phi\f$.
	 *  The routine allocates a \c Tbl  (size: \c mg->nzone ) to store 
	 *  the result (partial integral) in each domain and returns a pointer 
	 *  to it.
	 */
	virtual Tbl* integrale(const Cmp&) const = 0 ; 
	
    // PDE resolution :
    // --------------
    public:
	/** Computes the solution of a scalar Poisson equation.
	 * 
	 *   @param source [input] source \f$\sigma\f$ of the Poisson equation 
	 *	    \f$\Delta u = \sigma\f$.
	 *   @param par [input/output] possible parameters to control the
	 *   resolution of the Poisson equation. See the actual implementation 
	 *   in the derived class of \c Map for documentation. 
	 *   @param uu [input/output] solution \e u  with the boundary condition 
	 *	    \e u =0 at spatial infinity. 
	 */
	virtual void poisson(const Cmp& source, Param& par, Cmp& uu) const = 0 ;
	virtual void poisson_falloff(const Cmp& source, Param& par, Cmp& uu, int k_falloff) const = 0 ;

	/** Computes the solution of a scalar Poisson equation.
	 *   The regularized source
	 *
	 *   @param source [input] source \f$\sigma\f$ of the Poisson equation 
	 *	    \f$\Delta u = \sigma\f$.
	 *   @param k_div [input] regularization degree of the procedure
	 *   @param nzet [input] number of domains covering the star
	 *   @param unsgam1 [input] parameter \f$1/(\gamma-1)\f$ where \f$\gamma\f$
	 *          denotes the adiabatic index.
	 *   @param par [input/output] possible parameters to control the
	 *   resolution of the Poisson equation. See the actual implementation 
	 *   in the derived class of \c Map for documentation.
	 *   @param uu [input/output] solution \e u  with the boundary condition 
	 *	    \e u =0 at spatial infinity.
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
	 * \f$ a \Delta\psi + {\bf b} \cdot \nabla \psi = \sigma\f$.
	 * 
	 * @param source [input] source \f$\sigma\f$ of the above equation
	 * @param aa [input] factor \e a  in the above equation
	 * @param bb [input] vector \e \b b in the above equation
	 * @param par [input/output] possible parameters to control the
	 *   resolution of the equation. See the actual implementation 
	 *   in the derived class of \c Map for documentation. 
	 * @param psi [input/output] solution \f$\psi\f$ which satisfies
	 *	\f$\psi(0)=0\f$.  
	 */ 
	virtual void poisson_compact(const Cmp& source, const Cmp& aa, 
				     const Tenseur& bb, const Param& par, 
				     Cmp& psi) const = 0 ;

	/** Computes the solution of an angular Poisson equation.
	 * The angular Poisson equation is \f$\Delta_{\theta\varphi} u = \sigma\f$,
	 * where \f$\Delta_{\theta\varphi} u := \frac{\partial^2 u}
	 *  {\partial \theta^2} + \frac{1}{\tan \theta} \frac{\partial u}
	 *  {\partial \theta} +\frac{1}{\sin^2 \theta}\frac{\partial^2 u}
	 *  {\partial \varphi^2}\f$.
	 * 
	 *   @param source [input] source \f$\sigma\f$ of the Poisson equation 
	 *	    \f$\Delta_{\theta\varphi} u = \sigma\f$.
	 *   @param par [input/output] possible parameters to control the
	 *   resolution of the Poisson equation. See the actual implementation 
	 *   in the derived class of \c Map for documentation. 
	 *   @param uu [input/output] solution \e u  
	 */
	virtual void poisson_angu(const Scalar& source, Param& par, 
					Scalar& uu) const = 0 ;

	
    public:
	/**
	 * Function intended to be used by \c Map::poisson_vect 
	 * and \c Map::poisson_vect_oohara . It constructs the sets of 
	 * parameters used for each scalar Poisson equation from the one for 
	 * the vectorial one.
	 * 
	 * @param para [input] : the \c Param  used for the resolution of 
	 * the vectorial Poisson equation : \\
	 * \c para.int()  maximum number of iteration.\\
	 * \c para.double(0)  relaxation parameter.\\
	 * \c para.double(1)  required precision. \\
	 * \c para.tenseur_mod()  source of the vectorial part at the previous 
	 * step.\\
	 * \c para.cmp_mod()  source of the scalar part at the previous 
	 * step.
	 * 
	 * @param i [input] number of the scalar Poisson equation that is being 
	 * solved (values from 0 to 2 for the componants of the vectorial part
	 * and 3 for the scalar one).
	 * 
	 * @return the pointer on the parameter set used for solving the
	 * scalar Poisson equation labelled by \e i .
	 */
	virtual Param* donne_para_poisson_vect (Param& para, int i) const = 0;
	
	/**
	 * Computes the solution of a Poisson equation from the domain 
	 * \c num_front+1 .
	 * imposing a boundary condition at the boundary between the domains 
	 * \c num_front  and \c num_front+1 .
	 * 
	 * @param source [input] : source of the equation.
	 * @param limite [input] : \c limite[num_front]  contains the angular 
	 * function being the boudary condition.
	 * @param raccord [input] : 1 for the Dirichlet problem and 2 for 
	 * the Neumann one and 3 for Dirichlet-Neumann.
	 * @param num_front [input] : index of the boudary at which the boundary 
	 * condition has to be imposed.
	 * @param pot [output] : result.
	 * @param fact_dir [input] : Valeur by which we multiply the quantity 
	 * we solve. (in the case of Dirichlet-Neumann boundary condition.)
	 * @param fact_neu [input] : Valeur by which we multiply the radial 
	 * derivative of the quantity we solve.	(in the case of 
	 * Dirichlet-Neumann boundary condition.)
	 */
	virtual void poisson_frontiere (const Cmp& source,const Valeur& limite,
					int raccord, int num_front, Cmp& pot, 
					double = 0., double = 0.) const = 0 ;
	
	virtual void poisson_frontiere_double (const Cmp& source, const Valeur& lim_func,
			const Valeur& lim_der, int num_zone, Cmp& pot) const = 0 ;
	

	/**
	 * Computes the solution of a Poisson equation in the shell,
	 * imposing a boundary condition at the surface of the star 
	 * 
	 * @param source [input] : source of the equation.
	 * @param limite [input] : \c limite[num_front]  contains the angular 
	 * function being the boudary condition.
	 * @param par [input] : parameters of the computation.
	 * @param pot [output] : result.
	 */
	virtual void poisson_interne (const Cmp& source, const Valeur& limite,
			 Param& par, Cmp& pot) const = 0 ;


	/** Computes the solution of a 2-D Poisson equation.
	 *  The 2-D Poisson equation writes
	 *  \f[
	 *	{\partial^2 u\over\partial r^2} + 
	 *	    {1\over r} {\partial u \over \partial r} + 
	 *	    {1\over r^2} {\partial^2 u\over\partial \theta^2} = 
	 *		\sigma \ . 
	 *  \f] 
	 *
	 *   @param source_mat [input] Compactly supported part of 
	 *	    the source \f$\sigma\f$ of the 2-D Poisson equation (typically
	 *	    matter terms)
	 *   @param source_quad [input] Non-compactly supported part of 
	 *	    the source \f$\sigma\f$ of the 2-D Poisson equation (typically
	 *	    quadratic terms)
	 *   @param par [input/output] possible parameters to control the
	 *   resolution of the Poisson equation. See the actual implementation 
	 *   in the derived class of \c Map for documentation. 
	 *   @param uu [input/output] solution \e u  with the boundary condition 
	 *	    \e u =0 at spatial infinity. 
	 */
	virtual void poisson2d(const Cmp& source_mat, const Cmp& source_quad, 
			       Param& par, Cmp& uu) const = 0 ;

	/** Performs one time-step integration of the d'Alembert scalar equation
	 *   @param par [input/output] possible parameters to control the
	 *   resolution of the Poisson equation. See the actual implementation 
	 *   in the derived class of \c Map for documentation. Note that, 
	 *   at least, param must contain the time step as first \c double  
	 *   parameter.
	 *   @param fJp1 [output] solution \f$f^{J+1}\f$ at time \e J+1
	 *   with boundary conditions of outgoing radiation (not exact!)
	 *   @param fJ [input] solution \f$f^J\f$ at time \e J 
	 *   @param fJm1 [input] solution \f$f^{J-1}\f$ at time \e J-1 
	 *   @param source [input] source \f$\sigma\f$ of the d'Alembert equation 
	 *	    \f$\diamond u = \sigma\f$.
	 */
	virtual void dalembert(Param& par, Scalar& fJp1, const Scalar& fJ, 
			       const Scalar& fJm1, const Scalar& source) const = 0 ;

    // Friend functions : 
    // ----------------
	friend ostream& operator<<(ostream& , const Map& ) ;	///< Operator <<
};
ostream& operator<<(ostream& , const Map& ) ; 



			//------------------------------------//
			//          class Map_radial          //
			//------------------------------------//



/**
 * Base class for pure radial mappings. \ingroup (map)
 * 
 * A pure radial mapping is a mapping of the type \f$r=R(\xi, \theta', \phi')\f$, 
 * \f$\theta=\theta'\f$, \f$\phi=\phi'\f$. 
 * The class \c Map_radial is an abstract class. Effective implementations 
 * of radial mapping are performed by the derived class \c Map_af  and 
 * \c Map_et .
 * 
 * 
 */

class Map_radial : public Map {

    // Data : 
    // ----

    // 0th order derivatives of the mapping
    // - - - - - - - - - - - - - - - - - - 
    public:
	/**
	 * \f$\xi/R\f$ in the nucleus; \\
	 * \e 1/R  in the non-compactified shells; \\
	 * \f$(\xi-1)/U\f$ in the compactified outer domain.
	 */
	Coord xsr ;	    
	
    // 1st order derivatives of the mapping
    // - - - - - - - - - - - - - - - - - - 
    public:
	/**
	 * \f$1/(\partial R/\partial\xi) = \partial \xi /\partial r\f$ in the nucleus
	 *   and in the non-compactified shells; \\
	 * \f$-1/(\partial U/\partial\xi) = - \partial \xi /\partial u\f$ in the 
	 *   compactified outer domain.
	 */
	Coord dxdr ;	    
	
	/**
	 * \f$\partial R/\partial\theta'\f$ in the nucleus
	 *   and in the non-compactified shells; \\
	 * \f$-\partial U/\partial\theta'\f$ in the 
	 *   compactified external domain (CED).
	 */
	Coord drdt ; 

	/**
	 * \f${1\over\sin\theta} \partial R/\partial\varphi'\f$ in the nucleus
	 *   and in the non-compactified shells; \\
	 * \f$-{1\over\sin\theta}\partial U/\partial\varphi'\f$ in the 
	 *   compactified external domain (CED).
	 */
	Coord stdrdp ; 

	/**
	 * \f$1/R \times (\partial R/\partial\theta')\f$ in the nucleus
	 *   and in the non-compactified shells; \\
	 * \f$-1/U \times (\partial U/\partial\theta)\f$ in the 
	 *   compactified outer domain.
	 */
	Coord srdrdt ;	    

	/**
	 * \f$1/(R\sin\theta) \times (\partial R/\partial\varphi')\f$ in the nucleus
	 *   and in the non-compactified shells; \\
	 * \f$-1/(U\sin\theta) \times (\partial U/\partial\varphi')\f$ in the 
	 *   compactified outer domain.
	 */
	Coord srstdrdp ;	    

	/**
	 * \f$1/R^2 \times (\partial R/\partial\theta')\f$ in the nucleus
	 *   and in the non-compactified shells; \\
	 * \f$-1/U^2 \times (\partial U/\partial\theta')\f$ in the 
	 *   compactified outer domain.
	 */
	Coord sr2drdt ;	    
	
	/**
	 * \f$1/(R^2\sin\theta) \times (\partial R/\partial\varphi')\f$ in the nucleus
	 *   and in the non-compactified shells; \\
	 * \f$-1/(U^2\sin\theta) \times (\partial U/\partial\varphi')\f$ in the 
	 *   compactified outer domain.
	 */
	Coord sr2stdrdp ;   

    // 2nd order derivatives of the mapping
    // - - - - - - - - - - - - - - - - - - 
    public:
	/**
	 * \f$\partial^2 R/\partial\xi^2\f$ in the nucleus
	 *   and in the non-compactified shells; \\
	 * \f$-\partial^2 U/\partial\xi^2 \f$ in the 
	 *   compactified outer domain.
	 */
	Coord d2rdx2 ;	    
	
	/**
	 * \f$1/R^2 \times [ 1/\sin(\theta)\times \partial /\partial\theta'
	 *   (\sin\theta \partial R /\partial\theta') + 1/\sin^2\theta
	 *   \partial^2 R /\partial{\varphi'}^2] \f$ in the nucleus
	 *   and in the non-compactified shells; \\
	 * \f$- 1/U^2 \times [ 1/\sin(\theta)\times \partial /\partial\theta'
	 *   (\sin\theta \partial U /\partial\theta') + 1/\sin^2\theta
	 *   \partial^2 U /\partial{\varphi'}^2] \f$ in the 
	 *   compactified outer domain.
	 */
	Coord lapr_tp ;	    
	

	/**
	 * \f$\partial^2 R/\partial\xi\partial\theta'\f$ in the nucleus
	 *   and in the non-compactified shells; \\
	 * \f$-\partial^2 U/\partial\xi\partial\theta'\f$ in the 
	 *   compactified outer domain.
	 */
	Coord d2rdtdx ;	    
	
	/**
	 * \f$1/\sin\theta \times \partial^2 R/\partial\xi\partial\varphi'\f$ 
	 * in the nucleus and in the non-compactified shells; \\
	 * \f$-1/\sin\theta \times \partial^2 U/\partial\xi\partial\varphi' \f$ 
	 * in the compactified outer domain.
	 */
	Coord sstd2rdpdx ;  


	/**
	 * \f$1/R^2 \partial^2 R/\partial{\theta'}^2\f$ in the nucleus
	 *   and in the non-compactified shells; \\
	 * \f$-1/U^2 \partial^2 U/\partial{\theta'}^2\f$ in the 
	 *   compactified outer domain.
	 */
	Coord sr2d2rdt2 ;
	

    // Constructors, destructor : 
    // ------------------------

    protected:
	/// Constructor from a grid (protected to make \c Map_radial an abstract class) 
	Map_radial(const Mg3d& mgrid ) ;		    
	Map_radial(const Map_radial& mp) ;   ///< Copy constructor
	Map_radial (const Mg3d&, FILE* ) ; ///< Constructor from a file (see \c sauve(FILE* ) )

    public: 
	virtual ~Map_radial() ;		    ///< Destructor

    // Memory management
    // -----------------
    protected:
	virtual void reset_coord() ;  ///< Resets all the member \c Coords	
    // Modification of the mapping
    // ---------------------------
    public:
	/// Assignment to an affine mapping. 
	virtual void operator=(const Map_af& ) = 0 ;
		
    // Outputs
    // -------
    public:  
 	virtual void sauve(FILE* ) const ;	///< Save in a file
   
    // Extraction of information
    // -------------------------
	/** Returns the value of the radial coordinate \e r  for a given
	 *  \f$\xi\f$ and a given collocation point in \f$(\theta', \phi')\f$ 
	 *   in a given domain. 
	 *	@param l [input] index of the domain
	 *	@param xi [input] value of \f$\xi\f$
	 *	@param j [input] index of the collocation point in \f$\theta'\f$
	 *	@param k [input] index of the collocation point in \f$\phi'\f$
	 *	@return value of \f$r=R_l(\xi, {\theta'}_j, {\phi'}_k)\f$
	 */
	virtual double val_r_jk(int l, double xi, int j, int k) const = 0 ; 
	
	/** Computes the domain index \e l  and the value of \f$\xi\f$ corresponding
	 * to a point of arbitrary \e r  but collocation values of \f$(\theta, \phi)\f$ 
	 *	@param rr [input] value of \e r 
	 *	@param j [input] index of the collocation point in \f$\theta\f$
	 *	@param k [input] index of the collocation point in \f$\phi\f$
	 *	@param par [input/output] parameters to control the
	 *	    accuracy of the computation
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of \f$\xi\f$
	 */
	virtual void val_lx_jk(double rr, int j, int k, const Param& par, 
			       int& l, double& xi) const = 0 ; 

	/// Comparison operator (egality)
	virtual bool operator==(const Map& ) const = 0;  

    // Values of a Cmp at the new grid points
    // --------------------------------------
	/** Recomputes the values of a \c Cmp at the collocation points 
	 *  after a change in the mapping. 
	 *  @param mp_prev [input] Previous value of the mapping. 
	 *  @param nzet [input] Number of domains where the computation 
	 *    must be done: the computation is performed for the domains
	 *    of index \f$0\le {\tt l} \le {\tt nzet-1}\f$; \c uu  is set
	 *    to zero in the other domains. 
	 *  @param uu [input/output] input : \c Cmp previously computed on 
	 *			     the mapping \c *mp_prev ; ouput :
	 *			     values of (logically) the same \c Cmp
	 *			     at the grid points defined by \c *this. 
	 */
	virtual void reevaluate(const Map* mp_prev, int nzet, Cmp& uu) const ; 

	/** Recomputes the values of a \c Cmp at the collocation points 
	 *  after a change in the mapping. 
	 *  Case where the \c Cmp is symmetric with respect to the plane y=0.
	 *  @param mp_prev [input] Previous value of the mapping. 
	 *  @param nzet [input] Number of domains where the computation 
	 *    must be done: the computation is performed for the domains
	 *    of index \f$0\le {\tt l} \le {\tt nzet-1}\f$; \c uu  is set
	 *    to zero in the other domains. 
	 *  @param uu [input/output] input : \c Cmp previously computed on 
	 *			     the mapping \c *mp_prev ; ouput :
	 *			     values of (logically) the same \c Cmp
	 *			     at the grid points defined by \c *this. 
	 */
	virtual void reevaluate_symy(const Map* mp_prev, int nzet, Cmp& uu) 
				    const ; 

    
    // Various linear operators
    // ------------------------
    public: 
	/** Multiplication by \e r  of a \c Scalar, the \c dzpuis  
	 *  of \c uu  is not changed.
	 */
	virtual void mult_r(Scalar& uu) const ; 

	/** Multiplication by \e r  of a \c Cmp. In the CED, 
	 *  there is only a decrement of \c dzpuis 
	 */
	virtual void mult_r(Cmp& ci) const ; 

	/**
	 * Multiplication by \e r  (in the compactified external domain only)
	 * of a \c Scalar
	 */
	virtual void mult_r_zec(Scalar& ) const ;

	/** Multiplication by \f$r\sin\theta\f$ of a \c Scalar
	 */
	virtual void mult_rsint(Scalar& ) const ; 

	/** Division by \f$r\sin\theta\f$ of a \c Scalar
	 */
	virtual void div_rsint(Scalar& ) const ; 

	/** Division by \e r  of a \c Scalar
	 */
	virtual void div_r(Scalar& ) const ; 

	/** Division by \e r  (in the  compactified external domain only)
	 * of a \c Scalar
	 */
	virtual void div_r_zec(Scalar& ) const ;

	/** Multiplication by \f$\cos\theta\f$ of a \c Scalar
	 */
	virtual void mult_cost(Scalar& ) const ; 

	/** Multiplication by \f$\sin\theta\f$ of a \c Scalar
	 */
	virtual void mult_sint(Scalar& ) const ; 

	/** Division by \f$\sin\theta\f$ of a \c Scalar
	 */
	virtual void div_sint(Scalar& ) const  ; 

	/** Division by \f$\tan\theta\f$ of a \c Scalar
	 */
	virtual void div_tant(Scalar& ) const  ; 

	/** Computes the Cartesian x component (with respect to 
	 *  \c bvect_cart) of a vector given
	 *  by its spherical components with respect to \c bvect_spher.
	 *  
	 *  @param v_r [input] \e r -component of the vector 
	 *  @param v_theta [input] \f$\theta\f$-component of the vector 
	 *  @param v_phi [input] \f$\phi\f$-component of the vector 
	 *  @param v_x [output] x-component of the vector 
	 */
	virtual void comp_x_from_spherical(const Cmp& v_r, const Cmp& v_theta, 
					   const Cmp& v_phi, Cmp& v_x) const ; 
	
	/** Computes the Cartesian y component (with respect to 
	 *  \c bvect_cart ) of a vector given
	 *  by its spherical components with respect to \c bvect_spher .
	 *  
	 *  @param v_r [input] \e r -component of the vector 
	 *  @param v_theta [input] \f$\theta\f$-component of the vector 
	 *  @param v_phi [input] \f$\phi\f$-component of the vector 
	 *  @param v_y [output] y-component of the vector 
	 */
	virtual void comp_y_from_spherical(const Cmp& v_r, const Cmp& v_theta, 
					   const Cmp& v_phi, Cmp& v_y) const ; 
	
	/** Computes the Cartesian z component (with respect to 
	 *  \c bvect_cart ) of a vector given
	 *  by its spherical components with respect to \c bvect_spher .
	 *  
	 *  @param v_r [input] \e r -component of the vector 
	 *  @param v_theta [input] \f$\theta\f$-component of the vector 
	 *  @param v_z [output] z-component of the vector 
	 */
	virtual void comp_z_from_spherical(const Cmp& v_r, const Cmp& v_theta, 
					   Cmp& v_z) const ; 
	
	/** Computes the Spherical r component (with respect to 
	 *  \c bvect_spher ) of a vector given
	 *  by its cartesian components with respect to \c bvect_cart .
	 *  
	 *  @param v_x [input] x-component of the vector 
	 *  @param v_y [input] y-component of the vector 
	 *  @param v_z [input] z-component of the vector 
	 *  @param v_r [output] \e r -component of the vector 
	 */
	virtual void comp_r_from_cartesian(const Cmp& v_x, const Cmp& v_y, 
					   const Cmp& v_z, Cmp& v_r) const ; 
	
	/** Computes the Spherical \f$\theta\f$ component (with respect to 
	 *  \c bvect_spher ) of a vector given
	 *  by its cartesian components with respect to \c bvect_cart .
	 *  
	 *  @param v_x [input] x-component of the vector 
	 *  @param v_y [input] y-component of the vector 
	 *  @param v_z [input] z-component of the vector 
	 *  @param v_t [output] \f$\theta\f$-component of the vector 
	 */
	virtual void comp_t_from_cartesian(const Cmp& v_x, const Cmp& v_y, 
					   const Cmp& v_z, Cmp& v_t) const ; 
	
	/** Computes the Spherical \f$\phi\f$ component (with respect to 
	 *  \c bvect_spher ) of a vector given
	 *  by its cartesian components with respect to \c bvect_cart .
	 *  
	 *  @param v_x [input] x-component of the vector 
	 *  @param v_y [input] y-component of the vector 
	 *  @param v_p [output] \f$\phi\f$-component of the vector 
	 */
	virtual void comp_p_from_cartesian(const Cmp& v_x, const Cmp& v_y, 
					   Cmp& v_p) const ; 
	
	/**
	 * Decreases by 1 the value of \c dzpuis  of a \c Scalar 
	 *  and changes accordingly its values in the 
	 *  compactified external domain (CED).
	 */
	virtual void dec_dzpuis(Scalar& ) const ; 

	/**
	 * Decreases by 2 the value of \c dzpuis  of a \c Scalar 
	 *  and changes accordingly its values in the  
	 *  compactified external domain (CED).
	 */
	virtual void dec2_dzpuis(Scalar& ) const ; 

	/**
	 * Increases by 1 the value of \c dzpuis  of a \c Scalar 
	 *  and changes accordingly its values in the  
	 *  compactified external domain (CED).
	 */
	virtual void inc_dzpuis(Scalar& ) const ; 
	
	/**
	 * Increases by 2 the value of \c dzpuis  of a \c Scalar 
	 *  and changes accordingly its values in the  
	 *  compactified external domain (CED).
	 */
	virtual void inc2_dzpuis(Scalar& ) const ; 
	

    // PDE resolution :
    // --------------
    public:
	/** Resolution of the elliptic equation 
	 * \f$ a \Delta\psi + {\bf b} \cdot \nabla \psi = \sigma\f$.
	 * 
	 * @param source [input] source \f$\sigma\f$ of the above equation
	 * @param aa [input] factor \e a  in the above equation
	 * @param bb [input] vector \e \b b in the above equation
	 * @param par [input/output] parameters of the iterative method of
	 *  resolution : \\
	 *  \c par.get_int(0)  : [input] maximum number of iterations \\
	 *  \c par.get_double(0)  : [input] required precision: the iterative
	 *	  method is stopped as soon as the relative difference between
	 *	 \f$\psi^J\f$ and \f$\psi^{J-1}\f$ is greater than 
	 *	 \c par.get_double(0) \\
	 *  \c par.get_double(1)  : [input] relaxation parameter \f$\lambda\f$ \\
	 *  \c par.get_int_mod(0)  : [output] number of iterations 
	 *				    actually used to get the solution.
	 * @param psi [input/output]: input : previously computed value of \f$\psi\f$
	 *	to start the iteration (if nothing is known a priori, 
	 *	\c psi  must be set to zero); 
	 *	output: solution \f$\psi\f$ which satisfies \f$\psi(0)=0\f$.  
	 */ 
	virtual void poisson_compact(const Cmp& source, const Cmp& aa, 
				     const Tenseur& bb, const Param& par, 
				     Cmp& psi) const ;
};


			//------------------------------------//
			//            class Map_af            //
			//------------------------------------//



/**
 * Affine radial mapping. \ingroup (map)
 * 
 * The affine radial mapping is the simplest one between the grid coordinates
 * \f$(\xi, \theta', \phi')\f$ and the physical coordinates \f$(r, \theta, \phi)\f$. 
 * It is defined by \f$\theta=\theta'\f$, \f$\phi=\phi'\f$ and 
 *  \li \f$r=\alpha \xi + \beta\f$, in non-compactified domains, 
 *  \li \f$ u={1\over r} = \alpha \xi + \beta\f$ in the (usually outermost) compactified
 *        domain, 
 * where \f$\alpha\f$ and \f$\beta\f$ are constant in each domain. 
 * 
 * 
 */

class Map_af : public Map_radial {

    // Data :
    // ----
    private:
	/// Array (size: \c mg->nzone ) of the values of \f$\alpha\f$ in each domain
	double* alpha ;	 
	/// Array (size: \c mg->nzone ) of the values of \f$\beta\f$ in each domain
	double* beta ;	  
	
    // Constructors, destructor : 
    // ------------------------
    public:
	/**
	 * Standard Constructor
	 * @param mgrille  [input] Multi-domain grid on which the mapping is defined
	 * @param r_limits [input] Array (size: number of domains + 1) of the
	 *			   value of \e r  at the boundaries of the various 
	 *			   domains : 
	 *			   \li \c r_limits[l] : inner boundary of the 
	 *				 domain no. \c l 
	 *			   \li \c r_limits[l+1] : outer boundary of the 
	 *				 domain no. \c l 
	 */
	Map_af(const Mg3d& mgrille, const double* r_limits) ;	
	/**
	 * Standard Constructor with Tbl
	 * @param mgrille  [input] Multi-domain grid on which the mapping is defined
	 * @param r_limits [input] Array (size: number of domains) of the
	 *			   value of \e r  at the boundaries of the various 
	 *			   domains : 
	 *			   \li \c r_limits[l] : inner boundary of the 
	 *				 domain no. \c l 
	 *			   \li \c r_limits[l+1] : outer boundary of the 
	 *				 domain no. \c l except in the last domain
	 * The last boundary is set to inifnity if the grid contains a compactified domain.
	 */
	Map_af(const Mg3d& mgrille, const Tbl& r_limits) ;	
	
	Map_af(const Map_af& ) ;      ///< Copy constructor
	Map_af(const Mg3d&, FILE* ) ; ///< Constructor from a file (see \c sauve(FILE*) )
	
	/** Constructor from a general mapping.
	 * 
	 *  If the input mapping belongs to the class \c Map_af , this
	 *  constructor does the same job as the copy constructor 
	 *  \c Map_af(const \c Map_af& \c ) .
	 * 
	 *  If the input mapping belongs to the class \c Map_et , this
	 *  constructor sets in each domain, the values of 
	 *  \f$\alpha\f$ and \f$\beta\f$ to those of the \c Map_et .
	 *  
	 */
	explicit Map_af(const Map& ) ;      

	virtual ~Map_af() ;	      ///< Destructor

    // Assignment
    // ----------
    public: 
	/// Assignment to another affine mapping. 
	virtual void operator=(const Map_af& ) ;
        
    // Memory management
    // -----------------
    private:
	/// Assignment of the building functions to the member \c Coords
	void set_coord() ;	

    // Extraction of information
    // -------------------------
    public:
	/// Returns the pointer on the array \c alpha 
	const double* get_alpha() const ; 

	/// Returns the pointer on the array \c beta 
	const double* get_beta() const ; 
	
	/**
	 *  Returns the value of the radial coordinate \e r  for a given
	 *  \f$(\xi, \theta', \phi')\f$ in a given domain. 
	 *	@param l [input] index of the domain
	 *	@param xi [input] value of \f$\xi\f$
	 *	@param theta [input] value of \f$\theta'\f$
	 *	@param pphi [input] value of \f$\phi'\f$
	 *	@return value of \f$r=R_l(\xi, \theta', \phi')\f$
	 */
	virtual double val_r(int l, double xi, double theta, double pphi) const ; 

	/**
	 * Computes the domain index \e l  and the value of \f$\xi\f$ corresponding
	 * to a point given by its physical coordinates \f$(r, \theta, \phi)\f$. 
	 *	@param rr [input] value of \e r 
	 *	@param theta [input] value of \f$\theta\f$
	 *	@param pphi [input] value of \f$\phi\f$
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of \f$\xi\f$
	 */
	virtual void val_lx(double rr, double theta, double pphi,
			    int& l, double& xi) const ; 
	
	/** Computes the domain index \e l  and the value of \f$\xi\f$ corresponding
	 * to a point given by its physical coordinates \f$(r, \theta, \phi)\f$. 
	 *	@param rr [input] value of \e r 
	 *	@param theta [input] value of \f$\theta\f$
	 *	@param pphi [input] value of \f$\phi\f$
	 *	@param par [] unused by the \c Map_af  version
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of \f$\xi\f$
	 */
	virtual void val_lx(double rr, double theta, double pphi, 
			    const Param& par, int& l, double& xi) const ; 
		
	/** Returns the value of the radial coordinate \e r  for a given
	 *  \f$\xi\f$ and a given collocation point in \f$(\theta', \phi')\f$ 
	 *   in a given domain. 
	 *	@param l [input] index of the domain
	 *	@param xi [input] value of \f$\xi\f$
	 *	@param j [input] index of the collocation point in \f$\theta'\f$
	 *	@param k [input] index of the collocation point in \f$\phi'\f$
	 *	@return value of \f$r=R_l(\xi, {\theta'}_j, {\phi'}_k)\f$
	 */
	virtual double val_r_jk(int l, double xi, int j, int k) const ; 
	
	/** Computes the domain index \e l  and the value of \f$\xi\f$ corresponding
	 * to a point of arbitrary \e r  but collocation values of \f$(\theta, \phi)\f$ 
	 *	@param rr [input] value of \e r 
	 *	@param j [input] index of the collocation point in \f$\theta\f$
	 *	@param k [input] index of the collocation point in \f$\phi\f$
	 *	@param par [] unused by the \c Map_af  version
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of \f$\xi\f$
	 */
	virtual void val_lx_jk(double rr, int j, int k, const Param& par, 
			       int& l, double& xi) const ; 

	/// Comparison operator (egality)
	virtual bool operator==(const Map& ) const ;  


    // Outputs
    // -------
    public:
	virtual void sauve(FILE* ) const ;	  ///< Save in a file
    
    private:
	virtual ostream& operator>>(ostream &) const ;    ///< Operator >>

    // Modification of the mapping
    // ---------------------------
    public:
	/** Sets a new radial scale.
	 *	@param lambda [input] factor by which the value of \e r  is to
	 *	    be multiplied
	 */
	virtual void homothetie(double lambda) ;
	
	/** Rescales the outer boundary of one domain.
	 *  The inner boundary is unchanged. The inner boundary 
	 *  of the next domain is changed to match the new outer
	 *  boundary. 
	 *	@param l [input] index of the domain
	 *	@param lambda [input] factor by which the value of 
	 *	    \f$R(\theta, \varphi)\f$ defining the outer boundary
	 *	    of the domain is to be multiplied. 
	 */
	virtual void resize(int l, double lambda) ; 

	/** Sets a new radial scale at the bondary between the nucleus and the
	 * first shell.
	 *	@param lambda [input] factor by which the value of \e r  is to
	 *	    be multiplied
	 */
	void homothetie_interne(double lambda) ;
	
	/** Adaptation of the mapping to a given scalar field.
	 */
	virtual void adapt(const Cmp& ent, const Param& par) ; 

	/// Modifies the value of \f$\alpha\f$ in domain no. \e l 
	void set_alpha(double alpha0, int l) ;

	/// Modifies the value of \f$\beta\f$ in domain no. \e l 
	void set_beta(double beta0, int l) ;  

    // Differential operators:
    // ----------------------
    public:
	/** Computes \f$\partial/ \partial r\f$ of a \c Cmp.
	 *  Note that in the  compactified external domain (CED), it computes
	 *  \f$-\partial/ \partial u = r^2 \partial/ \partial r\f$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of \c ci 
	 */
	virtual void dsdr(const Cmp& ci, Cmp& resu) const ;  	    

	/** Computes \f$1/r \partial/ \partial \theta\f$ of a \c Cmp.
	 *  Note that in the  compactified external domain (CED), it computes
	 *  \f$1/u \partial/ \partial \theta = r \partial/ \partial \theta\f$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of \c ci 
	 */
	virtual void srdsdt(const Cmp& ci, Cmp& resu) const ;  	    
	
	/** Computes \f$1/(r\sin\theta) \partial/ \partial \phi\f$ of a \c Cmp.
	 *  Note that in the compactified external domain (CED), it computes
	 *  \f$1/(u\sin\theta) \partial/ \partial \phi = 
	 *	    r/\sin\theta \partial/ \partial \phi\f$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of \c ci 
	 */
	virtual void srstdsdp(const Cmp& ci, Cmp& resu) const ;  	    

	/** Computes \f$\partial/ \partial r\f$ of a \c Scalar.
	 *  Note that in the  compactified external domain (CED), the \c dzpuis 
	 *  flag of the output is 2 if the input has \c dzpuis  = 0, and 
	 *  is increased by 1 in other cases.
	 *  @param uu [input] field to consider
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void dsdr(const Scalar& uu, Scalar& resu) const  ;  
	    
	/** Computes \f$\partial/ \partial r\f$ of a \c Scalar.
	 *  Note that in the  compactified external domain (CED), the \c dzpuis 
	 *  flag of the output is 2 if the input has \c dzpuis  = 0, and 
	 *  is increased by 1 in other cases.
	 *  @param uu [input] field to consider
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void dsdradial(const Scalar&, Scalar&) const  ;  	    

	/** Computes \f$1/r \partial/ \partial \theta\f$ of a \c Scalar.
	 *  Note that in the  compactified external domain (CED), the \c dzpuis 
	 *  flag of the output is 2 if the input has \c dzpuis  = 0, and 
	 *  is increased by 1 in other cases.
	 *  @param uu [input] field to consider
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void srdsdt(const Scalar& uu, Scalar& resu) const  ;  	    
	
	/** Computes \f$1/(r\sin\theta) \partial/ \partial \phi\f$ of a \c Scalar.
	 *  Note that in the  compactified external domain (CED), the \c dzpuis 
	 *  flag of the output is 2 if the input has \c dzpuis  = 0, and 
	 *  is increased by 1 in other cases.
	 *  @param uu [input] field to consider
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void srstdsdp(const Scalar& uu, Scalar& resu) const ;  	    
        
	/** Computes \f$\partial/ \partial \theta\f$ of a \c Scalar.
	 *  @param uu [input] scalar field 
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void dsdt(const Scalar& uu, Scalar& resu) const ;  	    

	/** Computes \f$1/\sin\theta \partial/ \partial \varphi\f$ of a \c Scalar.
	 *  @param uu [input] scalar field 
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void stdsdp(const Scalar& uu, Scalar& resu) const ;  	    

	/** Computes the Laplacian of a scalar field.
	 *   @param uu	[input]  Scalar field \e u  (represented as a \c Cmp)
	 *			 the Laplacian \f$\Delta u\f$ of which is to be computed
	 *   @param zec_mult_r [input] Determines the quantity computed in
	 *			 the  compactified external domain (CED) :  \\
	 *		    zec_mult_r = 0 : \f$\Delta u\f$	\\
	 *		    zec_mult_r = 2 : \f$r^2 \,  \Delta u\f$	\\
	 *		    zec_mult_r = 4 (default) : \f$r^4 \, \Delta u\f$	
	 *  @param lap [output] Laplacian of \e u 
	 */
	virtual void laplacien(const Cmp& uu, int zec_mult_r, Cmp& lap) const ; 
	
	
	
	/** Computes the angular Laplacian of a scalar field.
	 *   @param uu	[input]  Scalar field \e u  (represented as a \c Scalar)
	 *			 the Laplacian \f$\Delta u\f$ of which is to be computed
	 *   @param lap [output] Angular Laplacian of \e u  (see documentation 
	 *                       of \c Scalar
	 */
	virtual void lapang(const Scalar& uu, Scalar& lap) const ; 
	

	/** Computes the radial primitive which vanishes for \f$r\to \infty\f$.
	 *  i.e. the function 
	 *   \f$ F(r,\theta,\varphi) = \int_r^\infty f(r',\theta,\varphi) \, dr' \f$
	 *
	 *      @param uu [input] function \e f (must have a \c dzpuis = 2)
	 *      @param resu [input] function \e F
	 *      @param null_infty if true (default), the primitive is null
	 *      at infinity (or on the grid boundary). \e F is null at the
	 *      center otherwise
	 */ 
	virtual void primr(const Scalar& uu, Scalar& resu,
			   bool null_infty) const ;  	    

	 
	/** Computes the integral over all space of a \c Cmp.
	 *  The computed quantity is 
	 *    \f$\int u \, r^2 \sin\theta \,  dr\, d\theta \, d\phi\f$.
	 *  The routine allocates a \c Tbl  (size: \c mg->nzone ) to store 
	 *  the result (partial integral) in each domain and returns a pointer 
	 *  to it.
	 */
	virtual Tbl* integrale(const Cmp&) const ; 
	
	 
    // PDE resolution :
    // --------------
    public:
	/** Computes the solution of a scalar Poisson equation.
	 *   @param source [input] source \f$\sigma\f$ of the Poisson equation 
	 *	    \f$\Delta u = \sigma\f$.
	 *   @param par [] not used by this \c Map_af  version. 
	 *   @param uu [output] solution \e u  with the boundary condition 
	 *	    \e u =0 at spatial infinity. 
	 */
	virtual void poisson(const Cmp& source, Param& par, Cmp& uu) const ;

	virtual void poisson_falloff(const Cmp& source, Param& par, Cmp& uu, int k_falloff) const ;

	/** Computes the solution of a scalar Poisson equation.
	 *   The regularized source
         *   \f$\sigma_{\rm regu} = \sigma - \sigma_{\rm div}\f$
         *   is constructed and solved.
	 *   @param source [input] source \f$\sigma\f$ of the Poisson equation 
	 *	    \f$\Delta u = \sigma\f$.
	 *   @param k_div [input] regularization degree of the procedure
	 *   @param nzet [input] number of domains covering the star
	 *   @param unsgam1 [input] parameter \f$1/(\gamma-1)\f$ where \f$\gamma\f$
         *          denotes the adiabatic index.
	 *   @param par [] not used by this \c Map_af  version.
	 *   @param uu [output] solution \e u  with the boundary condition 
	 *	    \e u =0 at spatial infinity.
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

	/** Computes the solution of an angular Poisson equation.
	 * The angular Poisson equation is \f$\Delta_{\theta\varphi} u = \sigma\f$,
	 * where \f$\Delta_{\theta\varphi} u := \frac{\partial^2 u}
	 *  {\partial \theta^2} + \frac{1}{\tan \theta} \frac{\partial u}
	 *  {\partial \theta} +\frac{1}{\sin^2 \theta}\frac{\partial^2 u}
	 *  {\partial \varphi^2}\f$.
	 * 
	 *   @param source [input] source \f$\sigma\f$ of the Poisson equation 
	 *	    \f$\Delta_{\theta\varphi} u = \sigma\f$.
	 *   @param par There must be one integer to force the result to
	 *          be expressed in spherical harmonics; otherwise (empty
	 *          \c Param ), standard angular bases are used. 
	 *   @param uu [input/output] solution \e u  
	 */
	virtual void poisson_angu(const Scalar& source, Param& par, 
					Scalar& uu) const ;

	/**
	 * Internal function intended to be used by \c Map::poisson_vect 
	 * and \c Map::poisson_vect_oohara . It constructs the sets of 
	 * parameters used for each scalar Poisson equation from the one for 
	 * the vectorial one.
	 * 
	 * In the case of a \c Map_af  the result is not used and the function 
	 * only returns \c & \c par .
	 */
	virtual Param* donne_para_poisson_vect (Param& par, int i) const ;
	
	/**
	 * Solver of the Poisson equation with boundary condition for the 
	 * affine mapping case.
	 */
	virtual void poisson_frontiere (const Cmp&, const Valeur&, int, int, Cmp&, double = 0., double = 0.) const ;

	/**
	 * Solver of the Poisson equation with boundary condition for the 
	 * affine mapping case, cases with boundary conditions of both 
	 * Dirichlet and Neumann type (no condition imposed at infinity).
	 */
	virtual void poisson_frontiere_double (const Cmp& source, const Valeur& lim_func,
			const Valeur& lim_der, int num_zone, Cmp& pot) const  ;

	/**
	 * Computes the solution of a Poisson equation in the shell,
	 * imposing a boundary condition at the surface of the star 
	 * 
	 * @param source [input] : source of the equation.
	 * @param limite [input] : \c limite[num_front]  contains the angular 
	 * function being the boudary condition.
	 * @param par [input] : parameters of the computation.
	 * @param pot [output] : result.
	 */
	virtual void poisson_interne (const Cmp& source, const Valeur& limite,
			 Param& par, Cmp& pot) const ;

	/**
	 * Performs the surface integration of \c ci  on the sphere of 
	 * radius \c rayon .
	 */
	double integrale_surface (const Cmp& ci, double rayon) const ;
	
	/**
	 * Performs the surface integration of \c ci  on the sphere of 
	 * radius \c rayon .
	 */
	double integrale_surface (const Scalar& ci, double rayon) const ;
	
	double integrale_surface_falloff (const Cmp& ci) const ;

	/**
	 * Performs the surface integration of \c ci  at infinity.
	 * \c ci  must have \c dzpuis  =2.
	 */
	double integrale_surface_infini (const Cmp& ci) const ;
	
	/**
	 * Performs the surface integration of \c ci  at infinity.
	 * \c ci  must have \c dzpuis  =2.
	 */
	double integrale_surface_infini (const Scalar& ci) const ;
	
	/**
	 * General elliptic solver. The field is zero at infinity.
	 *
	 * @param params [input] : the operators and variables to be uses.
	 * @param so [input] : the source.
	 * @param uu [output] : the solution.
	 **/
	void sol_elliptic (Param_elliptic& params, 
			   const Scalar& so,  Scalar& uu) const ;

	/**
	 * General elliptic solver.
	 * The equation is not solved in the compactified domain.
	 *
	 * @param params [input] : the operators and variables to be uses.
	 * @param so [input] : the source.
	 * @param uu [output] : the solution.
	 * @param val [input] : value at the last shell.
	 **/
	void sol_elliptic_no_zec (Param_elliptic& params, 
				  const Scalar& so, Scalar& uu, double val) const ;

	/**
	 * General elliptic solver.
	 * The equation is solved only in the compactified domain.
	 *
	 * @param params [input] : the operators and variables to be uses.
	 * @param so [input] : the source.
	 * @param uu [output] : the solution.
	 * @param val [input] : value at the inner boundary.
	 **/
	void sol_elliptic_only_zec (Param_elliptic& params, 
				  const Scalar& so, Scalar& uu, double val) const ;

	/**
	 * General elliptic solver.
	 * The equation is not solved in the compactified domain and the 
	 * matching is done with a solution of the type 
	 * \f$\frac{\sin \left( f r + \phi \right)}{r}\f$, minimizing the coefficient 
	 * with respect to \f$\phi\f$.
	 *
	 * @param params [input] : the operators and variables to be uses.
	 * @param so [input] : the source.
	 * @param uu [output] : the solution.
	 * @param freq [input] : the frequency \f$f\f$.
	 * @param nbr_phase [input] : number of points for the maximization 
	 * over \f$\phi\f$.
	 * @param error [output] : coefficient of the oscillatory solution 
	 * in the external domain.
	 * @param phase [output] : phase that minimizes \c coef .
	 **/
	void sol_elliptic_sin_zec (Param_elliptic& params, 
				  const Scalar& so, Scalar& uu, double freq, 
				   int nbr_phase, double& coef, 
				   double& phase) const ;
	/**
	 * General elliptic solver fixing the derivative at the origin 
	 * and relaxing the continuity of the first derivative at the 
	 * boundary between the nucleus and the first shell.
	 *
	 * @param val [input] : valeur of the derivative.
	 * @param params [input] : the operators and variables to be uses.
	 * @param so [input] : the source.
	 * @param uu [output] : the solution.
	 **/
	void sol_elliptic_fixe_der_zero (double val, 
					 Param_elliptic& params, 
					 const Scalar& so, Scalar& uu) const ;
	
	/** Computes the solution of a 2-D Poisson equation.
	 *  The 2-D Poisson equation writes
	 *  \f[
	 *	{\partial^2 u\over\partial r^2} + 
	 *	    {1\over r} {\partial u \over \partial r} + 
	 *	    {1\over r^2} {\partial^2 u\over\partial \theta^2} = 
	 *		\sigma \ . 
	 *  \f] 
	 *
	 *   @param source_mat [input] Compactly supported part of 
	 *	    the source \f$\sigma\f$ of the 2-D Poisson equation (typically
	 *	    matter terms)
	 *   @param source_quad [input] Non-compactly supported part of 
	 *	    the source \f$\sigma\f$ of the 2-D Poisson equation (typically
	 *	    quadratic terms)
	 *   @param par [output] Parameter which contains the constant
	 *			    \f$\lambda\f$ used to fulfill the virial
	 *			 identity GRV2 : \\ 
	 *  \c par.get_double_mod(0)  : [output] constant \c lambda 
	 *	    such that the source of the equation effectively solved
	 *	    is \c source_mat + \c lambda * \c source_quad . 
	 *   @param uu [input/output] solution \e u  with the boundary condition 
	 *	    \e u =0 at spatial infinity. 
	 */
	virtual void poisson2d(const Cmp& source_mat, const Cmp& source_quad, 
			       Param& par, Cmp& uu) const ;
	/**
	 * General elliptic solver in a 2D case. The field is zero at infinity.
	 *
	 * @param params [input] : the operators and variables to be uses.
	 * @param so [input] : the source.
	 * @param uu [output] : the solution.
	 **/
	void sol_elliptic_2d(Param_elliptic&, 
			     const Scalar&, Scalar&) const ;
	/**
	 * General elliptic solver in a pseudo 1d case. The field is zero at infinity.
	 *
	 * @param params [input] : the operators and variables to be uses.
	 * @param so [input] : the source.
	 * @param uu [output] : the solution.
	 **/
	void sol_elliptic_pseudo_1d(Param_elliptic&, 
			     const Scalar&, Scalar&) const ;

	/** Performs one time-step integration of the d'Alembert scalar equation
	 *   @param par [input/output] possible parameters to control the
	 *   resolution of the d'Alembert equation: \\
	 *   \c par.get_double(0)  : [input] the time step \e dt ,\\
	 *   \c par.get_int(0)  : [input] the type of boundary conditions
	 *   set at the outer boundary (0 : reflexion, 1 : Sommerfeld 
	 *   outgoing wave, valid only for \e l=0  components, 2 : Bayliss 
	 *   & Turkel outgoing wave, valid for \e l=0, 1, 2  components)\\
	 *   \c par.get_int_mod(0)  : [input/output] set to 0 at first
	 *   call, is used as a working flag after (must not be modified after
	 *   first call)\\
	 *   \c par.get_Cmp_mod(0)  : [input] (optional) if the wave 
	 *   equation is on a curved space-time, this is the potential in front
	 *   of the Laplace operator. It has to be updated at every time-step
	 *   (for a potential depending on time).\\
	 *   Note: there are many other working objects attached to this
	 *   \c Param , so one should not modify it.
	 *   @param fJp1 [output] solution \f$f^{J+1}\f$ at time \e J+1 
	 *   with boundary conditions defined by \c par.get_int(0) 
	 *   @param fJ [input] solution \f$f^J\f$ at time \e J 
	 *   @param fJm1 [input] solution \f$f^{J-1}\f$ at time \e J-1 
	 *   @param source [input] source \f$\sigma\f$ of the d'Alembert equation 
	 *	    \f$\diamond u = \sigma\f$.
	 */
	virtual void dalembert(Param& par, Scalar& fJp1, const Scalar& fJ, 
			       const Scalar& fJm1, const Scalar& source) const ;

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
    friend Mtbl* map_af_fait_drdt(const Map* ) ;
    friend Mtbl* map_af_fait_stdrdp(const Map* ) ;
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
     Mtbl* map_af_fait_drdt(const Map* ) ;
     Mtbl* map_af_fait_stdrdp(const Map* ) ;
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
			//            class Map_et            //
			//------------------------------------//



/**
 * Radial mapping of rather general form. \ingroup (map)
 *
 * This mapping relates the grid coordinates
 * \f$(\xi, \theta', \phi')\f$ and the physical coordinates \f$(r, \theta, \phi)\f$
 * in the following manner [see Bonazzola, Gourgoulhon & Marck,
 * \e Phys. \e Rev. \e D  \b 58 , 104020 (1998) for details]:
 * \f$\theta=\theta'\f$, \f$\phi=\phi'\f$ and
 *  \li \f$r = \alpha [\xi + A(\xi) F(\theta', \phi') + B(\xi) G(\theta', \phi')]
 *	       + \beta\f$ in non-compactified domains,
 *  \li \f$ u={1\over r} = \alpha [\xi + A(\xi) F(\theta', \phi')] + \beta\f$ in
 *	    the (usually outermost) compactified domain,
 * where \f$\alpha\f$ and \f$\beta\f$ are constant in each domain, \f$A(\xi)\f$ and
 * \f$B(\xi)\f$ are constant polynomials defined by
 *	\li \f$A(\xi) = 3 \xi^4 - 2 \xi^6\f$ and 
 *	      \f$B(\xi) = (5\xi^3 - 3\xi^5)/2\f$ in the nucleus (innermost domain
 *		    which contains \e r =0);
 *	\li \f$A(\xi) = (\xi^3 - 3\xi + 2)/4\f$ and
 *	      \f$B(\xi) = (-\xi^3 + 3\xi +2)/4\f$ in the other domains.
 * The functions \f$F(\theta', \phi')\f$ and \f$G(\theta', \phi')\f$ depend on the
 * domain under consideration and define the boundaries of this domain. 	
 *    
 * 
 */

class Map_et : public Map_radial {

    // Data :
    // ----
    private:
	/// Array (size: \c mg->nzone ) of the values of \f$\alpha\f$ in each domain
	double* alpha ;	  
	/// Array (size: \c mg->nzone ) of the values of \f$\beta\f$ in each domain
	double* beta ;	  

	/** Array (size: \c mg->nzone ) of \c Tbl  which stores the 
	 *  values of \f$A(\xi)\f$ in each domain
	 */
	Tbl** aa ;	  

	/** Array (size: \c mg->nzone ) of \c Tbl  which stores the 
	 *  values of \f$A'(\xi)\f$ in each domain
	 */
	Tbl** daa ;	  
	
	/** Array (size: \c mg->nzone ) of \c Tbl  which stores the 
	 *  values of \f$A''(\xi)\f$ in each domain
	 */
	Tbl** ddaa ;	  
	
	/// Values at the \c nr  collocation points of \f$A(\xi)/\xi\f$ in the nucleus
	Tbl aasx ; 
	 
	/// Values at the \c nr  collocation points of \f$A(\xi)/\xi^2\f$ in the nucleus
	Tbl aasx2 ; 
	 
	/** Values at the \c nr  collocation points of \f$A(\xi)/(\xi-1)\f$
	 *  in the outermost compactified domain
	 */
	Tbl zaasx ; 
	  
	/** Values at the \c nr  collocation points of \f$A(\xi)/(\xi-1)^2\f$
	 *  in the outermost compactified domain
	 */
	Tbl zaasx2 ; 
	  
	/** Array (size: \c mg->nzone ) of \c Tbl  which stores the 
	 *  values of \f$B(\xi)\f$ in each domain
	 */
	Tbl** bb ;	  

	/** Array (size: \c mg->nzone ) of \c Tbl  which stores the 
	 *  values of \f$B'(\xi)\f$ in each domain
	 */
	Tbl** dbb ;	  
	
	/** Array (size: \c mg->nzone ) of \c Tbl  which stores the 
	 *  values of \f$B''(\xi)\f$ in each domain
	 */
	Tbl** ddbb ;	  
	
	/// Values at the \c nr  collocation points of \f$B(\xi)/\xi\f$ in the nucleus
	Tbl bbsx ; 
	 
	/// Values at the \c nr  collocation points of \f$B(\xi)/\xi^2\f$ in the nucleus
	Tbl bbsx2 ; 
	 
	/** Values of the function \f$F(\theta', \phi')\f$ at the \c nt*np 
	 * angular collocation points in each domain.
	 * The \c Valeur \c ff is defined on the multi-grid \c mg->g_angu 
	 * (cf. class \c Mg3d ). 
	 */
	Valeur ff ; 
	  
	/** Values of the function \f$G(\theta', \phi')\f$ at the \c nt*np 
	 *   angular collocation points in each domain.
	 * The \c Valeur \c gg is defined on the multi-grid \c mg->g_angu 
	 * (cf. class \c Mg3d ). 
	 */
	Valeur gg ; 
	  
    public:
	/** \f$1/(\partial R/\partial \xi) R/\xi\f$ in the nucleus; \\
	 *  \f$1/(\partial R/\partial \xi) R/(\xi + \beta/\alpha)\f$ in the shells; \\
	 *  \f$1/(\partial U/\partial \xi) U/(\xi-1)\f$ in the outermost 
	 *  compactified domain.
	 */
	Coord rsxdxdr ;
	 
	/** \f$[ R/ (\alpha \xi + \beta) ]^2 (\partial R/\partial \xi) / \alpha\f$
	 *  in the nucleus and the shells; \\
	 *  \f$\partial U/\partial \xi / \alpha\f$ in the outermost compactified
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
	 *			   value of \e r  at the boundaries of the various 
	 *			   domains : 
	 *			   \li \c r_limits[l] : inner boundary of the 
	 *				 domain no. \c l 
	 *			   \li \c r_limits[l+1] : outer boundary of the 
	 *				 domain no. \c l 
	 */
	Map_et(const Mg3d& mgrille, const double* r_limits) ;	
	
	/**
	 * Constructor using the equation of the surface of the star.
	 * @param mgrille [input] Multi-domain grid on which the mapping is defined
	 * It must contains at least one shell.
	 * @param r_limits [input] Array (size: number of domains + 1) of the
	 *			   value of \e r  at the boundaries of the various 
	 *			   domains : 
	 *			   \li \c r_limits[l] : inner boundary of the 
	 *				 domain no. \c l 
	 *			   \li \c r_limits[l+1] : outer boundary of the 
	 *				 domain no. \c l 
	 * The first value is not used.
	 * @param tab [input] equation of the surface between the nucleus and the first
	 * shell in the form : \f${\rm tab}(k,j) = r(\phi_k,\theta_j)\f$, where
	 * \f$\phi_k\f$ and \f$\theta_j\f$ are the values of the angular colocation points.
	 */

	Map_et(const Mg3d& mgrille, const double* r_limits,const Tbl& tab);	
	Map_et(const Map_et& ) ;      ///< Copy constructor
	Map_et(const Mg3d&, FILE* ) ; ///< Constructor from a file (see \c sauve(FILE*) )

	virtual ~Map_et() ;	      ///< Destructor

    // Assignment
    // ----------
    public:
	/// Assignment to another \c Map_et  
	virtual void operator=(const Map_et& mp) ;
	
	/// Assignment to an affine mapping. 
	virtual void operator=(const Map_af& mpa) ;
	
	/// Assigns a given value to the function \f$F(\theta',\phi')\f$ 
	void set_ff(const Valeur& ) ; 
	/// Assigns a given value to the function \f$G(\theta',\phi')\f$ 
	void set_gg(const Valeur& ) ; 

    // Memory management
    // -----------------
    private:
	/// Assignement of the building functions to the member \c Coords
	void set_coord() ;		
    protected:
	/// Resets all the member \c Coords	
	virtual void reset_coord() ;
	
    private: 
	/// Construction of the polynomials \f$A(\xi)\f$ and \f$B(\xi)\f$
	void fait_poly() ; 	
	
    // Extraction of information
    // -------------------------
    public:
	/** Returns a pointer on the array \c alpha  (values of \f$\alpha\f$
	 *   in each domain)
	 */
	const double* get_alpha() const ; 

	/** Returns a pointer on the array \c beta  (values of \f$\beta\f$
	 *   in each domain)
	 */
	const double* get_beta() const ; 

	/// Returns a (constant) reference to the function \f$F(\theta',\phi')\f$
	const Valeur& get_ff() const ; 

	/// Returns a (constant) reference to the function \f$G(\theta',\phi')\f$
	const Valeur& get_gg() const ; 
	
	/**
	 *  Returns the value of the radial coordinate \e r  for a given
	 *  \f$(\xi, \theta', \phi')\f$ in a given domain. 
	 *	@param l [input] index of the domain
	 *	@param xi [input] value of \f$\xi\f$
	 *	@param theta [input] value of \f$\theta'\f$
	 *	@param pphi [input] value of \f$\phi'\f$
	 *	@return value of \f$r=R_l(\xi, \theta', \phi')\f$
	 */
	virtual double val_r(int l, double xi, double theta, double pphi) const ; 

	/**
	 * Computes the domain index \e l  and the value of \f$\xi\f$ corresponding
	 * to a point given by its physical coordinates \f$(r, \theta, \phi)\f$. 
	 *	@param rr [input] value of \e r 
	 *	@param theta [input] value of \f$\theta\f$
	 *	@param pphi [input] value of \f$\phi\f$
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of \f$\xi\f$
	 */
	virtual void val_lx(double rr, double theta, double pphi,
			    int& l, double& xi) const ; 
	
	/** Computes the domain index \e l  and the value of \f$\xi\f$ corresponding
	 * to a point given by its physical coordinates \f$(r, \theta, \phi)\f$. 
	 * This version enables to pass some parameters to control the
	 * accuracy of the computation. 
	 *	@param rr [input] value of \e r 
	 *	@param theta [input] value of \f$\theta\f$
	 *	@param pphi [input] value of \f$\phi\f$
	 *	@param par [input/output] parameters to control the
	 *	    accuracy of the computation:  \\
	 *  \c par.get_int(0)  : [input] maximum number of iterations in the 
	 *		    secant method to locate \f$\xi\f$ \\
	 *  \c par.get_int_mod(0)  : [output] effective number of iterations 
	 *				    used \\
	 *  \c par.get_double(0)  : [input] absolute precision in the secant 
	 *				method to locate \f$\xi\f$ 
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of \f$\xi\f$
	 */
	virtual void val_lx(double rr, double theta, double pphi, 
			    const Param& par, int& l, double& xi) const ; 
		
	/// Comparison operator (egality)
	virtual bool operator==(const Map& ) const ;  

	/** Returns the value of the radial coordinate \e r  for a given
	 *  \f$\xi\f$ and a given collocation point in \f$(\theta', \phi')\f$ 
	 *   in a given domain. 
	 *	@param l [input] index of the domain
	 *	@param xi [input] value of \f$\xi\f$
	 *	@param j [input] index of the collocation point in \f$\theta'\f$
	 *	@param k [input] index of the collocation point in \f$\phi'\f$
	 *	@return value of \f$r=R_l(\xi, {\theta'}_j, {\phi'}_k)\f$
	 */
	virtual double val_r_jk(int l, double xi, int j, int k) const ; 
	
	/** Computes the domain index \e l  and the value of \f$\xi\f$ corresponding
	 * to a point of arbitrary \e r  but collocation values of \f$(\theta, \phi)\f$ 
	 *	@param rr [input] value of \e r 
	 *	@param j [input] index of the collocation point in \f$\theta\f$
	 *	@param k [input] index of the collocation point in \f$\phi\f$
	 *	@param par [input/output] parameters to control the
	 *	    accuracy of the computation:  \\
	 *  \c par.get_int(0)  : [input] maximum number of iterations in the 
	 *		    secant method to locate \f$\xi\f$ \\
	 *  \c par.get_int_mod(0)  : [output] effective number of iterations 
	 *				    used \\
	 *  \c par.get_double(0)  : [input] absolute precision in the secant 
	 *				method to locate \f$\xi\f$ 
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of \f$\xi\f$
	 */
	virtual void val_lx_jk(double rr, int j, int k, const Param& par, 
			       int& l, double& xi) const ; 



    // Outputs
    // -------
    public:
	virtual void sauve(FILE* ) const ;	  ///< Save in a file
    
    private:
	virtual ostream& operator>>(ostream &) const ;    ///< Operator >>

    // Modification of the radial scale
    // --------------------------------
    public:
	/** Sets a new radial scale.
	 *	@param lambda [input] factor by which the value of \e r  is to
	 *	    be multiplied
	 */
	virtual void homothetie(double lambda) ;	

	/** Rescales the outer boundary of one domain.
	 *  The inner boundary is unchanged. The inner boundary 
	 *  of the next domain is changed to match the new outer
	 *  boundary. 
	 *	@param l [input] index of the domain
	 *	@param lambda [input] factor by which the value of 
	 *	    \f$R(\theta, \varphi)\f$ defining the outer boundary
	 *	    of the domain is to be multiplied. 
	 */
	virtual void resize(int l, double lambda) ; 

	/** Rescales the outer boundary of the outermost domain
	 *  in the case of non-compactified external domain.
	 *  The inner boundary is unchanged.
	 *  @param lambda [input] factor by which the value of the radius
	 *                        of the outermost domain is to be multiplied.
	 */
	void resize_extr(double lambda) ;

    // Modification of the mapping
    // ---------------------------
	/** Adaptation of the mapping to a given scalar field.
	 *  Computes the functions \f$F(\theta',\phi')\f$ and \f$G(\theta',\phi')\f$
	 *  as well as the factors \f$\alpha\f$ and \f$\beta\f$, so that the 
	 *  boundaries of some domains coincide with isosurfaces of the
	 *  scalar field \c ent .
	 *  @param ent [input] scalar field, the isosurfaces of which are
	 *     used to determine the mapping
	 *  @param par [input/output] parameters of the computation: \\
	 *   \c par.get_int(0)  : maximum number of iterations to locate
	 *     zeros by the secant method \\
	 *   \c par.get_int(1)  : number of domains where the adjustment
	 *     to the isosurfaces of \c ent  is to be performed \\
	 *   \c par.get_int(2)  : number of domains \c nz_search  where 
	 *      the isosurfaces will be searched : the routine scans the 
	 *	\c nz_search  innermost domains, starting from the domain
	 *	of index \c nz_search-1 . NB: the field \c ent  must
	 *	be continuous over these domains \\
	 *   \c par.get_int(3)  : 1 = performs the full computation, 
	 *			     0 = performs only the rescaling 
	 *				by the factor 
	 *				\c par.get_double_mod(0)  \\
	 *   \c par.get_int(4)  : theta index of the collocation point
	 *			     \f$(\theta_*, \phi_*)\f$ [using the notations
	 *    of Bonazzola, Gourgoulhon & Marck, \e Phys. \e Rev. \e D  \b 58 , 
	 *    104020 (1998)] defining an isosurface of \c ent  \\
	 *   \c par.get_int(5)  : phi index of the collocation point
	 *			     \f$(\theta_*, \phi_*)\f$ [using the notations
	 *    of Bonazzola, Gourgoulhon & Marck, \e Phys. \e Rev. \e D  \b 58 , 
	 *    104020 (1998)] defining an isosurface of \c ent  \\
	 *   \c par.get_int_mod(0)  [output] : number of iterations 
	 *	actually used in the secant method \\
	 *   \c par.get_double(0)  : required absolute precision in the
	 *	 determination of zeros by the secant method \\
	 *   \c par.get_double(1)  : factor by which the values of \f$\lambda\f$
	 *	    and \f$\mu\f$ [using the notations
	 *    of Bonazzola, Gourgoulhon & Marck, \e Phys. \e Rev. \e D  \b 58 , 
	 *    104020 (1998)] will be multiplied : 1 = regular mapping, 
	 *	    0 = contracting mapping \\
	 *   \c par.get_double(2)  : factor by which all the radial distances
	 *				will be multiplied \\
	 *   \c par.get_tbl(0)  : array of values of the field \c ent  to
	 *			     define the isosurfaces. 
	 *  
	 */
	virtual void adapt(const Cmp& ent, const Param& par)  ; 

    // Differential operators:
    // ----------------------
    public:
	/** Computes \f$\partial/ \partial r\f$ of a \c Cmp.
	 *  Note that in the  compactified external domain (CED), it computes
	 *  \f$-\partial/ \partial u = r^2 \partial/ \partial r\f$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of \c ci 
	 */
	virtual void dsdr(const Cmp& ci, Cmp& resu) const ;  	    

	/** Computes \f$1/r \partial/ \partial \theta\f$ of a \c Cmp.
	 *  Note that in the compactified external domain (CED), it computes
	 *  \f$1/u \partial/ \partial \theta = r \partial/ \partial \theta\f$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of \c ci 
	 */
	virtual void srdsdt(const Cmp& ci, Cmp& resu) const ;  	    
	
	/** Computes \f$1/(r\sin\theta) \partial/ \partial \phi\f$ of a \c Cmp.
	 *  Note that in the  compactified external domain (CED), it computes
	 *  \f$1/(u\sin\theta) \partial/ \partial \phi = 
	 *	    r/\sin\theta \partial/ \partial \phi\f$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of \c ci 
	 */
	virtual void srstdsdp(const Cmp& ci, Cmp& resu) const ;  	        

	/** Computes \f$\partial/ \partial r\f$ of a \c Scalar.
	 *  Note that in the  compactified external domain (CED), the \c dzpuis 
	 *  flag of the output is 2 if the input has \c dzpuis  = 0, and 
	 *  is increased by 1 in other cases.
	 *  @param uu [input] field to consider
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void dsdr(const Scalar& uu, Scalar& resu) const ;
	
	/** Computes \f$\partial/ \partial r\f$ of a \c Scalar if the description is affine and 
	 * \f$\partial/ \partial \ln r\f$ if it is logarithmic.
	 *  Note that in the  compactified external domain (CED), the \c dzpuis 
	 *  flag of the output is 2 if the input has \c dzpuis  = 0, and 
	 *  is increased by 1 in other cases.
	 *  @param uu [input] field to consider
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void dsdradial(const Scalar& uu, Scalar& resu) const ;  	    
	
	/** Computes \f$1/r \partial/ \partial \theta\f$ of a \c Scalar.
	 *  Note that in the  compactified external domain (CED), the \c dzpuis 
	 *  flag of the output is 2 if the input has \c dzpuis  = 0, and 
	 *  is increased by 1 in other cases.
	 *  @param uu [input] field to consider
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void srdsdt(const Scalar& uu, Scalar& resu) const ;  	    
	
	/** Computes \f$1/(r\sin\theta) \partial/ \partial \phi\f$ of a \c Scalar.
	 *  Note that in the  compactified external domain (CED), the \c dzpuis 
	 *  flag of the output is 2 if the input has \c dzpuis  = 0, and 
	 *  is increased by 1 in other cases.
	 *  @param uu [input] field to consider
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void srstdsdp(const Scalar& uu, Scalar& resu) const ;  	    
    
	/** Computes \f$\partial/ \partial \theta\f$ of a \c Scalar.
	 *  @param uu [input] scalar field 
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void dsdt(const Scalar& uu, Scalar& resu) const ;  	    

	/** Computes \f$1/\sin\theta \partial/ \partial \varphi\f$ of a \c Scalar.
	 *  @param uu [input] scalar field 
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void stdsdp(const Scalar& uu, Scalar& resu) const ;  	    

	/** Computes the Laplacian of a scalar field.
	 *   @param uu	[input]  Scalar field \e u  (represented as a \c Cmp)
	 *			 the Laplacian \f$\Delta u\f$ of which is to be computed
	 *   @param zec_mult_r [input] Determines the quantity computed in
	 *			 the  compactified external domain (CED) :  \\
	 *		    zec_mult_r = 0 : \f$\Delta u\f$	\\
	 *		    zec_mult_r = 2 : \f$r^2 \,  \Delta u\f$	\\
	 *		    zec_mult_r = 4 (default) : \f$r^4 \, \Delta u\f$	
	 *  @param lap [output] Laplacian of \e u 
	 */
	virtual void laplacien(const Cmp& uu, int zec_mult_r, Cmp& lap) const ; 
	
	
	/** Computes the angular Laplacian of a scalar field.
	 *   @param uu	[input]  Scalar field \e u  (represented as a \c Scalar)
	 *			 the Laplacian \f$\Delta u\f$ of which is to be computed
	 *   @param lap [output] Angular Laplacian of \e u  (see documentation 
	 *                       of \c Scalar
	 */
	virtual void lapang(const Scalar& uu, Scalar& lap) const ; 
	

	/** Computes the radial primitive which vanishes for \f$r\to \infty\f$.
	 *  i.e. the function 
	 *   \f$ F(r,\theta,\varphi) = \int_r^\infty f(r',\theta,\varphi) \, dr' \f$
	 *
	 *      @param uu [input] function \e f (must have a \c dzpuis = 2)
	 *      @param resu [input] function \e F
	 *      @param null_infty if true (default), the primitive is null
	 *      at infinity (or on the grid boundary). \e F is null at the
	 *      center otherwise
	 */ 
	virtual void primr(const Scalar& uu, Scalar& resu,
			   bool null_infty) const ;  	    

	 
	/** Computes the integral over all space of a \c Cmp.
	 *  The computed quantity is 
	 *    \f$\int u \, r^2 \sin\theta \,  dr\, d\theta \, d\phi\f$.
	 *  The routine allocates a \c Tbl  (size: \c mg->nzone ) to store 
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
	 * Gourgoulhon & Marck, \e Phys. \e Rev. \e D  \b 58 , 104020 (1998),  
	 * the Poisson equation \f$\Delta u = \sigma\f$ is re-written
	 * as \f$a \tilde\Delta u = \sigma + R(u)\f$,  where \f$\tilde\Delta\f$
	 * is the Laplacian in an affine mapping and \e R(u)  contains the
	 * terms generated by the deviation of the mapping \c *this
	 * from spherical symmetry. This equation is solved by iterations.
	 * At each step \e J  the equation effectively solved is 
	 *  \f$\tilde\Delta u^{J+1} = s^J\f$ where
	 * \f[
	 *   s^J = 1/a_l^{\rm max} \{ {\tt source} + R(u^J) + (a_l^{\rm max}-a)
	 *          [ \lambda s^{J-1} + (1-\lambda) s^{J-2} ] \} \ ,  
	 * \f]
	 * with \f$a_l^{\rm max} := \max(a)\f$ in domain no. \e l  and \f$\lambda\f$
	 * is a relaxation parameter. 
	 *  @param source [input] source \f$\sigma\f$ of the Poisson equation
	 *  @param par [input/output] parameters for the iterative method: \\
	 *  \c par.get_int(0)  : [input] maximum number of iterations \\
	 *  \c par.get_double(0)  : [input] relaxation parameter \f$\lambda\f$ \\
	 *  \c par.get_double(1)  : [input] required precision: the iterative
	 *	  method is stopped as soon as the relative difference between
	 *	 \f$u^J\f$ and \f$u^{J-1}\f$ is greater than \c par.get_double(1) \\
	 *  \c par.get_cmp_mod(0)  : [input/output] input : \c Cmp 
	 *		containing \f$s^{J-1}\f$ (cf. the above equation) to 
	 *		start the iteration (if nothing is known a priori, 
	 *		this \c Cmp must be set to zero); output: value
	 *		of \f$s^{J-1}\f$, to used in a next call to the routine \\
	 *  \c par.get_int_mod(0)  : [output] number of iterations 
	 *				    actually used to get the solution.
	 *
	 *  @param uu [input/output] input : previously computed value of \e u 
	 *	to start the iteration (term \e R(u) ) (if nothing is known a 
	 *	priori, \c uu  must be set to zero); output: solution \e u  
	 *	with the boundary condition \e u =0 at spatial infinity. 
	 */
	virtual void poisson(const Cmp& source, Param& par, Cmp& uu) const ;

	virtual void poisson_falloff(const Cmp& source, Param& par, Cmp& uu, int k_falloff) const ;

	/** Computes the solution of a scalar Poisson equation.
	 *   The regularized source
	 *   @param source [input] source \f$\sigma\f$ of the Poisson equation 
	 *	    \f$\Delta u = \sigma\f$.
	 *   @param k_div [input] regularization degree of the procedure
	 *   @param nzet [input] number of domains covering the star
	 *   @param unsgam1 [input] parameter \f$1/(\gamma-1)\f$ where \f$\gamma\f$
         *          denotes the adiabatic index.
	 *   @param par [input/output] parameters for the iterative method: \\
	 *  \c par.get_int(0)  : [input] maximum number of iterations \\
	 *  \c par.get_double(0)  : [input] relaxation parameter
	 *                             \f$\lambda\f$ \\
	 *  \c par.get_double(1)  : [input] required precision: the
	 *        iterative method is stopped as soon as the relative
	 *        difference between \f$u^J\f$ and \f$u^{J-1}\f$ is greater than
	 *        \c par.get_double(1) \\
	 *  \c par.get_cmp_mod(0)  : [input/output] input : \c Cmp 
	 *		containing \f$s^{J-1}\f$ (cf. the above equation) to 
	 *		start the iteration (if nothing is known a priori, 
	 *		this \c Cmp must be set to zero); output: value
	 *		of \f$s^{J-1}\f$, to used in a next call to the routine \\
	 *  \c par.get_int_mod(0)  : [output] number of iterations 
	 *				    actually used to get the solution.
	 *   @param uu [input/output] input : previously computed value of \e u 
	 *	to start the iteration (term \e R(u) ) (if nothing is known a 
	 *	priori, \c uu  must be set to zero); output: solution \e u  
	 *	with the boundary condition \e u =0 at spatial infinity.
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

	/** Computes the solution of an angular Poisson equation.
	 * The angular Poisson equation is \f$\Delta_{\theta\varphi}  u = \sigma\f$,
	 * where \f$\Delta_{\theta\varphi} u := \frac{\partial^2 u}
	 *  {\partial \theta^2} + \frac{1}{\tan \theta} \frac{\partial u}
	 *  {\partial \theta} +\frac{1}{\sin^2 \theta}\frac{\partial^2 u}
	 *  {\partial \varphi^2}\f$.
	 * 
	 *   @param source [input] source \f$\sigma\f$ of the Poisson equation 
	 *	    \f$\Delta_{\theta\varphi} u = \sigma\f$.
	 *   @param par [input/output] parameters to control the
	 *   resolution of the Poisson equation.  
	 *   @param uu [input/output] solution \e u  
	 */
	virtual void poisson_angu(const Scalar& source, Param& par, 
					Scalar& uu) const ;

	/**
	 * Internal function intended to be used by \c Map::poisson_vect 
	 * and \c Map::poisson_vect_oohara . It constructs the sets of 
	 * parameters used for each scalar Poisson equation from the one for 
	 * the vectorial one.
	 * 
	 * @param para [input] : the \c Param  used for the resolution of 
	 * the vectorial Poisson equation : \\
	 * \c para.int()  maximum number of iteration.\\
	 * \c para.double(0)  relaxation parameter.\\
	 * \c para.double(1)  required precision. \\
	 * \c para.tenseur_mod()  source of the vectorial part at the previous 
	 * step.\\
	 * \c para.cmp_mod()  source of the scalar part at the previous 
	 * step.
	 * 
	 * @param i [input] number of the scalar Poisson equation that is being 
	 * solved (values from 0 to 2 for the componants of the vectorial part
	 * and 3 for the scalar one).
	 * 
	 * @return the pointer on the parameter set used for solving the scalar 
	 * Poisson equation labelled by \e i .
	 */
	virtual Param* donne_para_poisson_vect (Param& para, int i) const ;
	
	/**
	 * Not yet implemented.
	 */
	virtual void poisson_frontiere (const Cmp&, const Valeur&, int, int, 
					Cmp&, double = 0., double = 0.) const ;
	virtual void poisson_frontiere_double (const Cmp& source, 
			const Valeur& lim_func, const Valeur& lim_der, 
			int num_zone, Cmp& pot) const  ;

	/**
	 * Computes the solution of a Poisson equation in the shell .
	 * imposing a boundary condition at the surface of the star 
	 * 
	 * @param source [input] : source of the equation.
	 * @param limite [input] : \c limite[num_front]  contains the angular 
	 * function being the boudary condition.
	 * @param par [input] : parameters of the computation.
	 * @param pot [output] : result.
	 */
	virtual void poisson_interne (const Cmp& source, const Valeur& limite,
			 Param& par, Cmp& pot) const ;


	/** Computes the solution of a 2-D Poisson equation.
	 *  The 2-D Poisson equation writes
	 *  \f[
	 *	{\partial^2 u\over\partial r^2} + 
	 *	    {1\over r} {\partial u \over \partial r} + 
	 *	    {1\over r^2} {\partial^2 u\over\partial \theta^2} = 
	 *		\sigma \ . 
	 *  \f] 
	 *
	 *   @param source_mat [input] Compactly supported part of 
	 *	    the source \f$\sigma\f$ of the 2-D Poisson equation (typically
	 *	    matter terms)
	 *   @param source_quad [input] Non-compactly supported part of 
	 *	    the source \f$\sigma\f$ of the 2-D Poisson equation (typically
	 *	    quadratic terms)
	 *   @param par [input/output] Parameters to control the resolution : \\ 
	 *  \c par.get_double_mod(0)  : [output] constant \c lambda 
	 *	    such that the source of the equation effectively solved
	 *	    is \c source_mat + lambda * source_quad , in order to
	 *	    fulfill the virial identity GRV2.  \\
	 *  If the theta basis is \c T_SIN_I , the following arguments
	 *  are required: \\
	 *  \c par.get_int(0)  : [input] maximum number of iterations \\
	 *  \c par.get_double(0)  : [input] relaxation parameter \\
	 *  \c par.get_double(1)  : [input] required precision: the iterative
	 *	  method is stopped as soon as the relative difference between
	 *	 \f$u^J\f$ and \f$u^{J-1}\f$ is greater than \c par.get_double(1) \\
	 *  \c par.get_cmp_mod(0)  : [input/output] input : \c Cmp 
	 *		containing \f$s^{J-1}\f$ to 
	 *		start the iteration (if nothing is known a priori, 
	 *		this \c Cmp must be set to zero); output: value
	 *		of \f$s^{J-1}\f$, to used in a next call to the routine \\
	 *  \c par.get_int_mod(0)  : [output] number of iterations 
	 *				    actually used to get the solution.
	 *	 
	 *   @param uu [input/output] solution \e u  with the boundary condition 
	 *	    \e u =0 at spatial infinity. 
	 */
	virtual void poisson2d(const Cmp& source_mat, const Cmp& source_quad, 
			       Param& par, Cmp& uu) const ;

	/**
	 * Not yet implemented.
	 */
	virtual void dalembert(Param& par, Scalar& fJp1, const Scalar& fJ, 
			       const Scalar& fJm1, const Scalar& source) const ;




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
    friend Mtbl* map_et_fait_drdt(const Map* ) ;
    friend Mtbl* map_et_fait_stdrdp(const Map* ) ;
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
     Mtbl* map_et_fait_drdt(const Map* ) ;
     Mtbl* map_et_fait_stdrdp(const Map* ) ;
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

			//------------------------------------//
			//            class Map_log            //
			//------------------------------------//

#define AFFINE 0
#define LOG 1

/**
 * Logarithmic radial mapping. \ingroup (map)
 * 
 * This mapping is a variation of the affine one.
 *
 * In each domain the description can be either affine (cf. Map_af documentation) or 
 * logarithmic. In that case (implemented only in the shells) we have
 *
 *  \li \f$\ln r=\alpha \xi + \beta\f$,
 * where \f$\alpha\f$ and \f$\beta\f$ are constant in each domain. 
 * 
 * 
 */

class Map_log : public Map_radial {

    // Data :
    // ----
    private:
	/// Array (size: \c mg->nzone ) of the values of \f$\alpha\f$ in each domain
	Tbl alpha ;	 
	/// Array (size: \c mg->nzone ) of the values of \f$\beta\f$ in each domain
	Tbl beta ;
	/** Array (size: \c mg->nzone ) of the type of variable in each domain.	  
	 * The possibles types are AFFINE and LOG.
	 **/
	Itbl type_var ;

 public:
	/**
	 * Same as dxdr if the domains where the description is affine and 
	 * \f$ \partial x / \partial \ln r\f$ where it is logarithmic.
	 */
	
	Coord dxdlnr ;

 private:
	void set_coord() ;

    // Constructors, destructor : 
    // ------------------------
    public:
	/**
	 * Standard Constructor
	 * @param mgrille  [input] Multi-domain grid on which the mapping is defined
	 * @param r_limits [input] Tbl (size: number of domains + 1) of the
	 *			   value of \e r  at the boundaries of the various 
	 *			   domains : 
	 *			   \li \c r_limits[l] : inner boundary of the 
	 *				 domain no. \c l 
	 *			   \li \c r_limits[l+1] : outer boundary of the 
	 *				 domain no. \c l 
	 * @param type_var [input] Array (size: number of domains) defining the type f mapping in each domain.
	 */
	Map_log (const Mg3d& mgrille, const Tbl& r_limits, const Itbl& typevar) ;	
	
	
	Map_log (const Map_log& ) ;      ///< Copy constructor
	Map_log (const Mg3d&, FILE* ) ; ///< Constructor from a file (see \c sauve(FILE*)

	virtual ~Map_log() ;	      ///< Destructor
	
	/// Returns \f$\alpha\f$ in the domain \c l
	double get_alpha (int l) const {return alpha(l) ;} ;
	/// Returns \f$\beta\f$ in the domain \c l
	double get_beta (int l) const {return beta(l) ;} ;
	/// Returns the type of description in the domain \c l
	int get_type (int l) const {return type_var(l) ;} ;
	
	/**
	 * General elliptic solver. The field is zero at infinity.
	 *
	 * @param params [input] : the operators and variables to be uses.
	 * @param so [input] : the source.
	 * @param uu [output] : the solution.
	 **/
	void sol_elliptic (Param_elliptic& params, 
			   const Scalar& so,  Scalar& uu) const ;

	/**
	 * General elliptic solver.
	 * The equation is not solved in the compactified domain.
	 *
	 * @param params [input] : the operators and variables to be uses.
	 * @param so [input] : the source.
	 * @param uu [output] : the solution.
	 * @param val [input] : value at the last shell.
	 **/
	void sol_elliptic_no_zec (Param_elliptic& params, 
			   const Scalar& so,  Scalar& uu, double) const ;

	
	virtual void sauve(FILE*) const ; ///< Save in a file

	/// Assignment to an affine mapping. 
	virtual void operator=(const Map_af& mpa) ;

	
	virtual ostream& operator>> (ostream&) const ; ///< Operator >>

	/**
	 *  Returns the value of the radial coordinate \e r  for a given
	 *  \f$(\xi, \theta', \phi')\f$ in a given domain. 
	 *	@param l [input] index of the domain
	 *	@param xi [input] value of \f$\xi\f$
	 *	@param theta [input] value of \f$\theta'\f$
	 *	@param pphi [input] value of \f$\phi'\f$
	 *	@return value of \f$r=R_l(\xi, \theta', \phi')\f$
	 */
	virtual double val_r (int l, double xi, double theta, double pphi) const ;

	/**
	 * Computes the domain index \e l  and the value of \f$\xi\f$ corresponding
	 * to a point given by its physical coordinates \f$(r, \theta, \phi)\f$. 
	 *	@param rr [input] value of \e r 
	 *	@param theta [input] value of \f$\theta\f$
	 *	@param pphi [input] value of \f$\phi\f$
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of \f$\xi\f$
	 */
	virtual void val_lx (double rr, double theta, double pphi, int& l, double& xi) const ;
	
	/** Computes the domain index \e l  and the value of \f$\xi\f$ corresponding
	 * to a point given by its physical coordinates \f$(r, \theta, \phi)\f$. 
	 *	@param rr [input] value of \e r 
	 *	@param theta [input] value of \f$\theta\f$
	 *	@param pphi [input] value of \f$\phi\f$
	 *	@param par [] unused by the \c Map_af  version
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of \f$\xi\f$
	 */
	virtual void val_lx (double rr, double theta, double pphi, const Param& par, int& l, double& xi) const ;

	
	virtual bool operator== (const Map&) const ;   /// < Comparison operator
	
	/** Returns the value of the radial coordinate \e r  for a given
	 *  \f$\xi\f$ and a given collocation point in \f$(\theta', \phi')\f$ 
	 *   in a given domain. 
	 *	@param l [input] index of the domain
	 *	@param xi [input] value of \f$\xi\f$
	 *	@param j [input] index of the collocation point in \f$\theta'\f$
	 *	@param k [input] index of the collocation point in \f$\phi'\f$
	 *	@return value of \f$r=R_l(\xi, {\theta'}_j, {\phi'}_k)\f$
	 */
	virtual double val_r_jk (int l, double xi, int j, int k) const ;
	
	/** Computes the domain index \e l  and the value of \f$\xi\f$ corresponding
	 * to a point of arbitrary \e r  but collocation values of \f$(\theta, \phi)\f$ 
	 *	@param rr [input] value of \e r 
	 *	@param j [input] index of the collocation point in \f$\theta\f$
	 *	@param k [input] index of the collocation point in \f$\phi\f$
	 *	@param par [] unused by the \c Map_af  version
	 *	@param l [output] value of the domain index
	 *	@param xi [output] value of \f$\xi\f$
	 */
	virtual void val_lx_jk (double rr, int j, int k, const Param& par, int& l, double& xi) const ;
	
	/** Computes \f$\partial/ \partial r\f$ of a \c Scalar.
	 *  Note that in the  compactified external domain (CED), it computes
	 *  \f$-\partial/ \partial u = r^2 \partial/ \partial r\f$.
	 *  @param ci [input] field to consider
	 *  @param resu [output] derivative of \c ci 
	 */
	virtual void dsdr (const Scalar& ci, Scalar& resu) const ;
	
	/** Computes \f$\partial/ \partial r\f$ of a \c Scalar if the description is affine and 
	 * \f$\partial/ \partial \ln r\f$ if it is logarithmic.
	 *  Note that in the  compactified external domain (CED), the \c dzpuis 
	 *  flag of the output is 2 if the input has \c dzpuis  = 0, and 
	 *  is increased by 1 in other cases.
	 *  @param uu [input] field to consider
	 *  @param resu [output] derivative of \c uu 
	 */
	virtual void dsdradial (const Scalar& uu, Scalar& resu) const ;

	virtual void homothetie (double) ; /// < Not implemented
	virtual void resize (int, double) ;/// < Not implemented
	virtual void adapt (const Cmp&, const Param&) ;/// < Not implemented
	virtual void dsdr (const Cmp&, Cmp&) const ;/// < Not implemented
	virtual void srdsdt (const Cmp&, Cmp&) const ;/// < Not implemented
	virtual void srstdsdp (const Cmp&, Cmp&) const ;/// < Not implemented
	virtual void srdsdt (const Scalar&, Scalar&) const ;/// < Not implemented
	virtual void srstdsdp (const Scalar&, Scalar&) const ;/// < Not implemented
	virtual void dsdt (const Scalar&, Scalar&) const ;/// < Not implemented
	virtual void stdsdp (const Scalar&, Scalar&) const ;/// < Not implemented
	virtual void laplacien (const Cmp&, int, Cmp&) const ;/// < Not implemented
	virtual void lapang (const Scalar&, Scalar&) const ;/// < Not implemented
	virtual Tbl* integrale (const Cmp&) const ;/// < Not implemented
	virtual void poisson (const Cmp&, Param&, Cmp&) const ;/// < Not implemented
	virtual void poisson_regular (const Cmp&, int, int, double, Param&, Cmp&, Cmp&, Cmp&, 
				      Tenseur&, Cmp&, Cmp&) const ;/// < Not implemented
	virtual void poisson_angu (const Scalar&, Param&, Scalar&) const ;/// < Not implemented
	virtual Param* donne_para_poisson_vect (Param&, int) const ;/// < Not implemented
	virtual void poisson_frontiere (const Cmp&, const Valeur&, int, int, Cmp&, double = 0., double = 0.) const ;/// < Not implemented
	virtual void poisson_frontiere_double (const Cmp&, const Valeur&, const Valeur&, int, Cmp&) const ;/// < Not implemented
	virtual void poisson_interne (const Cmp&, const Valeur&, Param&, Cmp&) const ;/// < Not implemented
	virtual void poisson2d (const Cmp&, const Cmp&, Param&, Cmp&) const ;/// < Not implemented
	virtual void dalembert (Param&, Scalar&, const Scalar&, const Scalar&, const Scalar&) const ;/// < Not implemented


	// Building functions for the Coord's
	// ----------------------------------
	friend Mtbl* map_log_fait_r(const Map* ) ;
	friend Mtbl* map_log_fait_tet(const Map* ) ;
	friend Mtbl* map_log_fait_phi(const Map* ) ;
	friend Mtbl* map_log_fait_sint(const Map* ) ;
	friend Mtbl* map_log_fait_cost(const Map* ) ;
	friend Mtbl* map_log_fait_sinp(const Map* ) ;
	friend Mtbl* map_log_fait_cosp(const Map* ) ;
	
	friend Mtbl* map_log_fait_x(const Map* ) ;
	friend Mtbl* map_log_fait_y(const Map* ) ;
	friend Mtbl* map_log_fait_z(const Map* ) ;
	
	friend Mtbl* map_log_fait_xa(const Map* ) ;
	friend Mtbl* map_log_fait_ya(const Map* ) ;
	friend Mtbl* map_log_fait_za(const Map* ) ;
	
	friend Mtbl* map_log_fait_xsr(const Map* ) ;
	friend Mtbl* map_log_fait_dxdr(const Map* ) ;
	friend Mtbl* map_log_fait_drdt(const Map* ) ;
	friend Mtbl* map_log_fait_stdrdp(const Map* ) ;
	friend Mtbl* map_log_fait_srdrdt(const Map* ) ;
	friend Mtbl* map_log_fait_srstdrdp(const Map* ) ;
	friend Mtbl* map_log_fait_sr2drdt(const Map* ) ;
	friend Mtbl* map_log_fait_sr2stdrdp(const Map* ) ;
	friend Mtbl* map_log_fait_d2rdx2(const Map* ) ;
	friend Mtbl* map_log_fait_lapr_tp(const Map* ) ;
	friend Mtbl* map_log_fait_d2rdtdx(const Map* ) ;
	friend Mtbl* map_log_fait_sstd2rdpdx(const Map* ) ;
	friend Mtbl* map_log_fait_sr2d2rdt2(const Map* ) ;
	friend Mtbl* map_log_fait_dxdlnr(const Map* ) ;
	
};

Mtbl* map_log_fait_r(const Map* ) ;
Mtbl* map_log_fait_tet(const Map* ) ;
Mtbl* map_log_fait_phi(const Map* ) ;
Mtbl* map_log_fait_sint(const Map* ) ;
Mtbl* map_log_fait_cost(const Map* ) ;
Mtbl* map_log_fait_sinp(const Map* ) ;
Mtbl* map_log_fait_cosp(const Map* ) ;

Mtbl* map_log_fait_x(const Map* ) ;
Mtbl* map_log_fait_y(const Map* ) ;
Mtbl* map_log_fait_z(const Map* ) ;

Mtbl* map_log_fait_xa(const Map* ) ;
Mtbl* map_log_fait_ya(const Map* ) ;
Mtbl* map_log_fait_za(const Map* ) ;

Mtbl* map_log_fait_xsr(const Map* ) ;
Mtbl* map_log_fait_dxdr(const Map* ) ;
Mtbl* map_log_fait_drdt(const Map* ) ;
Mtbl* map_log_fait_stdrdp(const Map* ) ;
Mtbl* map_log_fait_srdrdt(const Map* ) ;
Mtbl* map_log_fait_srstdrdp(const Map* ) ;
Mtbl* map_log_fait_sr2drdt(const Map* ) ;
Mtbl* map_log_fait_sr2stdrdp(const Map* ) ;
Mtbl* map_log_fait_d2rdx2(const Map* ) ;
Mtbl* map_log_fait_lapr_tp(const Map* ) ;
Mtbl* map_log_fait_d2rdtdx(const Map* ) ;
Mtbl* map_log_fait_sstd2rdpdx(const Map* ) ;
Mtbl* map_log_fait_sr2d2rdt2(const Map* ) ;

Mtbl* map_log_fait_dxdlnr (const Map*) ;

#endif
