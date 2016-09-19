/*
 *  Definition of Lorene class Cmp
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2002 Eric Gourgoulhon
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *   Copyright (c) 2000-2002 Jerome Novak
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


#ifndef __CMP_H_ 
#define __CMP_H_ 


/*
 * $Id$
 * $Log$
 * Revision 1.23  2016/09/19 15:26:22  j_novak
 * Correction of several bugs preventing the shared library compilation.
 *
 * Revision 1.22  2014/10/13 08:52:33  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.21  2014/10/06 15:09:39  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.20  2012/08/12 17:35:36  p_cerda
 * Magnetstar: adding new member to class Cmp
 *
 * Revision 1.19  2010/02/02 13:34:12  e_gourgoulhon
 * Marked DEPRECATED (in the documentation).
 *
 * Revision 1.18  2005/08/30 08:35:10  p_grandclement
 * Addition of the Tau version of the vectorial Poisson equation for the Tensors
 *
 * Revision 1.17  2005/08/02 06:09:58  m_saijo
 * Modified comment lines (div_r, multi_r, mult_rsint, div_rsint)
 *
 * Revision 1.16  2004/12/29 16:20:20  k_taniguchi
 * Addition of the function poisson_ylm.
 *
 * Revision 1.15  2004/11/30 20:38:10  k_taniguchi
 * Addition of the function poisson_falloff.
 *
 * Revision 1.14  2004/10/11 15:08:59  j_novak
 * The radial manipulation functions take Scalar as arguments, instead of Cmp.
 * Added a conversion operator from Scalar to Cmp.
 * The Cmp radial manipulation function make conversion to Scalar, call to the
 * Map_radial version with a Scalar argument and back.
 *
 * Revision 1.13  2004/03/31 11:21:02  f_limousin
 * Method Cmp::poisson_neumann_interne has been implemented to solve
 * the continuity equation for strange stars.
 *
 * Revision 1.12  2004/03/22 13:12:40  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.11  2004/03/01 09:54:58  j_novak
 * Suppression of the Cmp version of avance_dalembert (now only with Scalar's)
 *
 * Revision 1.10  2003/11/06 14:43:37  e_gourgoulhon
 * Gave a name to const arguments in certain method prototypes (e.g.
 * constructors) to correct a bug of DOC++.
 *
 * Revision 1.9  2003/09/24 20:52:37  e_gourgoulhon
 * Added constructor by conversion of a Scalar.
 *
 * Revision 1.8  2003/08/26 09:46:10  j_novak
 * Added the method multipole_spectrum
 *
 * Revision 1.7  2003/06/20 14:16:10  f_limousin
 * Add the function compare().
 *
 * Revision 1.6  2003/06/20 09:27:09  j_novak
 * Modif commentaires.
 *
 * Revision 1.5  2002/09/13 09:17:31  j_novak
 * Modif. commentaires
 *
 * Revision 1.4  2002/06/17 14:05:16  j_novak
 * friend functions are now also declared outside the class definition
 *
 * Revision 1.3  2002/05/17 12:08:46  e_gourgoulhon
 * Corrected error in the comment about dzpuis: multiplied --> divided
 *
 * Revision 1.2  2002/01/03 15:30:27  j_novak
 * Some comments modified.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.101  2001/10/29  15:36:03  novak
 * Ajout de Cmp::div_r()
 *
 * Revision 2.100  2001/10/16 10:03:57  novak
 * *** empty log message ***
 *
 * Revision 2.99  2001/08/31 14:52:10  novak
 * Back to 2.97 version 2.98 was useless
 *
 * Revision 2.97  2001/07/19 14:01:39  novak
 * new arguments for Cmp::avance_dalembert
 *
 * Revision 2.96  2001/05/29 16:09:40  eric
 * Modif commentaires (mise en conformite Doc++ 3.4.7).
 *
 * Revision 2.95  2001/05/26  15:07:20  eric
 * Ajout de operator% : multiplication de deux Cmp avec desaliasage
 *
 * Revision 2.94  2001/05/25  09:30:07  phil
 * ajout de filtre_phi
 *
 * Revision 2.93  2001/03/30  13:36:22  phil
 * ajout de raccord_externe
 *
 * Revision 2.92  2001/03/26  08:11:50  eric
 * Modif commentaires.
 *
 * Revision 2.91  2001/03/22  10:25:19  phil
 * modification prototypage de raccord_zec.C
 *
 * Revision 2.90  2001/02/12  18:08:10  phil
 * ajout de Cmp::fixe_decroissance
 *
 * Revision 2.89  2000/12/13  14:50:05  phil
 * changement nom variable dzpuis dans raccord_c1_zec
 *
 * Revision 2.88  2000/12/13  14:35:53  phil
 * *** empty log message ***
 *
 * Revision 2.87  2000/12/13  14:26:42  phil
 * *** empty log message ***
 *
 * Revision 2.86  2000/12/13  14:25:26  phil
 * vire commentaires des raccords (provisioire)
 *
 * Revision 2.85  2000/12/13  14:19:49  phil
 * modif commentaires
 *
 * Revision 2.84  2000/12/13  14:08:36  phil
 * ajout procedure raccord_c1_zec
 *
 * Revision 2.83  2000/12/04  16:48:47  novak
 * *** empty log message ***
 *
 * Revision 2.82  2000/12/04 15:06:15  novak
 * *** empty log message ***
 *
 * Revision 2.81  2000/11/15 13:24:28  phil
 * modification de asymptot
 *
 * Revision 2.80  2000/11/15  13:19:13  phil
 * *** empty log message ***
 *
 * Revision 2.79  2000/11/15  13:17:01  phil
 * *** empty log message ***
 *
 * Revision 2.78  2000/11/15  13:15:45  phil
 * gestion affichage dans asymptot
 *
 * Revision 2.77  2000/10/20  09:43:30  phil
 * changement commentaires
 *
 * Revision 2.76  2000/10/19  14:07:06  novak
 * Ajout de la fonction membre avance_dalembert (experimentale)
 *
 * Revision 2.75  2000/10/19 09:20:36  phil
 * *** empty log message ***
 *
 * Revision 2.74  2000/10/19  09:13:45  phil
 * ajout des fonctions :
 * filtre(int)
 * set_val_inf(double)
 * set_val_hor(double,int)
 *
 * Revision 2.73  2000/10/05  14:18:14  eric
 * La fonction check_poisson est rebaptisee test_poisson.
 *
 * Revision 2.72  2000/10/05  13:56:52  eric
 * *** empty log message ***
 *
 * Revision 2.71  2000/10/05  13:52:25  eric
 * Ajout de la fonction check_poisson.
 *
 * Revision 2.70  2000/09/13  12:21:44  eric
 * Modif commentaires.
 *
 * Revision 2.69  2000/09/13  12:11:48  eric
 * Ajout de la fonction allocate_all().
 *
 * Revision 2.68  2000/09/07  15:26:40  keisuke
 * Add a new argument Cmp& uu in Cmp::poisson_regular.
 *
 * Revision 2.67  2000/09/04  09:11:06  keisuke
 * Suppress Cmp::poisson_regular (version without parameter).
 *
 * Revision 2.66  2000/08/31  13:04:30  eric
 * Ajout des fonctions mult_rsint et div_rsint.
 *
 * Revision 2.65  2000/08/29  13:51:36  keisuke
 * *** empty log message ***
 *
 * Revision 2.64  2000/08/29  13:46:14  keisuke
 * Add the polar and azimuthal derivatives of the diverging potential
 * in Cmp::poisson_regular.
 * Modify the argumants of Cmp::poisson_regular.
 *
 * Revision 2.63  2000/08/28  15:48:22  keisuke
 * Modif Cmp::poisson_regular
 *
 * Revision 2.62  2000/08/28  15:43:11  keisuke
 * Modif Cmp::poisson_regular
 *
 * Revision 2.61  2000/08/04  12:09:58  eric
 * Ajout de l'operator()(int l) et de la fonction set(int l) pour
 * l'acces aux Tbl individuels.
 *
 * Revision 2.60  2000/08/04  09:18:05  keisuke
 * Transformation Cmp::poisson_regular_param en Cmp::poisson_regular
 *
 * Revision 2.59  2000/08/03  14:01:29  keisuke
 * Modif Cmp::poisson_regular et ajout de Cmp::poisson_regular_param
 *
 * Revision 2.58  2000/07/29  12:50:01  keisuke
 * Modif Cmp::poisson_regular
 *
 * Revision 2.57  2000/07/20  13:33:50  keisuke
 * Modif Cmp::poisson_regular
 *
 * Revision 2.56  2000/07/20  10:25:09  keisuke
 * Modif Cmp::poisson_regular
 *
 * Revision 2.55  2000/07/19  15:50:23  keisuke
 * Ajout de Cmp::poisson_regular
 *
 * Revision 2.54  2000/05/22  14:38:32  phil
 * ajout de dec_dzpuis et inc_dzpuis
 *
 * Revision 2.53  2000/04/27  15:18:57  phil
 * *** empty log message ***
 *
 * Revision 2.52  2000/03/28  17:44:41  phil
 * Cmp::raccord() -> Cmp::raccord(int)
 *
 * Revision 2.51  2000/03/28  17:31:31  phil
 * *** empty log message ***
 *
 * Revision 2.50  2000/03/28  17:25:35  phil
 * ajout de Cmp::raccord()
 *
 * Revision 2.49  2000/03/25  12:52:45  eric
 * Ajout de la fonction asymptot(int ).
 *
 * Revision 2.48  2000/03/20  13:33:31  phil
 * commentaires
 *
 * Revision 2.47  2000/03/17  17:32:54  phil
 * *** empty log message ***
 *
 * Revision 2.46  2000/03/17  17:07:14  phil
 * *** empty log message ***
 *
 * Revision 2.45  2000/03/17  16:56:00  phil
 * ajout de poisson_dirichlet et de son amie poisson_neumann
 *
 * Revision 2.44  2000/03/06  10:55:44  eric
 * Ajout des methodes import_symy et import_asymy.
 *
 * Revision 2.43  2000/02/28  16:29:48  eric
 * Ajout des fonctions import_gal, import_align, import_anti.
 *
 * Revision 2.42  2000/01/28  16:08:55  eric
 * Ajout des fonctions dz_nonzero et check_dzpuis.
 *
 * Revision 2.41  2000/01/07  16:28:15  eric
 * Suppression de la fonction membre gradient.
 *
 * Revision 2.40  1999/12/21  13:03:22  eric
 * Changement de prototype de la routine poisson avec Param& : la solution est
 *  desormais passee en argument (et non plus en valeur de retour)
 *  pour permettre l'initialisation de methodes de resolution iteratives.
 *
 * Revision 2.39  1999/12/21  10:06:52  eric
 * Il y a desormais deux versions de poisson: une sans Param et une
 * avec Param.
 *
 * Revision 2.38  1999/12/10  16:19:33  eric
 * Modif commentaires.
 *
 * Revision 2.37  1999/12/10  15:59:01  eric
 * Modif commentaires fonction set.
 *
 * Revision 2.36  1999/12/09  10:45:54  eric
 * Ajout du calcul d'integrale (membre p_integ et fonctions
 *   integrale et integrale_domains).
 *
 * Revision 2.35  1999/12/08  12:38:38  eric
 * Ajout de la fonction import.
 *
 * Revision 2.34  1999/12/07  14:53:13  eric
 * Changement ordre des arguments (phi,theta,r) --> (r,theta,phi)
 *   dans la routine val_point.
 *
 * Revision 2.33  1999/12/06  16:47:00  eric
 * Ajout de la fonction val_point.
 *
 * Revision 2.32  1999/12/02  17:59:11  phil
 * *** empty log message ***
 *
 * Revision 2.31  1999/12/02  14:28:46  eric
 * Reprototypage de la fonction poisson(): const.
 * Commentaires.
 *
 * Revision 2.30  1999/11/30  14:20:54  eric
 * Reprototypage des fonctions membres mult_r, mult_r_zec,
 *  dec2_dzpuis et inc2_dzpuis : Cmp --> void.
 *
 * Revision 2.29  1999/11/29  13:18:06  eric
 * Modif commentaires.
 *
 * Revision 2.28  1999/11/29  12:56:11  eric
 * Introduction des membres p_lap, ind_lap.
 * Changement prototype de la fonction laplacien.
 *
 * Revision 2.27  1999/11/26  14:22:54  eric
 * Ajout du membre dzpuis et des fonctions de manipulation associees.
 *
 * Revision 2.26  1999/11/25  16:27:00  eric
 * Reorganisation complete du calcul et stokage des derivees partielles.
 *
 * Revision 2.25  1999/11/23  16:21:32  eric
 * Suppression du membre statique Cmp_Zero.
 * Suppression du constructeur par defaut.
 *
 * Revision 2.24  1999/11/22  16:48:00  phil
 * Cmp_Zero est desormais public
 *
 * Revision 2.23  1999/11/22  16:34:17  eric
 * Ajout de l'element global Cmp_Zero.
 *
 * Revision 2.22  1999/11/22  15:41:42  eric
 * Ajout des operateurs set(l,k,j,i) et (l,k,j,i).
 * Ajout de la fonction annule(int l).
 *
 * Revision 2.21  1999/11/15  14:12:28  eric
 * Ajout des fonctions mathematiques cos, sin, ..., min, max, norme,...
 *
 * Revision 2.20  1999/11/12  17:08:10  eric
 * Ajout de la partie manquante de l'arithmetique.
 *
 * Revision 2.19  1999/10/28  09:36:56  eric
 * Ajout de la fonction affiche_seuil.
 *
 * Revision 2.18  1999/10/28  09:01:24  eric
 * Constructeur par lecture de fichier.
 * Ajout de la fonction annule(int, int).
 *
 * Revision 2.17  1999/10/27  16:46:23  phil
 * ajout de mult_r_zec
 *
 * Revision 2.16  1999/10/27  15:38:40  eric
 * Suppression du membre c.
 *
 * Revision 2.15  1999/10/27  08:42:40  eric
 * Introduction du membre Valeur va.
 * Le pointeur Valeur* c est desormais un membre prive constant qui pointe
 * sur va.
 * Suppression de la fonction nouveau(), ainsi que du constructeur par
 * defaut.
 *
 * Revision 2.14  1999/10/22  08:14:19  eric
 * Depoussierage.
 * Documentation.
 *
 * Revision 2.13  1999/10/19  14:40:51  phil
 * ajout de inc2_dzpuis()
 *
 * Revision 2.12  1999/09/16  13:16:47  phil
 * ajout de Cmp mult_r()
 *
 * Revision 2.11  1999/09/15  10:29:44  phil
 * ajout de dec2_dzpuis()
 *
 * Revision 2.10  1999/09/14  17:13:05  phil
 * ajout de Cmp operator*(double,const Cmp&)
 *
 * Revision 2.9  1999/09/14  13:45:27  phil
 * suppression de la divergence
 *
 * Revision 2.8  1999/09/14  12:50:31  phil
 * ajout de Cmp deriv(int) et de Cmp divergence()
 *
 * Revision 2.7  1999/09/07  16:08:04  phil
 * ajout de la fonction membre gradient
 *
 * Revision 2.6  1999/09/06  14:50:27  phil
 * ajout du laplacien
 *
 * Revision 2.5  1999/09/06  14:35:05  phil
 * ajout de poisson
 *
 * Revision 2.4  1999/03/03  11:13:46  hyc
 * *** empty log message ***
 *
 * Revision 2.3  1999/03/03  11:07:27  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

#include <cstdio>

#include "valeur.h"
#include "map.h"

namespace Lorene {
class Param ; 

/**
 * Component of a tensorial field *** DEPRECATED : use class \c Scalar instead ***. \ingroup (otens)
 */

class Cmp {

    // Data : 
    // -----
    private:
	const Map* mp ;	    ///< Reference mapping

	/// Logical state (\c ETATNONDEF , \c ETATQCQ  or \c ETATZERO ).
	int etat ;	    

	/**
	 * Power of \e r  by which the quantity represented by \c this  
	 * must be divided in the external compactified zone in order 
	 * to get the correct physical values
	 */
	int dzpuis ;	

    public:
	Valeur va ;		///< The numerical value of the \c Cmp     

    // Derived data : 
    // ------------
    private:
	/// Pointer on \f$\partial/\partial r\f$ of \c *this 
	mutable Cmp* p_dsdr ;	
	/// Pointer on \f$1/r \partial/\partial \theta\f$ of \c *this 
	mutable Cmp* p_srdsdt ;	
	/// Pointer on \f$1/(r\sin\theta) \partial/\partial \phi\f$ of \c *this 
	mutable Cmp* p_srstdsdp ;
	
	/** Pointer on \f$\partial/\partial x\f$ of \c *this ,
	 *  where \f$x=r\sin\theta \cos\phi\f$
	 */
	mutable Cmp* p_dsdx ;	

	/** Pointer on \f$\partial/\partial y\f$ of \c *this ,
	 *  where \f$y=r\sin\theta \sin\phi\f$
	 */
	mutable Cmp* p_dsdy ;	

	/** Pointer on \f$\partial/\partial z\f$ of \c *this ,
	 *  where \f$z=r\cos\theta\f$
	 */
	mutable Cmp* p_dsdz ;	

	/** Pointer on the Laplacian of \c *this 
	 */
	mutable Cmp* p_lap ;	

	/** Power of \e r  by which the last computed Laplacian has been 
	 *  multiplied in the external compactified domain.  
	 */
	mutable int ind_lap ; 

	/** Pointer on the space integral of \c *this  (values in each 
	 *  domain)
	 */
	mutable Tbl* p_integ ; 

    // Constructors - Destructor
    // -------------------------
	
    public:
	explicit Cmp(const Map& map) ;	///< Constructor from mapping
	explicit Cmp(const Map* p_map) ;	///< Constructor from mapping
	Cmp(const Cmp& a) ;		///< Copy constructor

	/// Constructor from a file (see \c sauve(FILE*) )
	Cmp(const Map&, const Mg3d&, FILE* ) ;    		

	~Cmp() ;			///< Destructor

    // Assignment
    // -----------
    public: 
	/// Assignment to another \c Cmp  defined on the same mapping
	void operator=(const Cmp& a) ;	

	void operator=(const Valeur& a) ; ///< Assignment to a \c Valeur 
	void operator=(const Mtbl& a) ;	 ///< Assignment to a \c Mtbl 
	void operator=(double ) ;	 ///< Assignment to a \c double 
	void operator=(int ) ;		 ///< Assignment to an \c int 
	    
	/** Assignment to another \c Cmp  defined on a different mapping.
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original \c Cmp . 
	 *	@param ci [input] \c Cmp  to be imported.
	 */
	void import(const Cmp& ci) ;	 

	/** Assignment to another \c Cmp  defined on a different mapping.
	 *  Case where the \c Cmp  is symmetric with respect to the plane y=0.
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original \c Cmp . 
	 *	@param ci [input] \c Cmp  to be imported.
	 */
	void import_symy(const Cmp& ci) ;	 

	/** Assignment to another \c Cmp  defined on a different mapping.
	 *  Case where the \c Cmp  is antisymmetric with respect to the 
	 *  plane y=0.
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original \c Cmp . 
	 *	@param ci [input] \c Cmp  to be imported.
	 */
	void import_asymy(const Cmp& ci) ;	 

	/** Assignment to another \c Cmp  defined on a different mapping.
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original \c Cmp . 
	 *	@param nzet [input] Number of domains of the destination
	 *			    mapping (i.e. \c this->mp ) where the 
	 *			    importation is performed: the assignment
	 *			    is done for the domains whose indices are
	 *			    between 0 and \c nzet-1 . In the other
	 *			    domains, \c *this  is set to zero. 
	 *	@param ci [input] \c Cmp  to be imported.
	 */
	void import(int nzet, const Cmp& ci) ;	 

	/** Assignment to another \c Cmp  defined on a different mapping.
	 *  Case where the \c Cmp  is symmetric with respect to the plane y=0.
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original \c Cmp . 
	 *	@param nzet [input] Number of domains of the destination
	 *			    mapping (i.e. \c this->mp ) where the 
	 *			    importation is performed: the assignment
	 *			    is done for the domains whose indices are
	 *			    between 0 and \c nzet-1 . In the other
	 *			    domains, \c *this  is set to zero. 
	 *	@param ci [input] \c Cmp  to be imported.
	 */
	void import_symy(int nzet, const Cmp& ci) ;	 

	/** Assignment to another \c Cmp  defined on a different mapping.
	 *  Case where the \c Cmp  is antisymmetric with respect to the 
	 *  plane y=0.
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original \c Cmp . 
	 *	@param nzet [input] Number of domains of the destination
	 *			    mapping (i.e. \c this->mp ) where the 
	 *			    importation is performed: the assignment
	 *			    is done for the domains whose indices are
	 *			    between 0 and \c nzet-1 . In the other
	 *			    domains, \c *this  is set to zero. 
	 *	@param ci [input] \c Cmp  to be imported.
	 */
	void import_asymy(int nzet, const Cmp& ci) ;	 

    private:
	/** Assignment to another \c Cmp  defined on a different mapping,
	 *  when the two mappings do not have a particular relative orientation.
	 *
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original \c Cmp . 
	 *	@param nzet [input] Number of domains of the destination
	 *			    mapping (i.e. \c this->mp ) where the 
	 *			    importation is performed: the assignment
	 *			    is done for the domains whose indices are
	 *			    between 0 and \c nzet-1 . In the other
	 *			    domains, \c *this  is set to zero. 
	 *	@param ci [input] \c Cmp  to be imported.
	 */
	void import_gal(int nzet, const Cmp& ci) ;	 

	/** Assignment to another \c Cmp  defined on a different mapping,
	 *  when the two mappings have aligned Cartesian axis. 
	 *
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original \c Cmp . 
	 *	@param nzet [input] Number of domains of the destination
	 *			    mapping (i.e. \c this->mp ) where the 
	 *			    importation is performed: the assignment
	 *			    is done for the domains whose indices are
	 *			    between 0 and \c nzet-1 . In the other
	 *			    domains, \c *this  is set to zero. 
	 *	@param ci [input] \c Cmp  to be imported.
	 */
	void import_align(int nzet, const Cmp& ci) ;	 

	/** Assignment to another \c Cmp  defined on a different mapping,
	 *  when the two mappings have anti-aligned Cartesian axis (i.e.
	 *  \f$x_1 = - x_2\f$,  \f$y_1 = - y_2\f$,  \f$z_1 = z_2\f$). 
	 *
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original \c Cmp . 
	 *	@param nzet [input] Number of domains of the destination
	 *			    mapping (i.e. \c this->mp ) where the 
	 *			    importation is performed: the assignment
	 *			    is done for the domains whose indices are
	 *			    between 0 and \c nzet-1 . In the other
	 *			    domains, \c *this  is set to zero. 
	 *	@param ci [input] \c Cmp  to be imported.
	 */
	void import_anti(int nzet, const Cmp& ci) ;	 

	/** Assignment to another \c Cmp  defined on a different mapping,
	 *  when the two mappings have aligned Cartesian axis. 
	 *  Case where the \c Cmp  is symmetric with respect to the plane y=0.
	 *
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original \c Cmp . 
	 *	@param nzet [input] Number of domains of the destination
	 *			    mapping (i.e. \c this->mp ) where the 
	 *			    importation is performed: the assignment
	 *			    is done for the domains whose indices are
	 *			    between 0 and \c nzet-1 . In the other
	 *			    domains, \c *this  is set to zero. 
	 *	@param ci [input] \c Cmp  to be imported.
	 */
	void import_align_symy(int nzet, const Cmp& ci) ;	 

	/** Assignment to another \c Cmp  defined on a different mapping,
	 *  when the two mappings have anti-aligned Cartesian axis (i.e.
	 *  \f$x_1 = - x_2\f$,  \f$y_1 = - y_2\f$,  \f$z_1 = z_2\f$). 
	 *  Case where the \c Cmp  is symmetric with respect to the plane y=0.
	 *
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original \c Cmp . 
	 *	@param nzet [input] Number of domains of the destination
	 *			    mapping (i.e. \c this->mp ) where the 
	 *			    importation is performed: the assignment
	 *			    is done for the domains whose indices are
	 *			    between 0 and \c nzet-1 . In the other
	 *			    domains, \c *this  is set to zero. 
	 *	@param ci [input] \c Cmp  to be imported.
	 */
	void import_anti_symy(int nzet, const Cmp& ci) ;	 

	/** Assignment to another \c Cmp  defined on a different mapping,
	 *  when the two mappings have aligned Cartesian axis. 
	 *  Case where the \c Cmp  is antisymmetric with respect to the 
	 *  plane y=0.
	 *
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original \c Cmp . 
	 *	@param nzet [input] Number of domains of the destination
	 *			    mapping (i.e. \c this->mp ) where the 
	 *			    importation is performed: the assignment
	 *			    is done for the domains whose indices are
	 *			    between 0 and \c nzet-1 . In the other
	 *			    domains, \c *this  is set to zero. 
	 *	@param ci [input] \c Cmp  to be imported.
	 */
	void import_align_asymy(int nzet, const Cmp& ci) ;	 

	/** Assignment to another \c Cmp  defined on a different mapping,
	 *  when the two mappings have anti-aligned Cartesian axis (i.e.
	 *  \f$x_1 = - x_2\f$,  \f$y_1 = - y_2\f$,  \f$z_1 = z_2\f$). 
	 *  Case where the \c Cmp  is antisymmetric with respect to the 
	 *  plane y=0.
	 *
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original \c Cmp . 
	 *	@param nzet [input] Number of domains of the destination
	 *			    mapping (i.e. \c this->mp ) where the 
	 *			    importation is performed: the assignment
	 *			    is done for the domains whose indices are
	 *			    between 0 and \c nzet-1 . In the other
	 *			    domains, \c *this  is set to zero. 
	 *	@param ci [input] \c Cmp  to be imported.
	 */
	void import_anti_asymy(int nzet, const Cmp& ci) ;	 



    // Access to individual elements
    // -----------------------------
    public:

	/** Read/write of the value in a given domain.
	 * NB: to gain in efficiency, the method \c del_deriv()  (to delete
	 *     the derived members) is not called by this function. It must
	 *     thus be invoqued by the user.  
	 *
	 * @param l [input] domain index
	 * @return Tbl containing the value of the field in domain \c l .
	 */ 
	Tbl& set(int l) {
	    assert(etat == ETATQCQ) ;
	    return va.set(l) ;
	};
	
	/** Read-only of the value in a given domain.
	 * @param l [input] domain index
	 * @return Tbl containing the value of the field in domain \c l .
	 */ 
	const Tbl& operator()(int l) const {
	    assert(etat == ETATQCQ) ;
	    return va(l) ;
	};


	/** Read/write of a particular element.
	 * NB: to gain in efficiency, the method \c del_deriv()  (to delete
	 *     the derived members) is not called by this function. It must
	 *     thus be invoqued by the user.  
	 *     
	 * @param l [input] domain index
	 * @param k [input] \f$\phi\f$ index
	 * @param j [input] \f$\theta\f$ index
	 * @param i [input] \e r  (\f$\xi\f$) index
	 */ 
	double& set(int l, int k, int j, int i) {
	    assert(etat == ETATQCQ) ;
	    return va.set(l, k, j, i) ;
	};
	
	
	/** Read-only of a particular element.
	 * @param l [input] domain index
	 * @param k [input] \f$\phi\f$ index
	 * @param j [input] \f$\theta\f$ index
	 * @param i [input] \e r  (\f$\xi\f$) index
	 */ 
	double operator()(int l, int k, int j, int i) const {
	    assert(etat != ETATNONDEF) ;
	    if (etat == ETATZERO) {
		double zero = 0. ;
		return zero ; 
	    }
	    else{ 	    
		return va(l, k, j, i) ;
	    }
	};

	/** Computes the value of the field represented by \c *this  at an
	*   arbitrary point \f$(r, \theta, \phi)\f$, by means of the spectral 
	*   expansion.
	*	 @param r [input] value of the coordinate \e r 
	*	 @param theta [input] value of the coordinate \f$\theta\f$
	*	 @param phi [input] value of the coordinate \f$\phi\f$
	*	 @return value at the point \f$(r, \theta, \phi)\f$ 
	*		 of the field represented by \c *this . 
	*/
	double val_point(double r, double theta, double phi) const ; 


    // Memory management
    // -----------------
    private:
	void del_t() ;		    ///< Logical destructor
	void del_deriv() ;	    ///< Logical destructor of the derivatives
	void set_der_0x0() ;	    ///< Sets the pointers for derivatives to 0x0

    public:

    /**
     * Sets the logical state to \c ETATNONDEF  (undefined). 
     * Calls the logical destructor of the \c Valeur \c va  and
     * deallocates the memory occupied by all the derivatives. 
     */
	void set_etat_nondef() ;   

    /**
     * Sets the logical state to \c ETATZERO (zero). 
     * Calls the logical destructor of the \c Valeur \c va  and
     * deallocates the memory occupied by all the derivatives. 
     */
	void set_etat_zero() ;	    
	
    /**
     * Sets the logical state to \c ETATQCQ  (ordinary state).
     * If the state is already \c ETATQCQ , this function does nothing.
     * Otherwise, it calls the logical destructor of the \c Valeur \c va  and
     * deallocates the memory occupied by all the derivatives.
     */
	void set_etat_qcq() ;	    

    /**
     * Sets the logical state to \c ETATQCQ  (ordinary state)
     *  and performs the memory allocation of all the 
     *  elements, down to the \c double  arrays of the \c Tbl s. 
     *  This function performs in fact recursive calls to \c set_etat_qcq() 
     *  on each element of the chain \c Cmp  ->
     *  \c Valeur  -> \c Mtbl  -> \c Tbl . 
     */
	void allocate_all() ; 

    /**
     * Sets the \c Cmp  to zero in a hard way. 
     * 1/ Sets the logical state to \c ETATQCQ , i.e. to an ordinary state.
     * 2/ Fills the \c Valeur \c va  with zeros. 
     * NB: this function must be used for debugging purposes only.
     * For other operations, the functions \c set_etat_zero()
     * or \c annule(int, int) must be perferred. 
     */
	void annule_hard() ;

    /**
     * Sets the \c Cmp  to zero in a given domain.
     *	@param l [input]  Index of the domain in which the \c Cmp 
     *			  will be set (logically) to zero.
     */
	void annule(int l) ; 

    /**
     * Sets the \c Cmp  to zero in several domains.
     *	@param l_min [input] The \c Cmp  will be set (logically) to zero
     *			     in the domains whose indices are in the range
     *			     \c [l_min,l_max].
     *	@param l_max [input] see the comments for \c l_min .
     * 
     * Note that \c annule(0,va.mg->get_nzone()-1) is equivalent to
     *	 \c set_etat_zero() .
     */
	void annule(int l_min, int l_max) ; 
    
    /**
     * Sets the \c n  lasts coefficients in \e r  to 0 in the external domain.
     */
	void filtre (int n) ;
    
    /**
     * Sets the \c n  lasts coefficients in \f$\Phi\f$ to 0 in the 
     * domain \c zone .
     */
	void filtre_phi (int n, int zone) ;
    
    /**
     * Sets the value of the \c Cmp  to \c val  at infinity. This is usefull
     * for dealing with undefined values. The external domain must be 
     * compactified.
     */
	void set_val_inf (double val) ;
    
    /**
     * Sets the value of the \c Cmp  to \c val  on the inner boudary of the
     * shell number \c zone .This is usefull
     * for dealing with undefined values.
     */
	void set_val_hor (double val, int zone) ;
    /**
     * Substracts all the components behaving like \f$r^{-n}\f$ in the external 
     * domain, with \e n  strictly lower than \c puis , so that \c *this  
     * decreases at least like \f$r^{\tt puis} \f$ at infinity.
     */
	void fixe_decroissance (int puis) ;

    /**
     * Gives the spectrum in terms of multipolar modes \e l .
     *  @return a \c Tbl  of size (nzone, lmax), where lmax is the
     *  maximal multipolar momentum over all domains. The \e l -th
     *  element contains the L1 norm of the \e l -th multipole 
     *  (i.e. a sum over all \e m  of the norms (coefficient space)
     *  of the component of a given \f$Y_l^m\f$.
     */
	Tbl multipole_spectrum () ;
	
    // Extraction of information
    // -------------------------
    public:
	/// Returns the logical state
	int get_etat() const {return etat;} ;	
	/// Returns the mapping
	const Map* get_mp() const {return mp;};	 
	/// Returns \c dzpuis 
	int get_dzpuis() const {return dzpuis;} ;
	
	/** Returns \c true  if the last domain is compactified and
	 *  \c *this  is not zero in this domain
	 */
	bool dz_nonzero() const ; 
	
	/** Returns \c false  if the last domain is compactified 
	 *  and \c *this  is not zero in this domain and \c dzpuis 
	 *  is not equal to \c dzi , otherwise return true. 
	 */
	bool check_dzpuis(int dzi) const ; 
	
	
    // Outputs
    // -------
    public:
	void sauve(FILE *) const ;	    ///< Save in a file
    
	/** Prints only the values greater than a given threshold.
	 *   @param ostr [input] Output stream used for the printing
	 *   @param type [input] Type of display : 0 = prints only the
	 *     coefficients,  1 = prints only the values in configuration 
	 *     space, 2 = prints both
	 *   @param precision [input] Number of printed digits (default: 4)
	 *   @param threshold [input] Value above which an array element is printed
	 *    (default: 1.e-7)
	 */
	void affiche_seuil(ostream& ostr, int type = 0, int precision = 4, 
			   double threshold = 1.e-7) const ;

	/// Display
	friend ostream& operator<<(ostream& , const Cmp & ) ;	


    // Member arithmetics
    // ------------------
    public:
	void operator+=(const Cmp &) ;		    ///< += Cmp
	void operator-=(const Cmp &) ;		    ///< -= Cmp
	void operator*=(const Cmp &) ;		    ///< *= Cmp

    // Manipulation of spectral bases
    // ------------------------------    
    /** Sets the spectral bases of the \c Valeur \c va  to the standard ones 
     *  for a scalar
     */
    void std_base_scal() ;	 


    // Differential operators and others
    // ---------------------------------
    public:
	/** Returns \f$\partial / \partial r\f$ of \c *this .
	 *  Note that in the external compactified domain (ZEC), it returns
	 *  instead \f$r^2 \partial/ \partial r\f$.
	 */
	const Cmp& dsdr() const ; 
	
	/** Returns \f$1/r \partial / \partial \theta\f$ of \c *this .
	 *  Note that in the external compactified domain (ZEC), it returns
	 *  instead \f$r \partial/ \partial \theta\f$.
	 */
	const Cmp& srdsdt() const ; 

	/** Returns \f$1/(r\sin\theta) \partial / \partial \phi\f$ of \c *this .
	 *  Note that in the external compactified domain (ZEC), it returns
	 *  instead \f$r/\sin\theta \partial/ \partial \phi\f$.
	 */
	const Cmp& srstdsdp() const ; 

	/** Returns \f$\partial/\partial x\f$ of \c *this ,
	 *  where \f$x=r\sin\theta \cos\phi\f$.
	 *  Note that in the external compactified domain (ZEC), it returns
	 *  instead \f$r^2 \partial/ \partial x\f$.
	 */
	const Cmp& dsdx() const ;	

	/** Returns \f$\partial/\partial y\f$ of \c *this ,
	 *  where \f$y=r\sin\theta \sin\phi\f$.
	 *  Note that in the external compactified domain (ZEC), it returns
	 *  instead \f$r^2 \partial/ \partial y\f$.
	 */
	const Cmp& dsdy() const ;	

	/** Returns \f$\partial/\partial z\f$ of \c *this ,
	 *  where \f$z=r\cos\theta\f$.
	 *  Note that in the external compactified domain (ZEC), it returns
	 *  instead \f$r^2 \partial/ \partial z\f$.
	 */
	const Cmp& dsdz() const ;	

	/** Returns \f$\partial/\partial x_i\f$ of \c *this ,
	 *  where \f$x_i = (x, y, z)\f$.
	 *  Note that in the external compactified domain (ZEC), it returns
	 *  instead \f$r^2 \partial/ \partial x_i\f$.
	 *  @param i [input] i=0 for \e x ,  i=1 for \e y , i=2 for \e z .
	 */
	const Cmp& deriv(int i) const ;	

	/** Returns the Laplacian of \c *this 
	 *   @param zec_mult_r [input] Determines the quantity computed in
	 *			 the external compactified domain (ZEC) 
	 *		(\e u  in the field represented by \c *this ) :  \\
	 *		    zec_mult_r = 0 : \f$\Delta u\f$	\\
	 *		    zec_mult_r = 2 : \f$r^2 \,  \Delta u\f$	\\
	 *		    zec_mult_r = 4 (default) : \f$r^4 \, \Delta u\f$	
	 */
	const Cmp& laplacien(int zec_mult_r = 4) const ; 

        /// Division by \e r  everywhere.
	void div_r() ;

        /// Multiplication by \e r  everywhere.
	void mult_r() ;

	/** Multiplication by \e r  in the external compactified domain (ZEC)
	 */
	void mult_r_zec() ;

	/// Multiplication by \f$r\sin\theta\f$
	void mult_rsint() ;

        /// Multiplication by \f\cos\theta\f$
	void mult_cost() ;

	/// Division by \f$r\sin\theta\f$
	void div_rsint() ;

	/** Decreases by 1 the value of \c dzpuis  and changes accordingly
	 *  the values of the \c Cmp  in the external compactified domain (ZEC).
	 */
	void dec_dzpuis() ; 

	/** Increases by the value of \c dzpuis  and changes accordingly
	 *  the values of the \c Cmp  in the external compactified domain (ZEC).
	 */
	void inc_dzpuis() ; 
	
	/** Decreases by 2 the value of \c dzpuis  and changes accordingly
	 *  the values of the \c Cmp  in the external compactified domain (ZEC).
	 */
	void dec2_dzpuis() ; 

	/** Increases by 2 the value of \c dzpuis  and changes accordingly
	 *  the values of the \c Cmp  in the external compactified domain (ZEC).
	 */
	void inc2_dzpuis() ; 

	void set_dzpuis(int ) ;  ///< Set a value to \c dzpuis 

	/** Computes the integral over all space of \c *this .
	 *  The computed quantity is (\e u  being the field represented by
	 *   \c *this )
	 *    \f$\int u \, r^2 \sin\theta \,  dr\, d\theta \, d\phi\f$.
	 *  Note that in the external compactified domain (ZEC), \c dzpuis  
	 *  must be 4 for the computation to take place. 
	 */
	double integrale() const ; 
	
	/** Computes the integral in each domain of \c *this .
	 *  The computed quantity is (\e u  being the field represented by
	 *   \c *this )
	 *    \f$\int u \, r^2 \sin\theta \,  dr\, d\theta \, d\phi\f$
	 *  in each domain. The result is returned a \c Tbl  on the 
	 *  various domains. 
	 *  Note that in the external compactified domain (ZEC), \c dzpuis  
	 *  must be 4 for the computation to take place. 
	 */
	const Tbl& integrale_domains() const ; 
	
	/** Asymptotic expansion at r = infinity. 
	 * 
	 *  Determines the coefficients \f$a_k(\theta, \phi)\f$ of the expansion
	 *  \f[
	 *	\sum_{k=0}^n {a_k(\theta, \phi) \over r^k}
	 *  \f] 
	 *  of \c *this  when \f$r \rightarrow \infty\f$. 
	 *
	 *	@param n order of the expansion
	 *	@param flag : output
	 *	@return Array of n+1 \c Valeur s on \c mg->angu  
	 *		describing the coefficients \f$a_k(\theta, \phi)\f$. 
	 *		This array is allocated by the routine. 
	 * 
	 */
	Valeur** asymptot(int n, const int flag = 0) const ; 
	

    // PDE resolution 
    // --------------
    public:
	/** Solves the scalar Poisson equation with \c *this  as a source.
	 *   The source \f$\sigma\f$ of the equation \f$\Delta u = \sigma\f$ is 
	 *   represented by the \c Cmp  \c *this . 
	 *   Note that \c dzpuis  must be equal to 2, 3 or 4, i.e. that the
	 *   quantity stored in \c *this  is in fact \f$r^2 \sigma\f$ or
	 *   \f$r^4 \sigma\f$ in the external compactified domain. 
	 *   The solution \e u  with the boundary condition \e u =0 at spatial
	 *   infinity is the returned \c Cmp . 
	 */
	Cmp poisson() const ;
	
	/**  Same as Poisson with a Tau method
	 */
	Cmp poisson_tau() const ;

	Cmp poisson_falloff(int k_falloff) const ;

	Cmp poisson_ylm(int nylm, double* intvec) const ;

	/** Solves the scalar Poisson equation with \c *this  as a source
	 *   (version with parameters to control the resolution).
	 *   The source \f$\sigma\f$ of the equation \f$\Delta u = \sigma\f$ is 
	 *   represented by the \c Cmp  \c *this . 
	 *   Note that \c dzpuis  must be equal to 2 or 4, i.e. that the
	 *   quantity stored in \c *this  is in fact \f$r^2 \sigma\f$ or
	 *   \f$r^4 \sigma\f$ in the external compactified domain. 
	 *   @param par [input/output] possible parameters
	 *   @param uu [input/output] solution \e u  with the boundary condition 
	 *   \e u =0 at spatial infinity. 
	 */
	void poisson(Param& par, Cmp& uu) const ;
	
	/**  Same as Poisson with a Tau method
	 */
	void poisson_tau(Param& par, Cmp& uu) const ;
	
	void poisson_falloff(Param& par, Cmp& uu, int k_falloff) const ;

	void poisson_ylm(Param& par, Cmp& uu, int nylm, double* intvec) const ;

	/**
	 * Is identicall to \c Cmp::poisson() . The regularity condition at the 
	 * origin is replace by a boundary condition of the Dirichlet type.
	 * 
	 * @param limite [input] : angular function. The boundary condition is 
	 * given by \c limite[num] .
	 * @param num [input] : index of the boudary at which the condition is to 
	 * be fullfilled.
	 * 
	 * More precisely we impose the solution is equal to \c limite[num]  at the
	 * boundary between the domains \c num  and \c num+1  (the latter one being 
	 * a shell).
	 * 
	 */
	Cmp poisson_dirichlet (const Valeur& limite, int num) const ;
	
	/**
	 * Idem as \c Cmp::poisson_dirichlet , the boundary condition being on 
	 * the radial derivative of the solution.
	 */
	Cmp poisson_neumann (const Valeur&, int) const ;

	/**
	 * Idem as \c Cmp::poisson_neumann , the boundary condition is on 
	 * the radial derivative of the solution. But in this method, the
	 * poisson equation is solved in the shell only. We have so to
	 * impose a boundary condition on the surface of the star. 
	 * This is used for example to solve the continuity equation 
	 * for the fluid in the star. 
	 */
	Cmp poisson_neumann_interne (const Valeur&, Param& par, Cmp& resu) const ;
	Cmp poisson_frontiere_double   (const Valeur&, const Valeur&, int) const ;

	/** Solves the scalar Poisson equation with \c *this  as a source
	 *   (version with parameters to control the resolution).
	 *   The source \f$\sigma\f$ of the equation \f$\Delta u = \sigma\f$ is 
	 *   represented by the \c Cmp  \c *this . 
	 *   The regularized source
	 *   \f$\sigma_{\rm regu} = \sigma - \sigma_{\rm div}\f$
	 *   is constructed and solved.
	 *   Note that \c dzpuis  must be equal to 2 or 4, i.e. that the
	 *   quantity stored in \c *this  is in fact \f$r^2 \sigma\f$ or
	 *   \f$r^4 \sigma\f$ in the external compactified domain.
	 *   @param k_div [input] regularization degree of the procedure
	 *   @param nzet [input] number of domains covering the star
	 *   @param unsgam1 [input] parameter \f$1/(\gamma-1)\f$ where \f$\gamma\f$
	 *          denotes the adiabatic index
	 *   @param par [input/output] possible parameters
	     @param uu [input/output] solution
	 *   @param uu_regu [output] solution of the regular part of
	 *          the source.
	 *   @param uu_div [output] solution of the diverging part of
	 *          the source.
	 *   @param duu_div [output] derivative of the diverging potential.
	 *   @param source_regu [output] regularized source
	 *   @param source_div [output] diverging part of the source
	 */
	void poisson_regular(int k_div, int nzet, double unsgam1, Param& par,
			     Cmp& uu, Cmp& uu_regu, Cmp& uu_div,
			     Tenseur& duu_div,
			     Cmp& source_regu, Cmp& source_div) const ;

	/** Checks if a Poisson equation with \c *this  as a source
	 *  has been correctly solved.
	 * 
	 *  @param uu [input] Solution \e u  of the Poisson equation
	 *		      \f$\Delta u = \sigma\f$,  \f$\sigma\f$ being 
	 *		      represented by the \c Cmp  \c *this .
	 * 
	 *  @param ostr [input/output] Output stream used for displaying
	 *		\c err .
	 *
	 *  @param detail [input] \li if \c true  displays \c err(0,*) , 
	 *		    \c err(1,*) and \c err(2,*) 
	 *		\li if \c false (default),  displays only 
	 *		the relative error \c err(0,*). 
	 *  
	 *  @return 2-D \c Tbl  \c err decribing the errors in each 
	 *	    domain: 
	 *	\li \c err(0,l) :  Relative error in domain no. \c l , 
	 *	    defined as the maximum value of 
	 *	    \f$|\Delta u - \sigma|\f$ in that domain divided by \e m , 
	 *	    where \e m  is the maximum value of \f$|\sigma|\f$ 
	 *	    over all domains if \c dzpuis = 0} or \f$\sigma\f$ is
	 *	    zero in the external compactified domain (ECD). If 
	 *	    \c dzpuis != 0} and \f$\sigma\f$ does not vanish in the 
	 *	    ECD, the value of \e m  used in the
	 *	    non-compactified domains is the maximum value over
	 *	    these domains, whereas the value of \e m  used in the
	 *	    external compactified domain is the maximum value
	 *	    on that particular domain. 
	 *	\li \c err(1,l) :   Maximum value of the absolute error
	 *			\f$|\Delta u - \sigma|\f$ in domain no. \c l  
	 *	\li \c err(2,l) :   Maximum value of \f$|\sigma|\f$ in domain 
	 *			    no. \c l  
	 */
	Tbl test_poisson(const Cmp& uu, ostream& ostr, 
					bool detail = false) const ;  	
	/**
	 * Performs the \f$C^n\f$ matching of the nucleus with respect to the 
	 * first shell.
	 */
	void raccord(int n) ;
	
	/**
	 * Performs the \f$C^1\f$ matching of the external domain with respect to
	 * the last shell using function like \f$\frac{1}{r^i}\f$ with 
	 * \f${\tt puis} \leq i \leq {\tt puis+nbre}\f$ for each spherical harmonics 
	 * with \f$l \leq {\tt lmax}\f$.
	 */
	void raccord_c1_zec (int puis, int nbre, int lmax) ;
	/**
	 * Matching of the external domain with the outermost shell
	 */
	void raccord_externe (int puis, int nbre, int lmax) ;
};
ostream& operator<<(ostream& , const Cmp & ) ;	

/**
 * \defgroup cmp_m Cmp Mathematics
 * \ingroup (otens)
 * @{
 */

Cmp operator+(const Cmp& ) ;			///< + Cmp
Cmp operator-(const Cmp& ) ;			///< \c - Cmp
Cmp operator+(const Cmp&, const Cmp &) ;	///< Cmp + Cmp
Cmp operator+(const Cmp&, double ) ;		///< Cmp + double
Cmp operator+(double, const Cmp& ) ;		///< double + Cmp 
Cmp operator+(const Cmp&, int ) ;		///< Cmp + int
Cmp operator+(int, const Cmp& ) ;		///< int + Cmp 
Cmp operator-(const Cmp &, const Cmp &) ;	///< Cmp - Cmp
Cmp operator-(const Cmp&, double ) ;		///< Cmp - double
Cmp operator-(double, const Cmp& ) ;		///< double - Cmp 
Cmp operator-(const Cmp&, int ) ;		///< Cmp - int
Cmp operator-(int, const Cmp& ) ;		///< int - Cmp 
Cmp operator*(const Cmp &, const Cmp &) ;	///< Cmp * Cmp
Cmp operator%(const Cmp &, const Cmp &) ;	///< Cmp * Cmp with desaliasing
Cmp operator*(const Cmp&, double ) ;		///< Cmp * double
Cmp operator*(double, const Cmp &) ;		///< double * Cmp
Cmp operator*(const Cmp&, int ) ;		///< Cmp * int
Cmp operator*(int, const Cmp& ) ;		///< int * Cmp 
Cmp operator/(const Cmp &, const Cmp &) ;	///< Cmp / Cmp
Cmp operator/(const Cmp&, double ) ;		///< Cmp / double
Cmp operator/(double, const Cmp &) ;		///< double / Cmp
Cmp operator/(const Cmp&, int ) ;		///< Cmp / int
Cmp operator/(int, const Cmp &) ;		///< int / Cmp

Cmp sin(const Cmp& ) ;		///< Sine
Cmp cos(const Cmp& ) ;		///< Cosine
Cmp tan(const Cmp& ) ;		///< Tangent
Cmp asin(const Cmp& ) ;		///< Arcsine
Cmp acos(const Cmp& ) ;		///< Arccosine
Cmp atan(const Cmp& ) ;		///< Arctangent
Cmp exp(const Cmp& ) ;		///< Exponential
Cmp log(const Cmp& ) ;		///< Neperian logarithm
Cmp log10(const Cmp& ) ;	///< Basis 10 logarithm
Cmp sqrt(const Cmp& ) ;		///< Square root
Cmp racine_cubique (const Cmp& ) ;		///< Cube root
Cmp pow(const Cmp& , int ) ;	///< Power \f${\tt Cmp} ^{\tt int}\f$
Cmp pow(const Cmp& , double ) ; ///< Power \f${\tt Cmp} ^{\tt double}\f$
Cmp abs(const Cmp& ) ;		///< Absolute value

/**
 * Maximum values of a \c Cmp  in each domain.
 * @return 1-D \c Tbl  of size the number of domains, the elements of which 
 *	   are the set of the maximum values in each domain.  
 */
Tbl max(const Cmp& ) ;   

/**
 * Minimum values of a \c Cmp  in each domain.
 * @return 1-D \c Tbl  of size the number of domains, the elements of which 
 *	   are the set of the minimum values in each domain.  
 */
Tbl min(const Cmp& ) ;   

/**
 * Sums of the absolute values of all the values of the \c Cmp  
 * in each domain.
 * @return 1-D \c Tbl  of size the number of domains, the elements of which 
 *	   are the set of the sums of the absolute values in each domain.  
 */
Tbl norme(const Cmp& ) ;   

/**
 * Relative difference between two \c Cmp  (norme version).
 * @return 1-D \c Tbl  of size the number of domains, the elements of which 
 *	   are \c norme[a(l)-b(l)]/norme[b(l)]  if \c b(l)!=0  and
 *	   \c norme[a(l)-b(l)]  if  \c b(l)=0 ,  where \c a(l)  and 
 *	   \c b(l)  denote symbolically the values of \c a  and \c b  
 *	   in domain no. \c l . 
 */
Tbl diffrel(const Cmp& a, const Cmp& b) ; 

/**
 * Relative difference between two \c Cmp  (max version).
 * @return 1-D \c Tbl  of size the number of domains, the elements of which 
 *	   are \c max[abs(a(l)-b(l))]/max[abs(b(l))]  if \c b(l)!=0  and
 *	   \c max[abs(a(l)-b(l))]  if  \c b(l)=0 ,  where \c a(l)  and 
 *	   \c b(l)  denote symbolically the values of \c a  and \c b  
 *	   in domain no. \c l . 
 */
Tbl diffrelmax(const Cmp& a, const Cmp& b) ; 

/** @}*/

}
#endif
