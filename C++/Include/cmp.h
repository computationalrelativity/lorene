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

#include <stdio.h>

#include "valeur.h"
#include "map.h"

class Param ; 

/**
 * Component of a tensorial field.
 * 
 * @version #$Id$#
 * 
 */

class Cmp {

    // Data : 
    // -----
    private:
	const Map* mp ;	    /// Reference mapping

	/// Logical state ({\tt ETATNONDEF}, {\tt ETATQCQ} or {\tt ETATZERO}).
	int etat ;	    

	/**
	 * Power of {\it r} by which the quantity represented by {\tt this} 
	 * must be divided in the external compactified zone in order 
	 * to get the correct physical values
	 */
	int dzpuis ;	

    public:
	Valeur va ;		/// The numerical value of the {\tt Cmp}    

    // Derived data : 
    // ------------
    private:
	/// Pointer on $\partial/\partial r$ of {\tt *this}
	mutable Cmp* p_dsdr ;	
	/// Pointer on $1/r \partial/\partial \theta$ of {\tt *this}
	mutable Cmp* p_srdsdt ;	
	/// Pointer on $1/(r\sin\theta) \partial/\partial \phi$ of {\tt *this}
	mutable Cmp* p_srstdsdp ;
	
	/** Pointer on $\partial/\partial x$ of {\tt *this},
	 *  where $x=r\sin\theta \cos\phi$
	 */
	mutable Cmp* p_dsdx ;	

	/** Pointer on $\partial/\partial y$ of {\tt *this},
	 *  where $y=r\sin\theta \sin\phi$
	 */
	mutable Cmp* p_dsdy ;	

	/** Pointer on $\partial/\partial z$ of {\tt *this},
	 *  where $z=r\cos\theta$
	 */
	mutable Cmp* p_dsdz ;	

	/** Pointer on the Laplacian of {\tt *this}
	 */
	mutable Cmp* p_lap ;	

	/** Power of {\it r} by which the last computed Laplacian has been 
	 *  multiplied in the external compactified domain.  
	 */
	mutable int ind_lap ; 

	/** Pointer on the space integral of {\tt *this} (values in each 
	 *  domain)
	 */
	mutable Tbl* p_integ ; 

    // Constructors - Destructor
    // -------------------------
	
    public:
	explicit Cmp(const Map& ) ;	/// Constructor from mapping
	explicit Cmp(const Map* ) ;	/// Constructor from mapping
	Cmp(const Cmp& ) ;		/// Copy constructor

	/// Constructor from a file (see {\tt sauve(FILE* )})
	Cmp(const Map&, const Mg3d&, FILE* ) ;    		

	~Cmp() ;			/// Destructor

    // Assignment
    // -----------
    public: 
	/// Assignment to another {\tt Cmp} defined on the same mapping
	void operator=(const Cmp&) ;	

	void operator=(const Valeur& ) ; /// Assignment to a {\tt Valeur}
	void operator=(const Mtbl& ) ;	 /// Assignment to a {\tt Mtbl}
	void operator=(double ) ;	 /// Assignment to a {\tt double}
	void operator=(int ) ;		 /// Assignment to an {\tt int}
	    
	/** Assignment to another {\tt Cmp} defined on a different mapping.
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original {\tt Cmp}. 
	 *	@param ci [input] {\tt Cmp} to be imported.
	 */
	void import(const Cmp& ci) ;	 

	/** Assignment to another {\tt Cmp} defined on a different mapping.
	 *  Case where the {\tt Cmp} is symmetric with respect to the plane y=0.
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original {\tt Cmp}. 
	 *	@param ci [input] {\tt Cmp} to be imported.
	 */
	void import_symy(const Cmp& ci) ;	 

	/** Assignment to another {\tt Cmp} defined on a different mapping.
	 *  Case where the {\tt Cmp} is antisymmetric with respect to the 
	 *  plane y=0.
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original {\tt Cmp}. 
	 *	@param ci [input] {\tt Cmp} to be imported.
	 */
	void import_asymy(const Cmp& ci) ;	 

	/** Assignment to another {\tt Cmp} defined on a different mapping.
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original {\tt Cmp}. 
	 *	@param nzet [input] Number of domains of the destination
	 *			    mapping (i.e. {\tt this->mp}) where the 
	 *			    importation is performed: the assignment
	 *			    is done for the domains whose indices are
	 *			    between 0 and {\tt nzet-1}. In the other
	 *			    domains, {\tt *this} is set to zero. 
	 *	@param ci [input] {\tt Cmp} to be imported.
	 */
	void import(int nzet, const Cmp& ci) ;	 

	/** Assignment to another {\tt Cmp} defined on a different mapping.
	 *  Case where the {\tt Cmp} is symmetric with respect to the plane y=0.
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original {\tt Cmp}. 
	 *	@param nzet [input] Number of domains of the destination
	 *			    mapping (i.e. {\tt this->mp}) where the 
	 *			    importation is performed: the assignment
	 *			    is done for the domains whose indices are
	 *			    between 0 and {\tt nzet-1}. In the other
	 *			    domains, {\tt *this} is set to zero. 
	 *	@param ci [input] {\tt Cmp} to be imported.
	 */
	void import_symy(int nzet, const Cmp& ci) ;	 

	/** Assignment to another {\tt Cmp} defined on a different mapping.
	 *  Case where the {\tt Cmp} is antisymmetric with respect to the 
	 *  plane y=0.
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original {\tt Cmp}. 
	 *	@param nzet [input] Number of domains of the destination
	 *			    mapping (i.e. {\tt this->mp}) where the 
	 *			    importation is performed: the assignment
	 *			    is done for the domains whose indices are
	 *			    between 0 and {\tt nzet-1}. In the other
	 *			    domains, {\tt *this} is set to zero. 
	 *	@param ci [input] {\tt Cmp} to be imported.
	 */
	void import_asymy(int nzet, const Cmp& ci) ;	 

    private:
	/** Assignment to another {\tt Cmp} defined on a different mapping,
	 *  when the two mappings do not have a particular relative orientation.
	 *
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original {\tt Cmp}. 
	 *	@param nzet [input] Number of domains of the destination
	 *			    mapping (i.e. {\tt this->mp}) where the 
	 *			    importation is performed: the assignment
	 *			    is done for the domains whose indices are
	 *			    between 0 and {\tt nzet-1}. In the other
	 *			    domains, {\tt *this} is set to zero. 
	 *	@param ci [input] {\tt Cmp} to be imported.
	 */
	void import_gal(int nzet, const Cmp& ci) ;	 

	/** Assignment to another {\tt Cmp} defined on a different mapping,
	 *  when the two mappings have aligned Cartesian axis. 
	 *
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original {\tt Cmp}. 
	 *	@param nzet [input] Number of domains of the destination
	 *			    mapping (i.e. {\tt this->mp}) where the 
	 *			    importation is performed: the assignment
	 *			    is done for the domains whose indices are
	 *			    between 0 and {\tt nzet-1}. In the other
	 *			    domains, {\tt *this} is set to zero. 
	 *	@param ci [input] {\tt Cmp} to be imported.
	 */
	void import_align(int nzet, const Cmp& ci) ;	 

	/** Assignment to another {\tt Cmp} defined on a different mapping,
	 *  when the two mappings have anti-aligned Cartesian axis (i.e.
	 *  $x_1 = - x_2$,  $y_1 = - y_2$,  $z_1 = z_2$). 
	 *
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original {\tt Cmp}. 
	 *	@param nzet [input] Number of domains of the destination
	 *			    mapping (i.e. {\tt this->mp}) where the 
	 *			    importation is performed: the assignment
	 *			    is done for the domains whose indices are
	 *			    between 0 and {\tt nzet-1}. In the other
	 *			    domains, {\tt *this} is set to zero. 
	 *	@param ci [input] {\tt Cmp} to be imported.
	 */
	void import_anti(int nzet, const Cmp& ci) ;	 

	/** Assignment to another {\tt Cmp} defined on a different mapping,
	 *  when the two mappings have aligned Cartesian axis. 
	 *  Case where the {\tt Cmp} is symmetric with respect to the plane y=0.
	 *
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original {\tt Cmp}. 
	 *	@param nzet [input] Number of domains of the destination
	 *			    mapping (i.e. {\tt this->mp}) where the 
	 *			    importation is performed: the assignment
	 *			    is done for the domains whose indices are
	 *			    between 0 and {\tt nzet-1}. In the other
	 *			    domains, {\tt *this} is set to zero. 
	 *	@param ci [input] {\tt Cmp} to be imported.
	 */
	void import_align_symy(int nzet, const Cmp& ci) ;	 

	/** Assignment to another {\tt Cmp} defined on a different mapping,
	 *  when the two mappings have anti-aligned Cartesian axis (i.e.
	 *  $x_1 = - x_2$,  $y_1 = - y_2$,  $z_1 = z_2$). 
	 *  Case where the {\tt Cmp} is symmetric with respect to the plane y=0.
	 *
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original {\tt Cmp}. 
	 *	@param nzet [input] Number of domains of the destination
	 *			    mapping (i.e. {\tt this->mp}) where the 
	 *			    importation is performed: the assignment
	 *			    is done for the domains whose indices are
	 *			    between 0 and {\tt nzet-1}. In the other
	 *			    domains, {\tt *this} is set to zero. 
	 *	@param ci [input] {\tt Cmp} to be imported.
	 */
	void import_anti_symy(int nzet, const Cmp& ci) ;	 

	/** Assignment to another {\tt Cmp} defined on a different mapping,
	 *  when the two mappings have aligned Cartesian axis. 
	 *  Case where the {\tt Cmp} is antisymmetric with respect to the 
	 *  plane y=0.
	 *
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original {\tt Cmp}. 
	 *	@param nzet [input] Number of domains of the destination
	 *			    mapping (i.e. {\tt this->mp}) where the 
	 *			    importation is performed: the assignment
	 *			    is done for the domains whose indices are
	 *			    between 0 and {\tt nzet-1}. In the other
	 *			    domains, {\tt *this} is set to zero. 
	 *	@param ci [input] {\tt Cmp} to be imported.
	 */
	void import_align_asymy(int nzet, const Cmp& ci) ;	 

	/** Assignment to another {\tt Cmp} defined on a different mapping,
	 *  when the two mappings have anti-aligned Cartesian axis (i.e.
	 *  $x_1 = - x_2$,  $y_1 = - y_2$,  $z_1 = z_2$). 
	 *  Case where the {\tt Cmp} is antisymmetric with respect to the 
	 *  plane y=0.
	 *
	 *  This assignment is performed point to point by means of the
	 *  spectral expansion of the original {\tt Cmp}. 
	 *	@param nzet [input] Number of domains of the destination
	 *			    mapping (i.e. {\tt this->mp}) where the 
	 *			    importation is performed: the assignment
	 *			    is done for the domains whose indices are
	 *			    between 0 and {\tt nzet-1}. In the other
	 *			    domains, {\tt *this} is set to zero. 
	 *	@param ci [input] {\tt Cmp} to be imported.
	 */
	void import_anti_asymy(int nzet, const Cmp& ci) ;	 



    // Access to individual elements
    // -----------------------------
    public:

	/** Read/write of the value in a given domain.
	 * NB: to gain in efficiency, the method {\tt del\_deriv()} (to delete
	 *     the derived members) is not called by this function. It must
	 *     thus be invoqued by the user.  
	 *
	 * @param l [input] domain index
	 * @return Tbl containing the value of the field in domain {\tt l}.
	 */ 
	Tbl& set(int l) {
	    assert(etat == ETATQCQ) ;
	    return va.set(l) ;
	};
	
	/** Read-only of the value in a given domain.
	 * @param l [input] domain index
	 * @return Tbl containing the value of the field in domain {\tt l}.
	 */ 
	const Tbl& operator()(int l) const {
	    assert(etat == ETATQCQ) ;
	    return va(l) ;
	};


	/** Read/write of a particular element.
	 * NB: to gain in efficiency, the method {\tt del\_deriv()} (to delete
	 *     the derived members) is not called by this function. It must
	 *     thus be invoqued by the user.  
	 *     
	 * @param l [input] domain index
	 * @param k [input] $\phi$ index
	 * @param j [input] $\theta$ index
	 * @param i [input] {\it r} ($\xi$) index
	 */ 
	double& set(int l, int k, int j, int i) {
	    assert(etat == ETATQCQ) ;
	    return va.set(l, k, j, i) ;
	};
	
	
	/** Read-only of a particular element.
	 * @param l [input] domain index
	 * @param k [input] $\phi$ index
	 * @param j [input] $\theta$ index
	 * @param i [input] {\it r} ($\xi$) index
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

	/** Computes the value of the field represented by {\tt *this} at an
	*   arbitrary point $(r, \theta, \phi)$, by means of the spectral 
	*   expansion.
	*	 @param r [input] value of the coordinate {\it r}
	*	 @param theta [input] value of the coordinate $\theta$
	*	 @param phi [input] value of the coordinate $\phi$
	*	 @return value at the point $(r, \theta, \phi)$ 
	*		 of the field represented by {\tt *this}. 
	*/
	double val_point(double r, double theta, double phi) const ; 


    // Memory management
    // -----------------
    private:
	void del_t() ;		    /// Logical destructor
	void del_deriv() ;	    /// Logical destructor of the derivatives
	void set_der_0x0() ;	    /// Sets the pointers for derivatives to 0x0

    public:

    /**
     * Sets the logical state to {\tt ETATNONDEF} (undefined). 
     * Calls the logical destructor of the {\tt Valeur va} and
     * deallocates the memory occupied by all the derivatives. 
     */
	void set_etat_nondef() ;   

    /**
     * Sets the logical state to {\tt ETATZERO} (zero). 
     * Calls the logical destructor of the {\tt Valeur va} and
     * deallocates the memory occupied by all the derivatives. 
     */
	void set_etat_zero() ;	    
	
    /**
     * Sets the logical state to {\tt ETATQCQ} (ordinary state).
     * If the state is already {\tt ETATQCQ}, this function does nothing.
     * Otherwise, it calls the logical destructor of the {\tt Valeur va} and
     * deallocates the memory occupied by all the derivatives.
     */
	void set_etat_qcq() ;	    

    /**
     * Sets the logical state to {\tt ETATQCQ} (ordinary state)
     *  and performs the memory allocation of all the 
     *  elements, down to the {\tt double} arrays of the {\tt Tbl}s. 
     *  This function performs in fact recursive calls to {\tt set\_etat\_qcq()}
     *  on each element of the chain {\tt Cmp} ->
     *  {\tt Valeur} -> {\tt Mtbl} -> {\tt Tbl}. 
     */
	void allocate_all() ; 

    /**
     * Sets the {\tt Cmp} to zero in a hard way. 
     * 1/ Sets the logical state to {\tt ETATQCQ}, i.e. to an ordinary state.
     * 2/ Fills the {\tt Valeur va} with zeros. 
     * NB: this function must be used for debugging purposes only.
     * For other operations, the functions {\tt set\_etat\_zero()}
     * or {\tt annule(int, int)} must be perferred. 
     */
	void annule_hard() ;

    /**
     * Sets the {\tt Cmp} to zero in a given domain.
     *	@param l [input]  Index of the domain in which the {\tt Cmp}
     *			  will be set (logically) to zero.
     */
	void annule(int l) ; 

    /**
     * Sets the {\tt Cmp} to zero in several domains.
     *	@param l_min [input] The {\tt Cmp} will be set (logically) to zero
     *			     in the domains whose indices are in the range
     *			     {\tt [l\_min, l\_max]}.
     *	@param l_max [input] see the comments for {\tt l\_min}.
     * 
     * Note that {\tt annule(0, va.mg->get\_nzone()-1)} is equivalent to
     *	 {\tt set\_etat\_zero()}.
     */
	void annule(int l_min, int l_max) ; 
    
    /**
     * Sets the {\tt n} lasts coefficients in {\it r} to 0 in the external domain.
     */
	void filtre (int n) ;
    
    /**
     * Sets the {\tt n} lasts coefficients in $\Phi$ to 0 in the 
     * domain {\tt zone}.
     */
	void filtre_phi (int n, int zone) ;
    
    /**
     * Sets the value of the {\tt Cmp} to {\tt val} at infinity. This is usefull
     * for dealing with undefined values. The external domain must be 
     * compactified.
     */
	void set_val_inf (double val) ;
    
    /**
     * Sets the value of the {\tt Cmp} to {\tt val} on the inner boudary of the
     * shell number {\tt zone}.This is usefull
     * for dealing with undefined values.
     */
	void set_val_hor (double val, int zone) ;
    /**
     * Substracts all the components behaving like $r^{-n}$ in the external 
     * domain, with {\it n} strictly lower than {\tt puis}, so that {\tt *this} 
     * decreases at least like $r^{\tt puis}$ at infinity.
     */
	void fixe_decroissance (int puis) ;
	
    // Extraction of information
    // -------------------------
    public:
	int get_etat() const {return etat;} ;	/// Returns the logical state
	const Map* get_mp() const {return mp;};	 /// Returns the mapping
	int get_dzpuis() const {return dzpuis;} ; /// Returns {\tt dzpuis}
	
	/** Returns {\tt true} if the last domain is compactified and
	 *  {\tt *this} is not zero in this domain
	 */
	bool dz_nonzero() const ; 
	
	/** Returns {\tt false} if the last domain is compactified 
	 *  and {\tt *this} is not zero in this domain and {\tt dzpuis}
	 *  is not equal to {\tt dzi}, otherwise return true. 
	 */
	bool check_dzpuis(int dzi) const ; 
	
	
    // Outputs
    // -------
    public:
	void sauve(FILE *) const ;	    /// Save in a file
    
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
	void operator+=(const Cmp &) ;		    /// += Cmp
	void operator-=(const Cmp &) ;		    /// -= Cmp
	void operator*=(const Cmp &) ;		    /// *= Cmp

    // Manipulation of spectral bases
    // ------------------------------    
    /** Sets the spectral bases of the {\tt Valeur va} to the standard ones 
     *  for a scalar
     */
    void std_base_scal() ;	 


    // Differential operators and others
    // ---------------------------------
    public:
	/** Returns $\partial / \partial r$ of {\tt *this}.
	 *  Note that in the external compactified domain (ZEC), it returns
	 *  instead $r^2 \partial/ \partial r$.
	 */
	const Cmp& dsdr() const ; 
	
	/** Returns $1/r \partial / \partial \theta$ of {\tt *this}.
	 *  Note that in the external compactified domain (ZEC), it returns
	 *  instead $r \partial/ \partial \theta$.
	 */
	const Cmp& srdsdt() const ; 

	/** Returns $1/(r\sin\theta) \partial / \partial \phi$ of {\tt *this}.
	 *  Note that in the external compactified domain (ZEC), it returns
	 *  instead $r/\sin\theta \partial/ \partial \phi$.
	 */
	const Cmp& srstdsdp() const ; 

	/** Returns $\partial/\partial x$ of {\tt *this},
	 *  where $x=r\sin\theta \cos\phi$.
	 *  Note that in the external compactified domain (ZEC), it returns
	 *  instead $r^2 \partial/ \partial x$.
	 */
	const Cmp& dsdx() const ;	

	/** Returns $\partial/\partial y$ of {\tt *this},
	 *  where $y=r\sin\theta \sin\phi$.
	 *  Note that in the external compactified domain (ZEC), it returns
	 *  instead $r^2 \partial/ \partial y$.
	 */
	const Cmp& dsdy() const ;	

	/** Returns $\partial/\partial z$ of {\tt *this},
	 *  where $z=r\cos\theta$.
	 *  Note that in the external compactified domain (ZEC), it returns
	 *  instead $r^2 \partial/ \partial z$.
	 */
	const Cmp& dsdz() const ;	

	/** Returns $\partial/\partial x_i$ of {\tt *this},
	 *  where $x_i = (x, y, z)$.
	 *  Note that in the external compactified domain (ZEC), it returns
	 *  instead $r^2 \partial/ \partial x_i$.
	 *  @param i [input] i=0 for {\it x},  i=1 for {\it y}, i=2 for {\it z}.
	 */
	const Cmp& deriv(int i) const ;	

	/** Returns the Laplacian of {\tt *this}
	 *   @param zec_mult_r [input] Determines the quantity computed in
	 *			 the external compactified domain (ZEC) 
	 *		({\it u} in the field represented by {\tt *this}) :  \\
	 *		    zec\_mult\_r = 0 : $\Delta u$	\\
	 *		    zec\_mult\_r = 2 : $r^2 \,  \Delta u$	\\
	 *		    zec\_mult\_r = 4 (default) : $r^4 \, \Delta u$	
	 */
	const Cmp& laplacien(int zec_mult_r = 4) const ; 

	void div_r() ;    /// Division by {\it r} everywhere.

	void mult_r() ;   /// Multiplication by {\it r} everywhere.

	/** Multiplication by {\it r} in the external compactified domain (ZEC)
	 */
	void mult_r_zec() ;
	
	void mult_rsint() ;   /// Multiplication by $r\sin\theta$

	void div_rsint() ;    /// Division by $r\sin\theta$

	/** Decreases by 1 the value of {\tt dzpuis} and changes accordingly
	 *  the values of the {\tt Cmp} in the external compactified domain (ZEC).
	 */
	void dec_dzpuis() ; 

	/** Increases by the value of {\tt dzpuis} and changes accordingly
	 *  the values of the {\tt Cmp} in the external compactified domain (ZEC).
	 */
	void inc_dzpuis() ; 
	
	/** Decreases by 2 the value of {\tt dzpuis} and changes accordingly
	 *  the values of the {\tt Cmp} in the external compactified domain (ZEC).
	 */
	void dec2_dzpuis() ; 

	/** Increases by 2 the value of {\tt dzpuis} and changes accordingly
	 *  the values of the {\tt Cmp} in the external compactified domain (ZEC).
	 */
	void inc2_dzpuis() ; 

	void set_dzpuis(int ) ;  /// Set a value to {\tt dzpuis}

	/** Computes the integral over all space of {\tt *this}.
	 *  The computed quantity is ({\it u} being the field represented by
	 *   {\tt *this})
	 *    $\int u \, r^2 \sin\theta \,  dr\, d\theta \, d\phi$.
	 *  Note that in the external compactified domain (ZEC), {\tt dzpuis} 
	 *  must be 4 for the computation to take place. 
	 */
	double integrale() const ; 
	
	/** Computes the integral in each domain of {\tt *this}.
	 *  The computed quantity is ({\it u} being the field represented by
	 *   {\tt *this})
	 *    $\int u \, r^2 \sin\theta \,  dr\, d\theta \, d\phi$
	 *  in each domain. The result is returned a {\tt Tbl} on the 
	 *  various domains. 
	 *  Note that in the external compactified domain (ZEC), {\tt dzpuis} 
	 *  must be 4 for the computation to take place. 
	 */
	const Tbl& integrale_domains() const ; 
	
	/** Asymptotic expansion at r = infinity. 
	 * 
	 *  Determines the coefficients $a_k(\theta, \phi)$ of the expansion
	 *  \begin{equation}
	 *	\sum_{k=0}^n {a_k(\theta, \phi) \over r^k}
	 *  \end{equation} 
	 *  of {\tt *this} when $r \rightarrow \infty$. 
	 *
	 *	@param n order of the expansion
	 *	@param flag : output
	 *	@return Array of {\tt n}+1 {\tt Valeur}s on {\tt mg->angu} 
	 *		describing the coefficients $a_k(\theta, \phi)$. 
	 *		This array is allocated by the routine. 
	 * 
	 */
	Valeur** asymptot(int n, const int flag = 0) const ; 
	
	/// Function to compare the values of two Cmp

	void compare(FILE* fich, const char* name_i) ;
	void compare(const Cmp& comp, const char* name, int ii = -1
		     , int jj = -1) ;
 



    // PDE resolution 
    // --------------
    public:
	/** Solves the scalar Poisson equation with {\tt *this} as a source.
	 *   The source $\sigma$ of the equation $\Delta u = \sigma$ is 
	 *   represented by the {\tt Cmp} {\tt *this}. 
	 *   Note that {\tt dzpuis} must be equal to 2 or 4, i.e. that the
	 *   quantity stored in {\tt *this} is in fact $r^2 \sigma$ or
	 *   $r^4 \sigma$ in the external compactified domain. 
	 *   The solution {\it u} with the boundary condition {\it u}=0 at spatial
	 *   infinity is the returned {\tt Cmp}. 
	 */
	Cmp poisson() const ;

	/** Solves the scalar Poisson equation with {\tt *this} as a source
	 *   (version with parameters to control the resolution).
	 *   The source $\sigma$ of the equation $\Delta u = \sigma$ is 
	 *   represented by the {\tt Cmp} {\tt *this}. 
	 *   Note that {\tt dzpuis} must be equal to 2 or 4, i.e. that the
	 *   quantity stored in {\tt *this} is in fact $r^2 \sigma$ or
	 *   $r^4 \sigma$ in the external compactified domain. 
	 *   @param par [input/output] possible parameters
	 *   @param uu [input/output] solution {\it u} with the boundary condition 
	 *   {\it u}=0 at spatial infinity. 
	 */
	void poisson(Param& par, Cmp& uu) const ;
	
	/**
	 * Is identicall to {\tt Cmp::poisson()}. The regularity condition at the 
	 * origin is replace by a boundary condition of the Dirichlet type.
	 * 
	 * @param limite [input] : angular function. The boundary condition is 
	 * given by {\tt limite[num]}.
	 * @param num [input] : index of the boudary at which the condition is to 
	 * be fullfilled.
	 * 
	 * More precisely we impose the solution is equal to {\tt limite[num]} at the
	 * boundary between the domains {\tt num} and {\tt num+1} (the latter one being 
	 * a shell).
	 * 
	 */
	Cmp poisson_dirichlet (const Valeur& limite, int num) const ;
	
	/**
	 * Idem as {\tt Cmp::poisson\_neumann}, the boundary condition being on 
	 * the radial derivative of the solution.
	 */
	Cmp poisson_neumann   (const Valeur&, int) const ;
	Cmp poisson_frontiere_double   (const Valeur&, const Valeur&, int) const ;

	/** Solves the scalar Poisson equation with {\tt *this} as a source
	 *   (version with parameters to control the resolution).
	 *   The source $\sigma$ of the equation $\Delta u = \sigma$ is 
	 *   represented by the {\tt Cmp} {\tt *this}. 
	 *   The regularized source
	 *   $\sigma_{\rm regu} = \sigma - \sigma_{\rm div}$
	 *   is constructed and solved.
	 *   Note that {\tt dzpuis} must be equal to 2 or 4, i.e. that the
	 *   quantity stored in {\tt *this} is in fact $r^2 \sigma$ or
	 *   $r^4 \sigma$ in the external compactified domain.
	 *   @param k_div [input] regularization degree of the procedure
	 *   @param nzet [input] number of domains covering the star
	 *   @param unsgam1 [input] parameter $1/(\gamma-1)$ where $\gamma$
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

	/** Checks if a Poisson equation with {\tt *this} as a source
	 *  has been correctly solved.
	 * 
	 *  @param uu [input] Solution {\it u} of the Poisson equation
	 *		      $\Delta u = \sigma$,  $\sigma$ being 
	 *		      represented by the {\tt Cmp} {\tt *this}.
	 * 
	 *  @param ostr [input/output] Output stream used for displaying
	 *		{\tt err}.
	 *
	 *  @param detail [input] if {\tt true} displays {\tt err(0, *)}, 
	 *		    {\tt err(1, *)} and {\tt err(2, *)} \\
	 *		if {\tt false} (default),  displays only 
	 *		the relative error {\tt err(0, *)}. 
	 *  
	 *  @return 2-D {\tt Tbl} {\tt err} decribing the errors in each 
	 *	    domain: \\
	 *	{\tt err(0, l) : } Relative error in domain no. {\tt l}, 
	 *	    defined as the maximum value of 
	 *	    $|\Delta u - \sigma|$ in that domain divided by {\it M}, 
	 *	    where {\it M} is the maximum value of $|\sigma|$ 
	 *	    over all domains if {\tt dzpuis = 0} or $\sigma$ is
	 *	    zero in the external compactified domain (ECD). If 
	 *	    {\tt dzpuis != 0} and $\sigma$ does not vanish in the 
	 *	    ECD, the value of {\it M} used in the
	 *	    non-compactified domains is the maximum value over
	 *	    these domains, whereas the value of {\it M} used in the
	 *	    external compactified domain is the maximum value
	 *	    on that particular domain. \\
	 *	{\tt err(1, l) : }  Maximum value of the absolute error
	 *			$|\Delta u - \sigma|$ in domain no. {\tt l} \\
	 *	{\tt err(2, l) : }  Maximum value of $|\sigma|$ in domain 
	 *			    no. {\tt l} 
	 */
	Tbl test_poisson(const Cmp& uu, ostream& ostr, 
					bool detail = false) const ;  

	/** Performs one time-step integration (from $t=j \to j+1$) of the 
	 *   scalar d'Alembert equation with {\tt *this} being the value of 
	 *   the function {\it f} at time {\it j}.
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
	 *   {\tt Param}, so one should not modify it.\\
	 *   There should be one {\tt Param} for each wave equation to be 
	 *   solved. 
	 *   @param fJm1 [input] solution $f^{J-1}$ at time {\it J-1}
	 *   @param source [input] source $\sigma$ of the d'Alembert equation 
	 *	    $\diamond u = \sigma$.
	 *   @return solution $f^{J+1}$ at time {\it J+1}
	 *   with boundary conditions defined by {\tt par.get\_int(0)}.
	 */
	Cmp avance_dalembert(Param& par, const Cmp& fjm1, const Cmp& source) 
	  const ;
	
	/**
	 * Performs the $C^n$ matching of the nucleus with respect to the 
	 * first shell.
	 */
	void raccord(int n) ;
	
	/**
	 * Performs the $C^1$ matching of the external domain with respect to
	 * the last shell using function like $\frac{1}{r^i}$ with 
	 * ${\tt puis} \leq i \leq {\tt puis+nbre}$ for each spherical harmonics 
	 * with $l \leq {\tt lmax}$.
	 */
	void raccord_c1_zec (int puis, int nbre, int lmax) ;
	/**
	 * Matching of the external domain with the outermost shell
	 */
	void raccord_externe (int puis, int nbre, int lmax) ;
};
ostream& operator<<(ostream& , const Cmp & ) ;	

// Prototypage de l'arithmetique
/**
 * @name Cmp mathematics
 */
    //@{
Cmp operator+(const Cmp& ) ;			/// + Cmp
Cmp operator-(const Cmp& ) ;			/// - Cmp
Cmp operator+(const Cmp&, const Cmp &) ;	/// Cmp + Cmp
Cmp operator+(const Cmp&, double ) ;		/// Cmp + double
Cmp operator+(double, const Cmp& ) ;		/// double + Cmp 
Cmp operator+(const Cmp&, int ) ;		/// Cmp + int
Cmp operator+(int, const Cmp& ) ;		/// int + Cmp 
Cmp operator-(const Cmp &, const Cmp &) ;	/// Cmp - Cmp
Cmp operator-(const Cmp&, double ) ;		/// Cmp - double
Cmp operator-(double, const Cmp& ) ;		/// double - Cmp 
Cmp operator-(const Cmp&, int ) ;		/// Cmp - int
Cmp operator-(int, const Cmp& ) ;		/// int - Cmp 
Cmp operator*(const Cmp &, const Cmp &) ;	/// Cmp * Cmp
Cmp operator%(const Cmp &, const Cmp &) ;	/// Cmp * Cmp with desaliasing
Cmp operator*(const Cmp&, double ) ;		/// Cmp * double
Cmp operator*(double, const Cmp &) ;		/// double * Cmp
Cmp operator*(const Cmp&, int ) ;		/// Cmp * int
Cmp operator*(int, const Cmp& ) ;		/// int * Cmp 
Cmp operator/(const Cmp &, const Cmp &) ;	/// Cmp / Cmp
Cmp operator/(const Cmp&, double ) ;		/// Cmp / double
Cmp operator/(double, const Cmp &) ;		/// double / Cmp
Cmp operator/(const Cmp&, int ) ;		/// Cmp / int
Cmp operator/(int, const Cmp &) ;		/// int / Cmp

Cmp sin(const Cmp& ) ;		/// Sine
Cmp cos(const Cmp& ) ;		/// Cosine
Cmp tan(const Cmp& ) ;		/// Tangent
Cmp asin(const Cmp& ) ;		/// Arcsine
Cmp acos(const Cmp& ) ;		/// Arccosine
Cmp atan(const Cmp& ) ;		/// Arctangent
Cmp exp(const Cmp& ) ;		/// Exponential
Cmp log(const Cmp& ) ;		/// Neperian logarithm
Cmp log10(const Cmp& ) ;	/// Basis 10 logarithm
Cmp sqrt(const Cmp& ) ;		/// Square root
Cmp racine_cubique (const Cmp& ) ;		/// Cube root
Cmp pow(const Cmp& , int ) ;	/// Power ${\tt Cmp}^{\tt int}$
Cmp pow(const Cmp& , double ) ; /// Power ${\tt Cmp}^{\tt double}$
Cmp abs(const Cmp& ) ;		/// Absolute value

/**
 * Maximum values of a {\tt Cmp} in each domain.
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are the set of the maximum values in each domain.  
 */
Tbl max(const Cmp& ) ;   

/**
 * Minimum values of a {\tt Cmp} in each domain.
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are the set of the minimum values in each domain.  
 */
Tbl min(const Cmp& ) ;   

/**
 * Sums of the absolute values of all the values of the {\tt Cmp} 
 * in each domain.
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are the set of the sums of the absolute values in each domain.  
 */
Tbl norme(const Cmp& ) ;   

/**
 * Relative difference between two {\tt Cmp} (norme version).
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are {\tt norme[a(l)-b(l)]/norme[b(l)]} if {\tt b(l)!=0} and
 *	   {\tt norme[a(l)-b(l)]} if  {\tt b(l)=0},  where {\tt a(l)} and 
 *	   {\tt b(l)} denote symbolically the values of {\tt a} and {\tt b} 
 *	   in domain no. {\tt l}. 
 */
Tbl diffrel(const Cmp& a, const Cmp& b) ; 

/**
 * Relative difference between two {\tt Cmp} (max version).
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are {\tt max[abs(a(l)-b(l))]/max[abs(b(l))]} if {\tt b(l)!=0} and
 *	   {\tt max[abs(a(l)-b(l))]} if  {\tt b(l)=0},  where {\tt a(l)} and 
 *	   {\tt b(l)} denote symbolically the values of {\tt a} and {\tt b} 
 *	   in domain no. {\tt l}. 
 */
Tbl diffrelmax(const Cmp& a, const Cmp& b) ; 

    //@}
#endif
