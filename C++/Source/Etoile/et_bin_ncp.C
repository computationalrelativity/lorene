/*
 *  Methods for class Et_bin_ncp
 *
 */

/*
 *   Copyright (c) 2002  Francois Limousin
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


char et_bin_ncp_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2003/02/04 16:55:39  f_limousin
 * Add several members and computational routines
 *
 * Revision 1.4  2003/01/20 17:13:26  j_novak
 * Modif des include <math.h> pour eviter les warning sous SGI.
 *
 * Revision 1.3  2003/01/20 17:07:05  f_limousin
 * add <math.h>
 *
 * Revision 1.2  2003/01/17 13:38:29  f_limousin
 * Add computational routines
 *
 * Revision 1.1  2002/12/09 10:46:50  f_limousin
 * Methods for class Et_bin_ncp.
 *
 *
 *
 *
 * $Header$
 *
 */

#include <math.h>

// Lorene headers
#include "et_bin_ncp.h"
#include "unites.h"

			    //--------------//
			    // Constructors //
			    //--------------//

// Standard constructor
// --------------------

Et_bin_ncp::Et_bin_ncp(Map& mp_i, int nzet_i, bool relat, const Eos& eos_i,
		       bool irrot, const Base_vect& ref_triad_i, const Metrique& flat0, const Tenseur_sym& source) 
             : Etoile_bin(mp_i, nzet_i, relat, eos_i, irrot, ref_triad_i),
               met_gamma(source),
	       flat(flat0),
               gtilde(mp_i, met_gamma, flat0),
	       gamma(met_gamma.determinant()),
	       loggamma(log(gamma)),
	       loggamma_auto(mp_i),
	       loggamma_comp(mp_i),
	       metgamma_auto(source),
	       metgamma_comp(source),
	       gtilde_auto(mp_i, met_gamma, flat0),
	       gtilde_comp(mp_i, met_gamma, flat0),
	       kcar_auto(mp_i),
	       kcar_comp(mp_i) {



  // All quantities are initialized to zero : 
  kcar_auto = 0 ;
  kcar_comp = 0 ;

  // The metric is initialized to the flat one
  gamma = 1 ; 
  loggamma = 0 ;
  loggamma_auto = 0 ;
  loggamma_comp = 0 ;

  gtilde = flat.cov() ;
  gtilde_auto = 0.5*flat.cov() ;
  gtilde_comp = 0.5*flat.cov() ; 
}

// Copy constructor
// ----------------

Et_bin_ncp::Et_bin_ncp(const Et_bin_ncp& et)
      	   : Etoile_bin(et),
	     met_gamma(et.met_gamma),
	     flat(et.flat),
             gtilde(et.gtilde),	     
	     gamma(et.gamma),
	     loggamma(et.loggamma),
	     loggamma_auto(et.loggamma_auto),
	     loggamma_comp(et.loggamma_comp),
	     metgamma_auto(et.metgamma_auto),
	     metgamma_comp(et.metgamma_comp),
	     gtilde_auto(et.gtilde_auto),
	     gtilde_comp(et.gtilde_comp),
	     kcar_auto(et.kcar_auto),
	     kcar_comp(et.kcar_comp) {}


// Constructor from a file
// -----------------------
Et_bin_ncp::Et_bin_ncp(Map& mp_i, const Eos& eos_i, const Base_vect& ref_triad_i,
		       const Metrique& flat0, FILE* fich)
  : Etoile_bin(mp_i, eos_i, ref_triad_i, fich),
    met_gamma(mp_i),
    flat(flat0),
    gtilde(mp_i, met_gamma, flat0),
    gamma(met_gamma.determinant()),
    loggamma(log(gamma)),
    loggamma_auto(mp_i),
    loggamma_comp(mp_i),
    metgamma_auto(mp_i),
    metgamma_comp(mp_i),
    gtilde_auto(mp_i, met_gamma, flat0),
    gtilde_comp(mp_i, met_gamma, flat0),
    kcar_auto(mp_i),
    kcar_comp(mp_i) {

}
  


			    //------------//
			    // Destructor //
			    //------------//

Et_bin_ncp::~Et_bin_ncp(){}



			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Et_bin_ncp
// --------------------------------
void Et_bin_ncp::operator=(const Et_bin_ncp& et) {

// Assignement of proper quantities of class Etoile_bin
met_gamma = et.met_gamma ; 
flat = et.flat ;
gtilde = et.gtilde ;
gamma = et.gamma ;
loggamma = et.loggamma ;
loggamma_auto = et.loggamma_auto ;
loggamma_comp = et.loggamma_comp ;
metgamma_auto = et.metgamma_auto ;
metgamma_comp = et.metgamma_comp ;
gtilde_auto = et.gtilde_auto ;
gtilde_comp = et.gtilde_comp ;
kcar_auto = et.kcar_auto ;
kcar_comp = et.kcar_comp ;

}



			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Et_bin_ncp::sauve(FILE* fich) const {
    
    Etoile_bin::sauve(fich) ; 
    
    met_gamma.sauve(fich) ; 
    flat.sauve(fich) ; 
    gtilde.sauve(fich) ; 
    gamma.sauve(fich) ;
    loggamma.sauve(fich) ;
    loggamma_auto.sauve(fich) ;
    loggamma_comp.sauve(fich) ;
    metgamma_auto.sauve(fich) ;
    metgamma_comp.sauve(fich) ;
    gtilde_auto.sauve(fich) ;
    gtilde_comp.sauve(fich) ;
    kcar_auto.sauve(fich) ;
    kcar_comp.sauve(fich) ;
}

// Printing
// --------

ostream& Et_bin_ncp::operator>>(ostream& ost) const {
    
    Etoile_bin::operator>>(ost) ; 
    
    ost << endl ; 
    ost << "Star in a binary system with non conformally flat metric" 
	<< endl ; 
    ost << "--------------------------------------------------------" 
	<< endl ; 

    ost << "metrique Gamma : " << met_gamma << endl ; 
    ost << "Flat : " << flat << endl ; 
    ost << "Gamma tilde : " << gtilde << endl ;
    ost << "Gamma tilde auto : " << gtilde_auto << endl ;
    ost << "Gamma tilde comp : " << gtilde_comp << endl ;
    ost << "Determinant gamma : " << gamma << endl ;
    ost << "log(determinant) : " << loggamma << endl ;
    ost << "log(determinant) auto : " << loggamma_auto << endl ;
    ost << "log(determinant) comp : " << loggamma_comp << endl ;

   ost << endl << "Central K^{ij} [c/km] : " << endl ; 
    ost << "  K^{xx} auto, comp : " 
	<< tkij_auto(0, 0)(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(0, 0)(0, 0, 0, 0) * km << endl ; 
    ost << "  K^{xy} auto, comp : " 
	<< tkij_auto(0, 1)(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(0, 1)(0, 0, 0, 0) * km << endl ; 
    ost << "  K^{xz} auto, comp : " 
	<< tkij_auto(0, 2)(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(0, 2)(0, 0, 0, 0) * km << endl ; 
    ost << "  K^{yy} auto, comp : " 
	<< tkij_auto(1, 1)(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(1, 1)(0, 0, 0, 0) * km << endl ; 
    ost << "  K^{yz} auto, comp : " 
	<< tkij_auto(1, 2)(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(1, 2)(0, 0, 0, 0) * km << endl ; 
    ost << "  K^{zz} auto, comp : " 
	<< tkij_auto(2, 2)(0, 0, 0, 0) * km  << "  "
	<< tkij_comp(2, 2)(0, 0, 0, 0) * km << endl ; 

    ost << endl << "Central K_{ij} K^{ij} [c^2/km^2] : " << endl ; 
    ost << "   K_{ij} K^{ij}  auto, comp : " 
	<< kcar_auto()(0, 0, 0, 0) * km*km  << "  "
	<< kcar_comp()(0, 0, 0, 0) * km*km << endl ; 
 
    return ost ; 
}
            		    //-------------------------//
			    //	Computational routines //
			    //-------------------------// 

Tenseur Et_bin_ncp::sprod(const Tenseur& t1, const Tenseur& t2) const {
     
  Tenseur* p_tens_metr  ;
  
   // Both indices should be contravariant or both covariant : 
    if (t1.get_type_indice(t1.get_valence()-1) == CON) {
      assert( t2.get_type_indice(0) == CON ) ;
      p_tens_metr = new Tenseur(met_gamma.cov()) ;
	}
    
    if (t1.get_type_indice(t1.get_valence()-1) == COV) {
      assert( t2.get_type_indice(0) == COV ) ;
      p_tens_metr = new Tenseur(met_gamma.con()) ;
       }


  assert ((t1.get_etat() != ETATNONDEF) && (t2.get_etat() != ETATNONDEF)) ;
    // Verifs :
    assert (t1.get_mp() == t2.get_mp()) ;
    
    // Contraction possible ?
     if ( (t1.get_valence() != 0) && (t2.get_valence() != 0) ) {
	    assert ( *(t1.get_triad()) == *(t2.get_triad()) ) ;
    }
    
    int val_res = t1.get_valence() + t2.get_valence() - 2;
    double poids_res = t1.get_poids() + t2.get_poids() ;
    poids_res = (fabs(poids_res) < 1.e-10 ? 0. : poids_res) ;
    const Metrique* met_res = 0x0 ;
    if (poids_res != 0.) {
      assert((t1.get_metric() != 0x0) || (t2.get_metric() != 0x0)) ;
      if (t1.get_metric() != 0x0) met_res = t1.get_metric() ;
      else met_res = t2.get_metric() ;
    }
    Itbl tipe(val_res) ;
    tipe.set_etat_qcq() ;

    for (int i=0 ; i<t1.get_valence() - 1 ; i++)
	tipe.set(i) = t1.get_type_indice(i) ;
    for (int i = t1.get_valence()-1 ; i<val_res ; i++)
	tipe.set(i) = t2.get_type_indice(i-t1.get_valence()+2) ;
	
    Tenseur res(*t1.get_mp(), val_res, tipe, t1.get_triad(), met_res, poids_res) ;

    // Cas particulier ou l'un des deux tenseurs est nul
    if ( (t1.get_etat() == ETATZERO) || (t2.get_etat() == ETATZERO) ) {
	res.set_etat_zero() ; 
	return res ; 
    }

    res.set_etat_qcq() ;
	
    Cmp work(t1.get_mp()) ;
    
    // Boucle sur les composantes de res :
	
    Itbl jeux_indice_t1(t1.get_valence()) ;
    Itbl jeux_indice_t2(t2.get_valence()) ;
    jeux_indice_t1.set_etat_qcq() ;
    jeux_indice_t2.set_etat_qcq() ;
    
    for (int ir=0 ; ir<res.get_n_comp() ; ir++) {    // Boucle sur les composantes
					       // du resultat 

	// Indices du resultat correspondant a la position ir : 
	Itbl jeux_indice_res = res.donne_indices(ir) ;

	// Premiers indices correspondant dans t1 : 
	for (int i=0 ; i<t1.get_valence() - 1 ; i++) {
	    jeux_indice_t1.set(i) = jeux_indice_res(i) ;
	}
	
	// Derniers indices correspondant dans t2 : 
	for (int i=1 ; i<t2.get_valence() ; i++) {
	    jeux_indice_t2.set(i) = jeux_indice_res(t1.get_valence()+i-2) ;
	}
	
	work.set_etat_zero() ;

	// Sommation sur le dernier indice de t1 et le premier de t2 : 
	
	for (int i=0 ; i<3 ; i++) {
	  for (int j=0 ; j<3 ; j++) {
	    jeux_indice_t1.set(t1.get_valence() - 1) = i ;
	    jeux_indice_t2.set(0) = j ;
	    work = work + (*p_tens_metr)(i,j)*t1(jeux_indice_t1)*t2(jeux_indice_t2) ;
	    }
	}
	res.set(jeux_indice_res) = work ;
    
    }	// fin de la boucle sur les composantes du resultat
    delete p_tens_metr ; 
    return res ;

}

