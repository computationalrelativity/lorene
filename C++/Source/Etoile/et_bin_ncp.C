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
 * Revision 1.11  2004/03/25 10:29:04  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.10  2003/12/05 14:50:26  j_novak
 * To suppress some warnings...
 *
 * Revision 1.9  2003/10/13 10:29:36  f_limousin
 * *** empty log message ***
 *
 * Revision 1.8  2003/06/20 14:26:51  f_limousin
 * Add many derivatives of the lapse, shift and 3-metric. Add a new argument conf_flat for the contructors and a new function equilibrium_spher().
 *
 * Revision 1.7  2003/03/03 19:20:29  f_limousin
 * Add new members : tensor stress and Cmp ssjm1_gtildeij. And modify the triad on which the tensors are calculated. Now this triad is the mapping triad.
 *
 * Revision 1.6  2003/02/12 18:44:26  f_limousin
 * Add members metgamma_auto and metgamma_comp.
 *
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
		       bool irrot, bool conf_flat0, const Base_vect& map_triad, const Metrique& flat0, const Tenseur_sym& source) 
             : Etoile_bin(mp_i, nzet_i, relat, eos_i, irrot, map_triad),
	       flat(flat0),
               gtilde(pow(flat0.determinant(),(-1./3.))*source, flat0),
	       gtilde_auto_con(mp_i, 2, CON, map_triad),
	       gtilde_auto(pow(flat0.determinant(),(-1./3.))*source
			   , flat0),
	       gtilde_comp(gtilde_auto),
	       dcov_logn_auto(mp_i, 1, COV, map_triad),
	       dcov_logn(mp_i, 1, COV, map_triad),
	       dcon_logn(mp_i, 1, CON, map_triad),
	       dcovdcov_logn_auto(mp_i, 2, COV, map_triad),
	       gamma(flat0.determinant()),
	       a_car_auto(flat0.determinant()),
	       a_car_comp(flat0.determinant()),
	       dcov_acar(mp_i, 1, COV, map_triad),
	       dcov_acar_auto(mp_i, 1, COV, map_triad),
	       dcovdcov_acar_auto(mp_i, 2, COV, map_triad),
	       hij(0*flat0.con()),
	       hij_auto(0*flat0.con()),
	       hij_comp(0*flat0.con()),
	       ricci_auto(mp_i, 2, COV, map_triad),
	       ricci_scal(mp_i),
	       qq(mp_i),
	       deltakij(shift * source),
	       deltakij_auto(shift * source),
	       kcar_auto(mp_i),
	       kcar_comp(mp_i),
	       kcar_con(mp_i, 2, CON, map_triad),
	       stress(mp_i, 2, CON, map_triad),
      	       ssjm1_a_car(mp_i),
	       ssjm1_hij00(mp_i),
	       ssjm1_hij10(mp_i),
	       ssjm1_hij20(mp_i),
	       ssjm1_hij11(mp_i),
	       ssjm1_hij21(mp_i),
	       ssjm1_hij22(mp_i),
	       decouple(mp_i),
	       conf_flat(conf_flat0)
{
 
  gtilde_auto_con.set_etat_qcq() ;


  for (int i=0; i<3; i++) 
     for (int j=i; j<3; j++) {
       if (i == j) {
	   gtilde_auto_con.set(i,i) = decouple   ; 
       }
       else {
	   gtilde_auto_con.set(i,j) = 0. ;
       }
     }
  
  // All quantities are initialized to zero : 
 

  dcov_logn_auto = 0 ;
  dcov_logn = 0 ;
  dcon_logn = 0 ;
  dcovdcov_logn_auto = 0 ;
  dcov_acar = 0 ;
  dcov_acar_auto = 0 ;
  dcovdcov_acar_auto = 0 ;
  ricci_auto = dcovdcov_logn_auto ;
  ricci_scal = 0 ;
  qq = 0 ;
  deltakij = 0 ;
  deltakij_auto = 0 ;
  kcar_auto = 0 ;
  kcar_comp = 0 ;
  kcar_con = 0 ;
  stress = 0 ;
  ssjm1_a_car = 0 ;
  ssjm1_hij00 = 0 ;
  ssjm1_hij10 = 0 ;
  ssjm1_hij20 = 0 ;
  ssjm1_hij11 = 0 ;
  ssjm1_hij21 = 0 ;
  ssjm1_hij22 = 0 ;
}

// Copy constructor
// ----------------

Et_bin_ncp::Et_bin_ncp(const Et_bin_ncp& et)
      	   : Etoile_bin(et),
	     flat(et.flat),
             gtilde(et.gtilde),	
	     gtilde_auto_con(et.gtilde_auto_con),
	     gtilde_auto(et.gtilde_auto),
	     gtilde_comp(et.gtilde_comp),
	     dcov_logn_auto(et.dcov_logn_auto),
	     dcov_logn(et.dcov_logn),
	     dcon_logn(et.dcon_logn),
	     dcovdcov_logn_auto(et.dcovdcov_logn_auto),
	     gamma(et.gamma),
	     a_car_auto(et.a_car_auto),
	     a_car_comp(et.a_car_comp),
	     dcov_acar(et.dcov_acar),
	     dcov_acar_auto(et.dcov_acar_auto),
	     dcovdcov_acar_auto(et.dcovdcov_acar_auto),
	     hij(et.hij),
	     hij_auto(et.hij_auto),
	     hij_comp(et.hij_comp),
	     ricci_auto(et.ricci_auto),
	     ricci_scal(et.ricci_scal),
	     qq(et.qq),
	     deltakij(et.deltakij),
	     deltakij_auto(et.deltakij_auto),
	     kcar_auto(et.kcar_auto),
	     kcar_comp(et.kcar_comp),
	     kcar_con(et.kcar_con),
	     stress(et.stress),
             ssjm1_a_car(et.ssjm1_a_car),
	     ssjm1_hij00(et.ssjm1_hij00),
	     ssjm1_hij10(et.ssjm1_hij10),
	     ssjm1_hij20(et.ssjm1_hij20),
	     ssjm1_hij11(et.ssjm1_hij11),
	     ssjm1_hij21(et.ssjm1_hij21),
	     ssjm1_hij22(et.ssjm1_hij22),
	     decouple(et.decouple),
	     conf_flat(et.conf_flat)
{
  set_der_0x0() ;    
}


// Constructor from a file  
// -----------------------

Et_bin_ncp::Et_bin_ncp(Map& mp_i, const Eos& eos_i, const Base_vect& 
		       map_triad, const Metrique& flat0, FILE* fich)
  : Etoile_bin(mp_i, eos_i, map_triad, fich),
    flat(flat0),
    gtilde(pow(flat0.determinant(),(-1./3.))*flat0.cov(), 
    	   flat0, flat0),
    gtilde_auto_con(mp_i, map_triad, fich),
    //gtilde_auto(gtilde_auto_con, flat0, flat0),
    gtilde_auto(gtilde),
    gtilde_comp(gtilde_auto),
    dcov_logn_auto(mp_i, 1, COV, map_triad),
    dcov_logn(mp_i, 1, COV, map_triad),
    dcon_logn(mp_i, 1, CON, map_triad),
    dcovdcov_logn_auto(mp_i, 2, COV, map_triad),
    gamma(mp_i, map_triad, fich, &flat0),
    a_car_auto(gamma),
    a_car_comp(gamma),
    dcov_acar(mp_i, 1, COV, map_triad),
    dcov_acar_auto(mp_i, 1, COV, map_triad),
    dcovdcov_acar_auto(mp_i, 2, COV, map_triad),
    hij(0*flat0.con()),
    hij_auto(0*flat0.con()),
    hij_comp(0*flat0.con()),
    ricci_auto(mp_i, 2, COV, map_triad),
    ricci_scal(mp_i),
    qq(mp_i),
    deltakij(shift * flat0.cov()),
    deltakij_auto(shift * flat0.cov()),
    kcar_auto(mp_i),
    kcar_comp(mp_i),
    kcar_con(mp_i, 2, CON, map_triad),
    stress(mp_i, 2, CON, map_triad),
    ssjm1_a_car(mp_i, *(mp_i.get_mg()), fich),
    ssjm1_hij00(mp_i, *(mp_i.get_mg()), fich),
    ssjm1_hij10(mp_i, *(mp_i.get_mg()), fich),
    ssjm1_hij20(mp_i, *(mp_i.get_mg()), fich),
    ssjm1_hij11(mp_i, *(mp_i.get_mg()), fich),
    ssjm1_hij21(mp_i, *(mp_i.get_mg()), fich),
    ssjm1_hij22(mp_i, *(mp_i.get_mg()), fich),
    decouple(mp_i)
{


  fread(&conf_flat, sizeof(bool), 1, fich) ;		
  
  gtilde_auto.set_con().set_etat_qcq() ;
  gtilde_comp.set_con().set_etat_qcq() ;
  
  for(int i=0; i<=2; i++) 
      for(int j=i; j<=2; j++) {
	gtilde_auto.set_con(i,j) = gtilde_auto_con(i,j) ;
      }

  cout << gtilde_auto_con(1,1)(0,0,0,0) << endl ;
  
  gtilde_auto.set_std_base() ;
  
   // All quantities are initialized to zero : 

  dcov_logn_auto = 0 ;
  dcov_logn = 0 ;
  dcon_logn = 0 ;
  dcovdcov_logn_auto = 0 ;
  dcov_acar = 0 ;
  dcov_acar_auto = 0 ;
  dcovdcov_acar_auto = 0 ;
  ricci_auto = dcovdcov_logn_auto ;
  ricci_scal = 0 ;
  qq = 0 ;
  deltakij = 0 ;
  deltakij_auto = 0 ;
  kcar_auto = 0 ;
  kcar_comp = 0 ;
  kcar_con = 0 ;
  stress = 0 ;

 }


			    //------------//
			    // Destructor //
			    //------------//

Et_bin_ncp::~Et_bin_ncp(){

    del_deriv() ; 

}



			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Et_bin_ncp
// --------------------------------
void Et_bin_ncp::operator=(const Et_bin_ncp& et) {

  // Assignment of quantities common to the derived classes of Etoile_bin
  Etoile_bin::operator=(et) ;	   


  // Assignement of proper quantities of class Et_bin_ncp
  flat = et.flat ;
  gtilde = et.gtilde ;
  gtilde_auto_con = et.gtilde_auto_con ;
  gtilde_auto = et.gtilde_auto ;
  gtilde_comp = et.gtilde_comp ;
  dcov_logn_auto = et.dcov_logn_auto ;
  dcov_logn = et.dcov_logn ;
  dcon_logn = et.dcon_logn ;
  dcovdcov_logn_auto = et.dcovdcov_logn_auto ;
  gamma = et.gamma ; 
  a_car_auto = et.a_car_auto ; 
  a_car_comp = et.a_car_comp ; 
  dcov_acar = et.dcov_acar ; 
  dcov_acar_auto = et.dcov_acar_auto ;
  dcovdcov_acar_auto = et.dcovdcov_acar_auto ;
  hij = et.hij ;
  hij_auto = et.hij_auto ;
  hij_comp = et.hij_comp ;
  ricci_auto = et.ricci_auto ;
  ricci_scal = et.ricci_scal ;
  qq = et.qq ;
  deltakij = et.deltakij ;
  deltakij_auto = et.deltakij_auto ;
  kcar_auto = et.kcar_auto ;
  kcar_comp = et.kcar_comp ;
  kcar_con = et.kcar_con ;
  stress = et.stress ;
  ssjm1_a_car = et.ssjm1_a_car ;
  ssjm1_hij00 = et.ssjm1_hij00 ;
  ssjm1_hij10 = et.ssjm1_hij10 ;
  ssjm1_hij20 = et.ssjm1_hij20 ;
  ssjm1_hij11 = et.ssjm1_hij11 ;
  ssjm1_hij21 = et.ssjm1_hij21 ;
  ssjm1_hij22 = et.ssjm1_hij22 ;
  decouple = et.decouple ;
  conf_flat = et.conf_flat ;
}



			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Et_bin_ncp::sauve(FILE* fich) const {
    
    Etoile_bin::sauve(fich) ; 

    gtilde_auto_con.sauve(fich) ;
    gamma.sauve(fich) ;
    
    ssjm1_a_car.sauve(fich) ;
    ssjm1_hij00.sauve(fich) ;
    ssjm1_hij10.sauve(fich) ;
    ssjm1_hij20.sauve(fich) ;
    ssjm1_hij11.sauve(fich) ;
    ssjm1_hij21.sauve(fich) ;
    ssjm1_hij22.sauve(fich) ;
    
    fwrite(&conf_flat, sizeof(bool), 1, fich) ;		


}

// Printing
// --------

ostream& Et_bin_ncp::operator>>(ostream& ost) const {
    
  using namespace Unites ;
    Etoile_bin::operator>>(ost) ; 
    
    ost << endl ; 
    ost << "Star in a binary system with non conformally flat metric" 
	<< endl ; 
    ost << "--------------------------------------------------------" 
	<< endl ; 

    /*  
    ost << "Flat : " << flat << endl ; 
    ost << "Gamma tilde : " << gtilde << endl ;
    ost << "Gamma tilde auto : " << gtilde_auto << endl ;
    ost << "Gamma tilde comp : " << gtilde_comp << endl ;
    ost << "Determinant gamma : " << gamma << endl ;
    ost << "log(determinant) : " << loggamma << endl ;

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
	<< kcar_comp()(0, 0, 0, 0) * km*km << endl ; */
 
    return ost ; 
}
            		    //-------------------------//
			    //	Computational routines //
			    //-------------------------// 

Tenseur Et_bin_ncp::sprod(const Tenseur& t1, const Tenseur& t2) const {
     
    Tenseur* p_tens_metr = 0x0 ;
  
   // Both indices should be contravariant or both covariant : 
    if (t1.get_type_indice(t1.get_valence()-1) == CON) {
      assert( t2.get_type_indice(0) == CON ) ;
      
      Tenseur met_gamma(pow(gamma, 1./3.)*gtilde.cov()) ;
      p_tens_metr = new Tenseur(met_gamma) ;
	}
    
    if (t1.get_type_indice(t1.get_valence()-1) == COV) {
      assert( t2.get_type_indice(0) == COV ) ;

      Tenseur met_gamma(pow(gamma, 1./3.)*gtilde.cov()) ;
      p_tens_metr = new Tenseur(met_gamma) ;
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


void Et_bin_ncp::relaxation(const Et_bin_ncp& star_jm1, double relax_ent, 
			    double relax_met, int mer, int fmer_met) {

    double relax_ent_jm1 = 1. - relax_ent ; 
    double relax_met_jm1 = 1. - relax_met ; 

    ent = relax_ent * ent + relax_ent_jm1 * star_jm1.ent ; 

    if ( (mer != 0) && (mer % fmer_met == 0)) {

	logn_auto = relax_met * logn_auto + relax_met_jm1 * star_jm1.logn_auto ;

	logn_auto_regu = relax_met * logn_auto_regu
	  + relax_met_jm1 * star_jm1.logn_auto_regu ;

	logn_auto_div = relax_met * logn_auto_div
	  + relax_met_jm1 * star_jm1.logn_auto_div ;

	a_car_auto = relax_met * a_car_auto 
	                  + relax_met_jm1 * star_jm1.a_car_auto ;
	
	shift_auto = relax_met * shift_auto 
					+ relax_met_jm1 * star_jm1.shift_auto ;
	for(int i=0; i<=2; i++) {
	  for(int j=i; j<=2; j++) {
	    gtilde_auto.set_con(i,j) = relax_met * (gtilde_auto.con())(i,j) 
	      + relax_met_jm1 * ((star_jm1.get_gtilde_auto()).con())(i,j) ;
	    
	    if (!conf_flat){
		hij_auto.set(i,j) = relax_met * hij_auto(i,j) 
		    + relax_met_jm1 * (star_jm1.get_hij_auto())(i,j) ;
	    }
	  }
	}
	
    }
    
    del_deriv() ; 
    
    equation_of_state() ; 

}


