//standard
#include <cstdlib>
#include <cmath>


// Headers Lorene :
#include "nbr_spx.h"
#include "tensor.h" 
#include "cmp.h" 
#include "utilitaires.h"

using namespace Lorene ;

int main() {

    int nz = 2 ;
    double R = 2. ;

    int nr0 = 9 ; int nt0 = 5; int np0 = 4 ;
    
    // echantillonnage en phi :
    int* np = new int [nz] ;
    for (int l=0 ; l<nz ; l++)
	np[l] = np0 ;
    int type_p = SYM ;
    
    // echantillonnage en theta :
    int* nt = new int [nz] ;
    for (int l=0 ; l<nz ; l++)
	nt[l] = nt0 ;
    int type_t = SYM ;

    int* type_r = new int[nz] ;
    type_r[0] = RARE ;
    for (int l=1 ; l<nz ; l++)
	type_r[l] = FIN ;

    type_r[nz-1] = UNSURR ;
    
   // Construction de la grille 

    // echantillonage en r :
    int* nr = new int [nz] ;
    for (int l=0 ; l<nz ; l++)
      nr[l] = nr0 ;
    nr[0] = nr0/2 + 1 ;
	
    Mg3d grille (nz, nr, type_r, nt, type_t, np, type_p) ;

    //Construction du mapping :
    double* bornes = new double[nz+1] ;
    for (int i=0 ; i<nz ; i++)
	bornes[i] = R*double(i)/double(nz - 1) ;
    bornes[nz] = __infinity ;
    Map_af mapping(grille, bornes) ;

    Scalar essai(mapping) ;

    essai = 1 ;

    Scalar deux = essai + 5 ;

    essai = deux ;

    FILE* fich = fopen("res.d", "w") ; 
    
    grille.sauve(fich) ; 
    mapping.sauve(fich) ; 
    essai.sauve(fich) ; 
    
    fclose(fich) ; 
    
    fich = fopen("res.d", "r") ; 

    Mg3d mg1(fich) ; 
    Map_af mp1(mg1, fich) ; 

    Scalar scalfich(mapping, grille, fich) ;

    fclose(fich) ; 
    
    Tensor num(mapping, 1, COV, mapping.get_bvect_spher()) ;

    Itbl ind(1) ;

    ind.set_etat_qcq() ;
    ind.set(0) = 1 ;

    num.set(ind) = essai ;

    cout << num ;
    
    Cmp cuu(mapping) ;
    Cmp theo(mapping) ;
    Coord& rr = mapping.r ;
    Coord& zz = mapping.z ;
	
    cuu = exp(-rr*rr) ;
    cuu.std_base_scal() ;
    essai = cuu ;

    theo = -2*zz*exp(-rr*rr) ;
    theo.va.set_base_r(0,R_CHEBI) ;
    theo.va.set_base_r(1,R_CHEBU) ;
    theo.va.set_base_t(T_COS_I) ;
    theo.va.set_base_p(P_COSSIN_P) ;
    theo.inc2_dzpuis() ;
    scalfich = theo ;

    Scalar resu = scalfich - essai.dsdz() ;
    resu.spectral_display("resu", 1.e-10) ;

    arrete() ; 

    Vector vv(mapping, COV, mapping.get_bvect_spher() ) ;
    vv = num ;

    cout << vv ;

    Sym_tensor tt(mapping, CON, mapping.get_bvect_cart() ) ;

    cout << tt ;

    fich = fopen("res2.d", "w") ; 
    
    grille.sauve(fich) ; 
    mapping.sauve(fich) ; 
    num.sauve(fich) ; 
    
    fclose(fich) ; 
    
    fich = fopen("res2.d", "r") ; 

    Mg3d mg2(fich) ; 
    Map_af mp2(mg2, fich) ; 

    Vector toto(mapping, mapping.get_bvect_spher(), fich) ;

    fclose(fich) ;

    
    delete [] bornes ;
    delete [] nr ;
    delete [] nt ;
    delete [] np ;
    delete [] type_r ;
    
    return EXIT_SUCCESS ; 
}    






