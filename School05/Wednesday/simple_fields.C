// Lorene headers
#include "metric.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"

int main() {

  // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
  
    int nz = 3 ; 	// Number of domains
    int nr = 33 ; 	// Number of collocation points in r in each domain
    int nt = 9 ; 	// Number of collocation points in theta in each domain
    int np = 8 ; 	// Number of collocation points in phi in each domain
    int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
    int symmetry_phi = SYM ; // no symmetry in phi
    bool compact = true ; // external domain is compactified
  
    // Multi-domain grid construction:
    
    Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;
	
    cout << mgrid << endl ; 

  
    // Setup of an affine mapping : grid --> physical space (Lorene class Map_af)
    // --------------------------------------------------------------------------

    // radial boundaries of each domain:
    double r_limits[] = {0., 1., 2., __infinity} ; // ?? que veut dire le []??-> collection de valeurs...//
    assert( nz == 3 ) ;  // since the above array describes only 3 domains
  
    Map_af map(mgrid, r_limits) ;   // Mapping construction
  	
    cout << map << endl ;  
       
    // Denomination of various coordinates associated with the mapping 
    // ---------------------------------------------------------------

    const Coord& r = map.r ;      
    const Coord& x = map.x ;      
    const Coord& y = map.y ;      
    const Coord& z = map.z;
   
// définition de la métrique plate fij:  (a partir de maintenant tout se fait en 3+1)//

   
    Sym_tensor fij(map, COV, map.get_bvect_spher()) ; 
    fij.set(1,1) = 1 ; 
    fij.set(1,2) = 0 ; 
    fij.set(1,3) = 0 ; 
    fij.set(2,2) = 1 ; 
    fij.set(2,3) = 0 ; 
    fij.set(3,3) = 1 ; 

    fij.std_spectral_base() ;  // Standar polynomial bases will be used 
                                // to perform the spectral expansions

     
    // déclaration des quantités 3+1 de la métrique de Kerr Schild (classe Tensor, Vector... //  
    
    float a = 1. ; // moment angulaire par unité de masse
    float M= 1. ;//masse du trou noir
    Scalar rho(map);
    Scalar H (map);
    Scalar N (map);//lapse
    Sym_tensor gamij (map, COV, map.get_bvect_spher());// 3-métrique
    Sym_tensor Kij (map, COV, map.get_bvect_spher());// 3-métrique
    Vector l (map, COV, map.get_bvect_spher()); // vecteur moment angulaire qui corrige la 3-metrique
    Vector beta (map, COV, map.get_bvect_spher());//shift
   
    // définition des variables incriminées..

  rho= sqrt ( (1/2) *( r*r - a*a)+ sqrt ((1/4)*( r*r - a*a)*( r*r - a*a)+ a*a*z*z));
  rho.set_domain(0)=1.;
    H= M* (rho*rho*rho)/(rho*rho*rho*rho+a*a*z*z);
    N = 1. / (sqrt( 1+2*H));
    l.set(1)= (rho*x + a*y)/( rho*rho + a*a) ;
    l.set(2)=(rho*y - a*x)/( rho*rho + a*a);
    l.set(3)= z/rho; 
    l.annule_domain(0);
    beta= 2*H* l;
    beta.annule_domain(0);
    for (int i=1; i<4; i++) { 
	for (int j=1; j<4; j++)
	    gamij.set (i, j)= fij(i,j)+2*H*l(i)*l(j);
    }
    gamij.std_spectral_base();
    gamij.annule_domain(0);
    Kij= (1/2*N)*(gamij.derive_lie(beta));// courbure extrinsèque

    Kij.std_spectral_base ();
    Kij.annule_domain(0);


 // construction de la métrique (ca semble obligatoire...notamment pour définir les dérivées covariantes.// 
 
	      Metric gamma (gamij);

// Tenseur de Ricci
  Sym_tensor Kii = Kij.up_down(gamma) ;
	      Sym_tensor Rij= gamma.ricci();

// trace de la courbure extrinsèque:
	      Scalar K= Kij.trace(gamma);

	      // Courbure, ou scalaire de Ricci:
	        Scalar R= gamma.ricci_scal();

// COurbure extrinseque contractée avec elle-meme: 
 
		Scalar KijKij= contract (Kij, 0,1, Kii,0,1);

	   
 
//****************************** 
	      // Il ne reste donc plus qu'à vérifier que les équations d'évolution (si on peut dire) sont vérifiées . Les expressions sont évaluées puis comparées à leur valeur théorique avec la commande max(abs))
// Il est important de noter également que les dérivées présentes dans les équations sont les dérivées covariantes associées à la 3-métrique gamma.//
 
	      Scalar Hamiltonian= R+ K*K - KijKij;
	      cout << "Ecart à la contrainte Hamiltonienne" << endl; 
              cout << max(abs(Hamiltonian)) << endl;
 
 

    return EXIT_SUCCESS ; 

}
  
