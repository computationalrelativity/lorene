#include "headcpp.h"
#include "math.h"

#include "tbl.h"
#include "matrice.h"
#include "leg.h"
#include "cheby.h"
		
		//*******************
		// First derivative
		//*******************
Tbl ope_der_cheb (const Tbl& so) {
  	// Verification of size :
	assert (so.get_ndim() == 1) ;
	int size = so.get_dim(0) ; // ??//
	
	// The result is set to zero
	Tbl res2 (size) ; //** pas typé????**
	res2.annule_hard() ;
	
	// the computation
	for (int n=0 ; n<size ; n++)
	     for (int p=n+1 ; p<size ; p++)
		if ((p+n)%2 == 1)
			res2.set(n) += 2*p*so(p) ;
	
	// Normalisation of the first coef :
	res2.set(0) /= 2. ;
	
	return res2 ;

}
// *********************************************************************//
		//*******************
		// Second derivative
		//*******************
Tbl ope_der_sec_cheb (const Tbl& so) {
	// Verification of size :
	assert (so.get_ndim() == 1) ;
	int size = so.get_dim(0) ;
	
	// The result is set to zero
	Tbl res (size) ;
	res.annule_hard() ;
	
	// the computation
	for (int n=0 ; n<size ; n++)
	     for (int p=n+2 ; p<size ; p++)
		if ((p+n)%2 == 0)
			res.set(n) += p*(p*p-n*n)*so(p) ;
	
	// Normalisation of the first coef :
	res.set(0) /= 2. ;
	
	return res ;
}

                           //******************
                           //  Transfer Column
                           //*****************

Tbl trans (const Tbl& so, double d, double e, double f) 
{ 
  	// Verification of size :
	assert (so.get_ndim() == 1) ;
	int size = so.get_dim(0) ;
	
	// The result is set to zero
	Tbl res (size) ;
	res.annule_hard() ;
	
	// the computation
	Tbl A (size);
        Tbl B (size);
	A=ope_der_sec_cheb(so);
        B=ope_der_cheb(so);


    for (int n=0; n<size ; n++)
	res.set(n) = d*A(n) + e*B(n) + f*so(n) ;
	    return res;
}
 
                 //******************
                 // Basis vectors
                 //*****************

Tbl base ( int i, int dim){
    assert (dim > i);
    Tbl res (dim);
    	res.annule_hard() ;
	res.set(i) = 1.;
	return res;
}

      
      

 

//***********************************************************************//
              // Computation of Matrix Lij//

Matrice transfer (int size, double a=1. , double b=-4. , double c=4.) {
     

// on applique maintenant le vecteur de transfert trans aux vecteurs de Tchebyshev pour obtenir la matrice de passage//


    Matrice res3 (size, size);
     res3.annule_hard() ; 
    Tbl column (size);
    column.annule_hard();
    for (int i=0; i<size; i++){
	      column= trans ( base(i, size), a, b, c);
              for (int j=0; j<size; j++)

		  res3.set(j,i)= column (j);
    }
	      return res3;
   
 } 


// COefficients de Chebyshev du membre de droite

Tbl droite (int taille) {
    Tbl res0 (taille);
    res0.annule_hard ();
    Tbl pts (taille);
    pts.annule_hard();
	pts= coloc_cheb(taille);
   
    for (int j=0; j<taille; j++)
        res0.set(j) = exp( pts(j))+ (-4*2.718/(1+2.718*2.718)) ; 
    Tbl res (taille);
    res.annule_hard();
    res= coef_cheb (res0);

    return res;
}
   

// introduction des conditions aux bords****
// Remplacement des deux dernières lignes des objets
 
Tbl droite2 (int taille2) {
    Tbl res (taille2);
    Tbl A (taille2);
    A.annule_hard();
    A= droite (taille2);
	res.annule_hard();
	for (int i=0; i<taille2; i++)
	res.set(i)= A(i);
        res.set(taille2-1)=0. ;
	res.set(taille2 - 2) = 0. ;

	    return res ; 
}

Matrice transfer2 (int size)  
  {  Matrice res (size, size);
    res.annule_hard();
    res= transfer (size);

    for (int i=0; i<size; i++) {
    res.set(size -1, i)=cheby (i, 1.);
    res.set(size -2, i)= cheby (i, -1.);}


return res;
}
    
// Décomposition de Tchebyshev de la solution exacte

Tbl exact (int taille) {
    Tbl res0 (taille);
    res0.annule_hard ();
    Tbl pts (taille);
    pts.annule_hard();
	pts= coloc_cheb(taille);
  

   
    for (int j=0; j<taille; j++)
        res0.set(j) = exp( pts(j))- (sinh (1.)/sinh(2.))*exp(2*pts(j))-2.718/(1+2.718*2.718) ;
    Tbl res (taille);
    res.annule_hard();
    res= coef_cheb (res0);
 
    cout << "la solution analytique du système a pour décomposition de Chebyshev" << endl;
    cout << res << endl;
 

    return res;
}
   
    

// Main code :

int main() {
    int size;
    cout << "Nombre de polynômes?" << endl;
    cin >>  size;
  

   cout << "La matrice de transfert (sans les conditions aux bords) est" << endl;
    cout << transfer (size);
 
   
    Matrice test (size, size);
    test.set_etat_qcq();
    test= transfer2 (size);
    
    Tbl sol (size);
    sol.annule_hard();
    sol= droite2 (size);
 
    Tbl res (size);
    res.annule_hard();

    test.set_lu();
    res= test.inverse(sol);
   
    cout <<"solution du systeme:" << endl;
    cout << res << endl;; 
   

    Tbl reste (size);
    reste.annule_hard();
     Tbl pass (size); 
     pass.annule_hard();
     pass= exact(size);
	 for (int i=0; i<size; i++)
	     reste.set(i)= res(i)- pass (i);
    cout << "comparaison avec la valeur analytique:vecteur des résidus" << endl;
       cout << reste << endl;

     return EXIT_SUCCESS ;
}

