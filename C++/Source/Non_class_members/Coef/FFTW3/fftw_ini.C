#include <fftw3.h>
#include "tbl.h"

namespace Lorene {

namespace {
  const int nmax = 50 ; //Maximal number of FFT sizes 
  int nworked = 0 ;
  Tbl* tab_tab[nmax] ;
  fftw_plan plan_fft[nmax] ;
  int nb_fft[nmax] ;
}

fftw_plan prepare_fft(int n, Tbl*& pg) {
  int index = -1 ;
  for (int i=0; ((i<nworked) && (index<0)); i++) 
    if (nb_fft[i] == n) index = i ; //Has the plan already been estimated?

  if (index <0) { //New plan needed
    index = nworked ;
    if (index >= nmax) {
      cout << "prepare_fft: " << endl ;
      cout << "too many plans!" << endl ;
      abort() ;
    }
    tab_tab[index] = new Tbl(n) ;
    Tbl& tab = (*tab_tab[index]) ;
    tab.set_etat_qcq() ;
    plan_fft[index] = 
      fftw_plan_r2r_1d(n, tab.t, tab.t, FFTW_R2HC, FFTW_ESTIMATE) ;
    nb_fft[index] = n ;
    nworked++ ;
  }
  assert((index>=0)&&(index<nmax)) ;
  pg = tab_tab[index] ;
  return plan_fft[index] ;
}

}
