/* Définition d'une structure metrique pour sauver les coefs de la metrique */

struct point
{
  double rho;
  double gamma;
  double alpha;
  double omega;
  double drhodr;
  double dgammadr;
  double dalphadr;
  double domegadr;
  double drhodtheta;
  double dgammadtheta;
  double dalphadtheta;
  double domegadtheta;
};
typedef struct point point;

int max( int i, int j)
{
  if ( i<=j ) return j;
  else return i;
}

int min( int i, int j)
{
  if ( i>=j ) return j;
  else return i;
}


double sign( double x )
{
  if ( x >= 0 ) return 1. ;
  else return -1.;
}
