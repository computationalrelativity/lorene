
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

int max( int i, int j);

int min( int i, int j);

double sign( double x );

void rk4(double y[], int n, double x, double h, double yout[]
         , const Etoile_rot& star
         , void (*derivs)(double, double [], double [] , const Etoile_rot&));


point interpol2( double r0, double theta, const Etoile_rot& star );

void geodesics( double x, double *y, double *dydx , const Etoile_rot& star );
