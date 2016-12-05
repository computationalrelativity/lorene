

/*
 *   Copyright (c) 2003 CHABBERT Jean-Philippe
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.3  2016/12/05 16:18:26  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.2  2014/10/06 15:12:51  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2003/02/07 17:31:52  jp_chabbert
 * First version with rotstar input data
 *
 *
 *
 *
 * $Header$
 *
 */

// C++ headers


// C headers
#include <cstdio>
#include <cmath>

#include "nrutil.h"

// Lorene headers
#include "etoile.h"
#include "eos.h"
#include "headcpp.h"
#include "main.h"

void init(char nfile[],double *y,int *pnmax, double *pre,
	  char *fmetric,char *fresult)
     /* lit le fichier contenant les conditions initiales */
{
  int i;
  FILE *pfile=fopen(nfile,"r");
  double angle,r;
  if (pfile == NULL) printf("Erreur d'ouverture du fichier %s\n",nfile);
  for (i=1;i<=7;i++) fscanf(pfile,"%lf",&y[i]);
  fscanf(pfile,"%i",pnmax);
  fscanf(pfile,"%lf",pre);
  fscanf(pfile,"%s",fmetric);
  fscanf(pfile,"%s",fresult);
  angle=y[6]*M_PI/180;
  r=y[2];
  y[6]=1/sqrt(tan(angle)*tan(angle)/(r*r)+1); /* rp */ 
  y[8]=sign(angle)*sqrt(1-y[6]*y[6]); /* phip */
    #ifdef DEBUG
  for (i=1;i<=8;i++) printf("%lf ",y[i]);
  #endif
  fclose(pfile); 
}

void calcul_dt( const Etoile_rot& star, double *y )
     /* Calcule tp0 tel que ds2=0 */
{
  double t,r,theta,phi,rp,thetap,phip,a,b,c,delta,s1,s2;
  point p;
  /* variables intermediaires */
  t=y[1];
  r=y[2];
  theta=y[3];
  phi=y[4];
  rp=y[6];
  thetap=y[7];
  phip=y[8];

  /* Lecture coefs metrique */
   p=interpol2(r,theta,star);

  a=-exp(p.gamma+p.rho)+p.omega*p.omega*exp(p.gamma-p.rho)*r*r*sin(theta)*sin(theta);
  b=-2*p.omega*phip*r*r*sin(theta)*sin(theta)*exp(p.gamma-p.rho);
  c=phip*phip*r*r*sin(theta)*sin(theta)*exp(p.gamma-p.rho)+exp(2*p.alpha)*(rp*rp+r*r*thetap*thetap);
  delta=b*b-4*a*c;
  s1=(-b+sqrt(delta))/(2*a);
  s2=(-b-sqrt(delta))/(2*a);
  printf("delta=%f\ns1=%f s2=%f\n",delta,s1,s2);
  y[5]=s2;
}

double ds2( const Etoile_rot& star, double y1[9] , double y2[9] )
     /* Calcule le ds^2  */
{
  point p;
  double dt2,dr2,dtheta2,dphi2,dt,dphi,out,r,theta;
  r=y1[2];
  theta=y1[3];
  p=interpol2(r,theta,star);
  dt=(y2[1]-y1[1]);
  dt2=dt*dt;
  dr2=(y2[2]-y1[2])*(y2[2]-y1[2]);
  dtheta2=(y2[3]-y1[3])*(y2[3]-y1[3]);
  dphi=(y2[4]-y1[4]);
  dphi2=dphi*dphi;
  out=-exp(p.gamma+p.rho)*dt2+exp(p.gamma-p.rho)*r*r*sin(theta)*sin(theta)*(dphi-p.omega*dt)*(dphi-p.omega*dt)+exp(2*p.alpha)*(dr2+r*r*dtheta2);
  return out;
}

int main (int argc, char **argv)
     /* main pour test intégration */
{
  FILE *pfile;
  /*FILE *ds2_file=fopen("ds2.verif","w");*/
  double y[9],yout[9];
  double x=0.,h=1e-3,re,r,dd=0.;
  int i=0,j,k=0,nmax;
  char fmetric[30],fresult[30];

  if ( argc != 2 ) printf("Usage: test fichier_init\n");
 
 /* Conditions initiales */
  init(argv[1],y,&nmax,&re,fmetric,fresult);
  pfile=fopen(fresult,"w");

  /* Chargement metrique */
  /*  met=load_potentials(fmetric); */ 

    FILE* fich = fopen(fmetric, "r") ;   // open binary file in readonly

    Mg3d mg(fich) ;
    Map_et mp(mg, fich) ;

    Eos* peos = Eos::eos_from_file(fich) ;

    Etoile_rot star(mp, *peos, fich) ;

    fclose(fich) ;

    star.update_metric() ;
    star.equation_of_state() ;
    star.hydro_euler() ;
    // re=star.ray_eq()*10.;
    cout << star << endl;

  calcul_dt(star,y);
  printf("y5=%f\n",y[5]);
  /* integration */
  r=re;
  /* for (i=0;i<nmax;i++) */
  while ( r <= 10*re && i < nmax )
    {
      rk4(y,8,x,h,yout,star,geodesics);
      dd=dd+ds2(star,y,yout);
      for(j=1;j<=8;j++) y[j]=yout[j];
      #ifdef DEBUG
      printf("%i %f %f %f %f %f \n",i,x,y[1],y[2],y[3],y[4]);
      #endif
      if ( k == 10 )
	{
	  fprintf(pfile,"%+028.20e %+028.20e %+028.20e %+028.20e %+028.20e \n",
		  x,y[1],y[2],y[3],y[4]);
	  k=0;
	}
      x=x+h;
      r=y[2];
      i++;
      k++;
    }
  dd=dd/i;
  fclose(pfile);
  /*fclose(ds2_file);*/
  printf("%i pas\nr = %f\nds2 moyen = %e\n",i,r,dd);
  return EXIT_SUCCESS;
}
