# Parameters for the computation of Kerr black hole in dirac gauge
#######################
3	nz: total number of domains
7	nt: number of points in theta (the same in each domain)
8	np: number of points in phi   (the same in each domain)
# Number of points in r and (initial) inner boundary of each domain:
9	0.	<-   nr	  &   min(r)  in domain 0  (nucleus)
9	1.	<-   nr	  &   min(r)  in domain 1
9	2.	<-   nr	  &   min(r)  in domain 1
0.99	parametre hh
1.	masse M
1.e-6 	seuil : Threshold on xsi relative change for ending the computation


