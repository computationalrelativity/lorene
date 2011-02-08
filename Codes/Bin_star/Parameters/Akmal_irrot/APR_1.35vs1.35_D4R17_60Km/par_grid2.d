# Multi-grid parameters 
#######################
4	nz: total number of domains
1	nzet: number of domains inside the star
13	nt: number of points in theta (the same in each domain)
12	np: number of points in phi   (the same in each domain)
# Number of points in r and (initial) inner boundary of each domain:
17	0.	<-   nr	  &   min(r)  in domain 0  (nucleus)  	
17	1.	<-   nr	  &   min(r)  in domain 1
17	1.5	<-   nr   &   min(r)  in domain 2
17	2.5	<-   nr   &   min(r)  in domain 3


