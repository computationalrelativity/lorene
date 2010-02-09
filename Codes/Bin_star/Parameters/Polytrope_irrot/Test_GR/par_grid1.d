# Multi-grid parameters
#######################
4	nz: total number of domains
1	nzet: number of domains inside the star
7	nt: number of points in theta (the same in each domain)
4	np: number of points in phi   (the same in each domain)
# Number of points in r and (initial) inner boundary of each domain:
9	0.	<-   nr	  &   min(r)  in domain 0  (nucleus)
9	1.	<-   nr	  &   min(r)  in domain 1
9	1.5	<-   nr   &   min(r)  in domain 2
9	2.5	<-   nr   &   min(r)  in domain 2


