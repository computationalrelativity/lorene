# Multi-grid parameters 
#######################
4	nz: total number of domains
1	nzet: number of domains inside the star
25	nt: number of points in theta (the same in each domain)
24	np: number of points in phi   (the same in each domain)
# Number of points in r and (initial) inner boundary of each domain:
33	0.	<-   nr	  &   min(r)  in domain 0  (nucleus)  	
33	1.	<-   nr	  &   min(r)  in domain 1
33	1.5	<-   nr   &   min(r)  in domain 2
33	2.5	<-   nr   &   min(r)  in domain 3


