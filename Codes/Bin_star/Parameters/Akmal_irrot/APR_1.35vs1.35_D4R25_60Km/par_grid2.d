# Multi-grid parameters 
#######################
4	nz: total number of domains
1	nzet: number of domains inside the star
17	nt: number of points in theta (the same in each domain)
16	np: number of points in phi   (the same in each domain)
# Number of points in r and (initial) inner boundary of each domain:
25	0.	<-   nr	  &   min(r)  in domain 0  (nucleus)  	
25	1.	<-   nr	  &   min(r)  in domain 1
25	1.5	<-   nr   &   min(r)  in domain 2
25	2.5	<-   nr   &   min(r)  in domain 3


