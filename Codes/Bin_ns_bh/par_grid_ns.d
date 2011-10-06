# Multi-grid parameters for the NS
##########################################
8	nz: total number of domains
1       number of domains inside the star
11	nt: number of points in theta (the same in each domain)
8	np: number of points in phi   (the same in each domain)
# Number of points in r and (initial) inner boundary of each domain:
13	0.	<-   nr	  &   min(r)  in domain 0  (nucleus)
13	1.5.	<-   nr	  &   min(r)  in domain 1
13	3 	<-   nr   &   min(r)  in domain 2
13      6.
13      12
13      24 
13      48
13      96  




