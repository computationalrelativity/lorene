# Multi-grid parameters for the NS
##########################################
8	nz: total number of domains
1
21	nt: number of points in theta (the same in each domain)
20	np: number of points in phi   (the same in each domain)
# Number of points in r and (initial) inner boundary of each domain:
33	0.	<-   nr	  &   min(r)  in domain 0  (nucleus)
33	1.5.	<-   nr	  &   min(r)  in domain 1
33	3 	<-   nr   &   min(r)  in domain 2
33      6.
33      12
33      24 
33      48
33      96  




