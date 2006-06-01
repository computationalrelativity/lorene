# Multi-grid parameters for the NS
##########################################
6	nz: total number of domains
1
21	nt: number of points in theta (the same in each domain)
20	np: number of points in phi   (the same in each domain)
# Number of points in r and (initial) inner boundary of each domain:
33	0.	<-   nr	  &   min(r)  in domain 0  (nucleus)
33	2.	<-   nr	  &   min(r)  in domain 1
33	4. 	<-   nr   &   min(r)  in domain 2
33      8.
33      16
33      32   




