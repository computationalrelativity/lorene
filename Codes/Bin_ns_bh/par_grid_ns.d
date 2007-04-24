# Multi-grid parameters for the NS
##########################################
8	nz: total number of domains
1       nzet number of domains describing the NS
13	nt: number of points in theta (the same in each domain)
12	np: number of points in phi   (the same in each domain)
# Number of points in r and (initial) inner boundary of each domain:
17	0.	<-   nr	  &   min(r)  for all the domain 0  
17      1.5
17	3. 
17      6.
17      12
17      24   
17      48
17      96

