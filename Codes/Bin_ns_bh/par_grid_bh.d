# Multi-grid parameters for the black hole
##########################################
9	nz: total number of domains
21	nt: number of points in theta (the same in each domain)
20	np: number of points in phi   (the same in each domain)
# Number of points in r and (initial) inner boundary of each domain:
33	0.	<-   nr	  &   min(r)  for all the domains
33      0.5
33      1
33      2
33	4
33      8
33      16
33      32
33      64

