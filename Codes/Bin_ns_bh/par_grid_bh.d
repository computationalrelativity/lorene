# Multi-grid parameters for the black hole
##########################################
9	nz: total number of domains
13	nt: number of points in theta (the same in each domain)
12	np: number of points in phi   (the same in each domain)
# Number of points in r and (initial) inner boundary of each domain:
17	0.	<-   nr	  &   min(r)  for all the domains
17      0.5
17      1
17      2
17	4
17      8
17      16
17      32
17      64

