########## PARAMETERS FOR CODE VISU_EVOL ########################
# The following line contains the root for the input data file names:
sigma
2       jmin: initial time step
200     jmax: final time step
20      jstep: distance between two successive time steps
# The following line contains the root for the OpenDX file names to be created:
mu
z       section_type: z (resp. x, y) for a plane z=a (resp. x=a, y=a) with a=const,  
1.       parameter a defining the section plane (see above)
-4. 4.    umin, umax: range of the plane coordinate u
-4. 4.    vmin, vmax: range of the plane coordinate v
200     nu : number of points in the u direction (uniform sampling)
200     nv : number of points in the v direction (uniform sampling)
