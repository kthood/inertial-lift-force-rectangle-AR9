README - INERTIAL MIGRATION RECTANGULAR CHANNEL AR=9 CODE
Kaitlyn Hood 5/23/2017
DOI: 10.5281/zenodo.833527

Cite as:
Kaitlyn Hood. (2017, July 21). kthood/inertial-lift-force-rectangle-AR9: First release of inertial lift force rectangle AR9. Zenodo. http://doi.org/10.5281/zenodo.833527

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
getLiftForce_AR9.m

REQUIRES: HoLeal_channel_aspectratio9_Re1_mesh8_04-05-2016.mat

INPUT: 	x 		- x coordinate of center of particle
		    y 		- y coordinate of center of particle
	    	rsph 	- radius of particle
		    U		  - maximum velocity of background flow
		    rho		- fluid density
		    L		  - side length of square channel
		
OUTPUT:	forcex	- lift force in the x-direction including order rsph^5
	    	forcey 	- lift force in the y-direction including order rsph^5
		
NOTES:	Results may be inaccurate for large particle radius (rsph > .2L)
		

Code computes the fourth and fifth rider terms of the (inertial) lift force of a particle in a square channel.  
Side length of square channel is L, and the coordinates are chosen:
	-4.5L < x < 4.5L, 	-.5L < y < .5L
The result is in dimensional coordinates.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
getLiftForceO4_AR9.m

REQUIRES: HoLeal_channel_aspectratio9_Re1_mesh8_04-05-2016.mat

INPUT: 	x 		- x coordinate of center of particle
		    y 		- y coordinate of center of particle
	    	rsph 	- radius of particle
	    	U		  - maximum velocity of background flow
	    	rho		- fluid density
	    	L		  - side length of square channel
		
OUTPUT:	forcex	- lift force in the x-direction including order rsph^4
		    forcey 	- lift force in the y-direction including order rsph^4
		
NOTES:	Results may be inaccurate for moderate to large particle radius
		

Code computes the fourth order term of the (inertial) lift force of a particle 
in a square channel.  
Side length of square channel is L, and the coordinates are chosen:
	-4.5L < x < 4.5L, 	-.5L < y < .5L
The result is in dimensional coordinates.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
getLiftImage_AR9.m

REQUIRES: HoLeal_channel_aspectratio9_Re1_mesh8_04-05-2016.mat

INPUT: 	x 		- x coordinate of center of particle
	    	y 		- y coordinate of center of particle
    		rsph 	- radius of particle
    		L		- side length of square channel

OUTPUT
  v1 - image velocity due to the stokeslet in the x-direction (2x1 vector)
  v2 - image velocity due to the stokeslet in the x-direction (2x1 vector)

NOTES:	Results may be inaccurate for moderate to large particle radius
		

Code computes the image velocities of the (inertial) lift force of a particle 
in a square channel.  
Side length of square channel is L, and the coordinates are chosen:
	-4.5L < x < 4.5L, 	-.5L < y < .5L
The result is in non-dimensional coordinates. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
liftforce_AR9_channel_plot.m

REQUIRES:  HoLeal_channel_aspectratio9_Re1_mesh8_04-05-2016.mat
	      	 getLiftForce_AR9.m

INPUT:

OUTPUT:		plot of the lift force field (.eps)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
animation_focusing_AR9.m

REQUIRES: 	getLiftForce_AR9.m
		        inertial_constants_AR2_n901_04-18-2016.mat

INPUT:

OUTPUT:		animation of particles moving according to the 				lift forces in a square channel (.gif)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
