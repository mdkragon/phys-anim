README.txt

We extended from Boom Blox. So, most of the controlls are the same as the boom blox controlls. The only performance difference is we added a few more parameters into the materials xml parser. These new parameters are not needed for compiling, but they are needed if you wish to use different parameters for the material.

Angry_Spheres is made by Tiju Thomas, who has let us use his code base. We imported our sound code base into his code, and made a few changes so the two would run togther smoother (mostly to dampen the threshold in velocity it would take for a sound to be generated). As a result of this, if an object is traveling at a relatively low velocity and collides, there is a chance that sound might not be generated.

To compile both of these sets of code, nothing special is needed to be done. We have made the code base self contained, and pressing run on either debug or release in vs2010 should get the code running. 

