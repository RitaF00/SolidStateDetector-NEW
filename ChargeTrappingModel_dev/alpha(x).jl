using Plots
using SpecialFunctions
using Distributions
using JSON
using ProgressMeter

gr()

"""
I want to model a distributin from wich i can take the dimension of the rasius of the Li dendrites as a function of depth.
In this way I can have a more realistic model of the RDR, where the radius of
the Li dendrites is not constant but depends on the depth.

I want also to built it as a funztion of the volume that can be coccupied and the number relative of centers.
"""