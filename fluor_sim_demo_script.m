%{

Demonstration of the fluor_sim function which models fluorophores (psf = 3d
Gaussians) in a TIRF field.  The parameters of the model are set inside the
function fluor_sim.

The sample data is the hemisphere of a 60 vertex polyhedron with
icosahedral symmetry.  The radius is made to be 80nm.  This models a half
completed clathrin coat.

Author: Nathan Willy (willy.2@osu.edu)

%}

close all

load('ccs_coords.mat')
img = fluor_sim(coord);

imagesc(img)
axis equal