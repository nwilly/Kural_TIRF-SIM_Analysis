%{

Create an image using fluorophores at given coordinates placed in a TIRF
field.  It models the psf of fluophores as a 3d gaussian.  The TIRF field
is modeled as a exponential decay.

inputs:
    coord: n x 3 matrix giving x,y,z coordinates of fluorophores.  z=0 
        is the image plane.  Units are nm.

outputs:
    img: simulated fluorescence image using the coordinates and the given
        parameters

Author: Nathan Willy (willy.2@osu.edu)
%}

function img = fluor_sim(coord)

psf_r_xy = 84;% lateral SD of psf 3d gaussian (nm)
psf_r_z = 500;% lateral SD of psf gaussian (nm)
psf_amp = 100;% bightness of the fluorophores (AU)
pixel = 30;% pixel size nm
window = 10;% img will be window x window
tirf_z = 50;% decay constant of tirf field

img = zeros(window);
img_shot = zeros(window);
img_bg = zeros(window);
img_noise = zeros(window);

for(x=1:window)
    for(y=1:window)
        x_nm = x*pixel;
        y_nm = y*pixel;
        for(i=1:size(coord,1))
            img(y,x) = img(y,x) + ...
                1*exp(-(coord(i,3)/tirf_z)^2)*psf_amp*...%excitation
                exp(-(((coord(i,3))/(psf_r_z/3))^2)/2)*...%z psf
                exp(-(((coord(i,1)-x_nm)/(psf_r_xy/3))^2)/2)*...%x psf
                exp(-(((coord(i,2)-y_nm)/(psf_r_xy/3))^2)/2);%y psf
        end
    end
end

            
        