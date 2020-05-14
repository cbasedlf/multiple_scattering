# multiple_scattering
 Scattering of an spherical wavefront over several thin scattering layers


The code generates a point source at a fixed {x,y,z} position. The source generates an spherical wave that propagates to the first scattering layer. Then the field id multiplied by a random phase (representing the thin scattering layer). After that, the field propagates to the next thin layer, and so on. After the last thin layer, the field is propagated once more, and the intensity is measured with a CCD.

Propagation is done through Fresnel propagation (near field)

You can choose the number of layers, the separation between them, and the properties and location of the light source.

You can see the memory effect when you move the source in {x,y,z}, or when you slightly change the wavefront wavelength. The way to control how much memory effect you have is to change the number of layers, the distance between them, and the feature size of the random phase that represents them (small feature size provides higher scattering angles, while big feature size resembles more and more the result of what would be a global phase change, which will do nothing)

It would be iteresting to think about anisotropy factor at some point. How to implement it in the most efficient way? Right now I think we have that effect with the parameters that we play with, but it is not direct to interpret how.
