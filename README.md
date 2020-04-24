# mumag
Very small µmagnetic solver written in plain C


**Features (and Limitations):**
* One tensor kernel performs all physics (Magnetostatic Field , Exchange and Uniaxial Anisotropy)
* Code doesn't support regions with different parameters
* Bidimensional. Just one layer
* Toy simulator to test integration algorithms
* RK4 solver gives accurate result for NIST's Standard Problem Number 4
* Save Magnetization Plot in PNG format
* Routine to write OVF files
* Demag kernel computed in high precission with Newell formula and cropped to `float` after FFT.

## Acknowledgment
* Thanks to Prof. Dr. Luis López Díaz (USAL) for the original Demagnetizing Tensor Routines in Fortran
* Thanks to Dr. Michele Voto for some tips and talking that helped me to speed up the convolution process
