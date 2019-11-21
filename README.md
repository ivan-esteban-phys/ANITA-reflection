# ANITA-reflection
Code to compute reflection of polarized radio pulses off an ice surface

If we have an incident flux of (partially) linearly polarized radio pulses phi(sin(epsilon), psi<sub>i</sub>, P<sub>i</sub>), where
* epsilon is the elevation from which the flux comes. epsilon<0 means that the pulse was reflected before observing it
* psi<sub>i</sub> is the initial polarization angle with respect to the horizontal
* P<sub>i</sub> is the initial degree of polarization

this repository containts a simple C++ code that computes the observed flux phi(sin(epsilon), psi<sub>f</sub>, P<sub>f</sub>). Here, psi<sub>f</sub> and P<sub>f</sub> are the observed polarization angle and degree of polarization taking into account reflection for epsilon<0. The program requires no external library. It takes as a command-line argument the incident degree polarization of radio waves (P<sub>i</sub>), and creates a file named "Output_P<sub>i</sub>.dat" with the following columns:

sin(epsilon) psi<sub>f</sub P<sub>f</sub> phi(sin(epsilon), psi<sub>f</sub>, P<sub>f</sub>)

An example compilation command would be
<code> g++ -std=c++11 reflect.cpp -o reflect.out </code>
and to execute the code with, e.g., 90% initial degree of polarization
<code> ./reflect.out 0.9 </code>

The following parameters can be tuned in lines 11-18:
* Height of the observer. The default is 38 km (corresponding to the ANITA balloon).
* Minimum value of sin(epsilon) to scan. The default is sin(-40 degrees).
* Maximum value of sin(epsilon) to scan. The default is 0.
* Number of bins in sin(epsilon). The default is 41.
* Number of bins in psi<sub>f</sub>. The default is 21.
* Number of bins in P<sub>f</sub>. The default is 11.
* Number of random values of psi<sub>i</sub> to generate in order to sample the final flux. The default is 1e7 in order to get small sampling fluctuations
* The initial distribution in polarization angle. The default is a flat distribution, but other distributions such as a Gaussian can be easily implemented with the <random> header in the C++11 standard library

The current hard-coded assumptions are:
* The incident flux is isotropic, i.e., flat in sin(epsilon). This can be easily overcome by weighting the final flux in the output table by the probability of an event coming from the corresponding value of sin(epsilon).
* All waves in the incident flux have the same degree of polarization P<sub>i</sub>.
