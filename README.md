# stochprop

Simulations of infrasonic propagation in the atmosphere typically utilize a single atmospheric 
specification describing the acoustic sound speed, ambient winds, and density as a function of 
altitude.  Due to the dynamic and sparsely sampled nature of the atmosphere, there is a notable 
amount of uncertainty in the atmospheric state at a given location and time so that a more robust 
analysis of infrasonic propagation requires inclusion of this uncertainty.  This Python library, 
stochprop, has been implemented using methods developed jointly by infrasound scientists at 
Los Alamos National Laboratory (LANL) and the University of Mississippi's National Center for 
Physical Acoustics (NCPA).  This software library includes methods to quantify variability in the 
atmospheric state, identify typical seasonal variability in the atmospheric state and generate 
suites of representative atmospheric states during a given season, as well as perform uncertainty 
analysis on a specified atmospheric state given some level of uncertainty.  These methods have 
been designed to interface between propagation modeling capabilities such as InfraGA/GeoAc and 
NCPAprop and signal analysis methods in the LANL InfraPy tool.


## Authorship
stochprop methods were developed by infrasound scientists at LANL as well as UMiss NCPA and has 
been implemented by the LANL Seismoacoustics (LANL-SA) Team.  

## Documentation
The complete documentation can be found in the docs/build/html/ directory of the package (start from index.html)

