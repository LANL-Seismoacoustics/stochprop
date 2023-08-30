
=============================
Overview
=============================

Simulations of infrasonic propagation in the atmosphere typically utilize a single atmospheric specification describing the acoustic sound speed, ambient winds, and density as a function of altitude.  Due to the dynamic and sparsely sampled nature of the atmosphere, there is a notable amount of uncertainty in the atmospheric state at a given location and time so that a more robust analysis of infrasonic propagation requires inclusion of this uncertainty.  This Python library, *stochprop*, has been implemented using methods developed jointly by infrasound scientists at Los Alamos National Laboratory (LANL) and the University of Mississippi's National Center for Physical Acoustics (NCPA).  This software library includes methods to quantify variability in the atmospheric state, identify typical seasonal variability in the atmospheric state and generate suites of representative atmospheric states during a given season, as well as perform uncertainty analysis on a specified atmospheric state given some level of uncertainty.  These methods have been designed to interface between propagation modeling capabilities such as InfraGA/GeoAc or NCPAprop and signal analysis methods in the LANL InfraPy tool.  



**License**

© 2023. Triad National Security, LLC. All rights reserved.

This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are.
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare.
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit.
others to do so.
 
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
