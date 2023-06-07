.. _quickstart:

==========
Quickstart
==========

Quickstart discussion

---------------------------
Identifying Seasonal Trends
---------------------------

Stuff...

    .. code-block:: console

        Usage: stochprop plot ess-ratio [OPTIONS]

          stochprop plot ess-ratio
          -----------------------
        
        Example Usage:
            stochprop plot ess-ratio --atmo-dir profs/ --results-path example

        Options:
          --atmo-dir TEXT        Directory of atmospheric specifications (required)
          --results-path TEXT    Output path and prefix
          --atmo-pattern TEXT    Specification file pattern (default: '*.dat')
          --atmo-format TEXT     Specification format (default: 'zTuvdp')
          --year-selection TEXT  Limit analysis to specific year(s) (default: None)
          --include-NS BOOLEAN   Option to include north/south analysis
          --show-title BOOLEAN   Option to display title text
          -h, --help             Show this message and exit.


Stuff...

    .. code:: none

        stochprop prop season-trends --atmo-dir profs/ --results-path example

Stuff...

    .. figure:: _static/_images/example.ess-ratio.png
        :width: 800px
        :align: center
        :figclass: align-center


Stuff...

    .. code-block:: none

        #####################################
        ##                                 ##
        ##            stochprop            ##
        ##      Visualization Methods      ##
        ##   ESS Ratio Seasonal Analysis   ##
        ##                                 ##
        #####################################

        Run summary:
        Source directory: profs/
        Specification pattern: *.dat
        Specification format: zTuvdp

            Loading profiles from profs/ with pattern: *.dat
                Extracted ground elevation: 0.165

        Computing effective sound speed ratio for each day-of-year...

        Eastward waveguide changes...
            Waveguide dissipates: April 10  (yday: 101, week: 14)
            Waveguide forms: April 11  (yday: 102, week: 15)
            Waveguide dissipates: April 12  (yday: 103, week: 15)
            Waveguide forms: September 23  (yday: 267, week: 38)

        Westward waveguide changes...
            Waveguide forms: May 02  (yday: 123, week: 18)
            Waveguide dissipates: May 04  (yday: 125, week: 18)
            Waveguide forms: May 11  (yday: 132, week: 19)
            Waveguide dissipates: August 29  (yday: 242, week: 35)

There is some excess variability that causes multiple instances of the stratospheric waveguide forming and dissipating, but the general result from this analysis is that the eastward waveguide forms in September and lasts until early April (weeks 38 -- 15) and the westward waveguide forms in early May and lasts until the end of August (weeks 19 - 35).  These seasonal trends will be utilized in constructing atmospheric statistics.


------------------------------
EOFs for Atmosphere Statistics
------------------------------


**Build EOFs**

Stuff...

    .. code:: none

        Usage: stochprop stats build-eofs [OPTIONS]

          stochprop stats build-eofs
          --------------------------
        
        Example Usage:
            stochprop stats build-eofs --atmo-dir profs/ --eofs-path eofs/example
            stochprop stats build-eofs --atmo-dir profs/ --eofs-path eofs/example_low_alt --max-alt 80.0 --eof-cnt 50
            stochprop stats build-eofs --atmo-dir profs/ --eofs-path eofs/example_winter --month-selection '10:12, 01:03'

        Options:
          --atmo-dir TEXT          Directory of atmospheric specifications (required)
          --eofs-path TEXT         EOF output path and prefix (required)
          --atmo-pattern TEXT      Specification file pattern (default: '*.dat')
          --atmo-format TEXT       Specification format (default: 'zTuvdp')
          --month-selection TEXT   Limit analysis to specific month(s) (default: None)
          --week-selection TEXT    Limit analysis to specific week(s) (default: None)
          --year-selection TEXT    Limit analysis to specific year(s) (default: None)
          --save-datetime BOOLEAN  Save date time info (default: False)
          --max-alt TEXT           Maximum altitude for trimming data (default: None)
          --eof-cnt INTEGER        Number of EOFs to store (default: 100)
          -h, --help               Show this message and exit.
        
Stuff...

    .. code:: none

        stochprop stats build-eofs --atmo-dir profs/ --eofs-path eofs/example_winter --week-selection '38:52,1:15'

Stuff...

    .. code:: none


        ##################################
        ##                              ##
        ##           stochprop          ##
        ##      Statistics Methods      ##
        ##   Build SVD to Define EOFs   ##
        ##                              ##
        ##################################


        Run summary:
          Source directory: profs/
          Specification pattern: *.dat
          Specification format: zTuvdp
          Limited weeks: ['38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15']
          EOF count: 100
          Output path: eofs/example_winter

            Loading profiles from profs/ with pattern: *.dat
                Weeks filter: ['38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15']
            Building EOFs using SVD...

Stuff...

    .. code:: none

        stochprop stats build-eofs --atmo-dir profs/ --eofs-path eofs/example_summer --week-selection '19:35'

        stochprop stats build-eofs --atmo-dir profs/ --eofs-path eofs/example_transition --week-selection '16:18,36:37'  --eof-cnt 35

Note that in the spring/fall analysis, there aren't enough atmospheric specifications in the 5 weeks defining the spring and fall transitions and the methods error out if more EOFs are requested than atmospheric specifications provided.  In more general analysis, sampling these weeks across multiple years provide sufficient atmospheric specification samples to produce a full 100 EOFs, but in this example the EOF count needs to be limited to 35.


**Visualize EOFs**

Discussion...

    .. code:: none

        stochprop plot eofs --eofs-path eofs/example_winter


Stuff...

    .. figure:: _static/_images/winter_eofs.png
        :width: 600px
        :align: center
        :figclass: align-center


Repeat for the summer season and show the first 10 EOFs (click the image to view it in a larger format).

    .. code:: none

        stochprop plot eofs --eofs-path eofs/example_summer --eof-cnt 10

    .. figure:: _static/_images/summer_eofs.png
        :width: 800px
        :align: center
        :figclass: align-center



**Analyze Fitting Accuracy**

Stuff...

    .. code:: none

        stochprop plot eof-fit --atmo-file profs/g2stxt_2011010118_39.1026_-84.5123.dat --eofs-path eofs/example_winter --eof-cnt 10

    .. figure:: _static/_images/eof-fit_10.png
        :width: 500px
        :align: center
        :figclass: align-center


More stuff...

    .. code:: none

        stochprop plot eof-fit --atmo-file profs/g2stxt_2011010118_39.1026_-84.5123.dat --eofs-path eofs/example_winter --eof-cnt 25

    .. figure:: _static/_images/eof-fit_25.png
        :width: 500px
        :align: center
        :figclass: align-center




**Sampling the Atmospheric Structure**


Stuff...

    .. code:: none 

        stochprop stats eof-coeffs --atmo-dir profs/ --eofs-path eofs/example_winter --coeff-path coeffs/example_winter --week-selection '38:52,1:15'

Use the coefficients to sample...

    .. code:: none

        stochprop stats sample-eofs --eofs-path eofs/example_winter --coeff-path coeffs/example_winter --sample-path samples/winter/example_winter --sample-cnt 50


Visualize the samples (need to write function)...

    .. code:: none 

        stochprop plot atmo-ensemble --atmo-dir samples/winter/ --atmo-pattern '*.met'

    .. figure:: _static/_images/winter_eof-samples.png
        :width: 500px
        :align: center
        :figclass: align-center


Stuff...


    .. code:: none

        stochprop plot atmo-ensemble --atmo-dir profs/ --atmo-pattern '*.dat' --week-selection '38:52,1:15'

    .. figure:: _static/_images/winter_g2s-atmos.png
        :width: 500px
        :align: center
        :figclass: align-center



**Check Seasonal Trends with EOFs (optional)**

Stuff...

    .. code:: none

        stochprop stats build-eofs --atmo-dir profs/ --eofs-path eofs/example_all  --max-alt 80.0

Stuff...


    .. code:: none

        stochprop stats eof-coeffs --atmo-dir profs/ --eofs-path eofs/example_all --run-all-weeks True --coeff-path coeffs/example_all --eof-cnt 50

Stuff...

    .. code:: none

        stochprop stats coeff-overlap --eofs-path eofs/example_all --coeff-path coeffs/example_all --eof-cnt 50

Visualize...

    .. code:: none

        stochprop plot coeff-overlap --overlap coeffs/example_all-overlap.npy


    .. figure:: _static/_images/example_all-seasonality.png
        :width: 300px
        :align: center
        :figclass: align-center


This identifies seasonal trends such that summer extends from week 20 -- 33 and winter covering weeks 38 -- 52 and 1 -- 14.  For comparison, the :code:`stochprop plot ess-ratio` analysis above identified similar seasonal trends though with a slightly longer summer (weeks 19 - 35).  While this additional analysis isn't overtly needed for this mid-latitude location, analysis of seasonal trends near the equatorial and polar regions is often elucidated by this EOF coefficient overlap analysis.


----------------------
Propagation Statistics
----------------------

**Constructing Propagation Statistics**

Stuff...

    .. code:: none

        stochprop prop build-pgm --atmos-dir samples/winter/ --output-path prop/winter/winter --src-loc '30.0, -120.0, 0.0' --verbose True --cpu-cnt 14 --clean-up False

Output...

    .. code:: none

        #####################################
        ##                                 ##
        ##            stochprop            ##
        ##       Propagation Methods       ##
        ##       Path Geometry Model       ##
        ##                                 ##
        ######################################


        Data IO summary:
          Atmospheric specifications directory: samples/winter/
          Specification pattern: *.met
          Model output path: prop/winter/winter

        infraGA/GeoAc parameters:
          Source location: 30.0, -120.0, 0.0
          Inclination angles (min, max, step): [2.0, 50.0, 2.0]
          Azimuth angles (min, max, step): [0.0, 360.0, 6.0]
          Bounces: 25
          Ground elevation: 0.0
          Range max: 1000.0
          Frequency: 0.5
          Clean up: True
          CPU count: 8

        Path Geometry Model (PGM) parameters:
          Range window: 50.0
          Range step: 10.0
          Azimuth bin count: 16
          Azimuth bin width: 30.0

        Generating ray paths for example_winter-00.met  

            ###########################################
            ####     Running infraga-accel-sph     ####
            ####            Propagation            ####
            ###########################################

        Interpolating atmosphere data in 'samples/winter//example_winter-00.met' using format 'zTuvdp'...
            Propagation region limits:
                latitude = -90, 90
                longitude = -180, 180
                altitutde = 0, 150

            User selected range maximum = 1000

        Parameter summary:
            inclination: 2, 50, 2
            azimuth: 0, 360, 6
            bounces: 25
            source location (lat, lon, alt): 30, -120, 0
            ground elevation: 0
            frequency: 0.5
            S&B atten coeff: 1
            write_atmo: false
            write_rays: false
            write_topo: false
            calc_amp: false
            threads: 8

        Calculating ray paths: (2, 28) degrees inclination range, 0 degrees azimuth.
        Calculating ray paths: (30, 50) degrees inclination range, 0 degrees azimuth.
        Calculating ray paths: (2, 28) degrees inclination range, 6 degrees azimuth.
        Calculating ray paths: (30, 50) degrees inclination range, 6 degrees azimuth.
        Calculating ray paths: (2, 28) degrees inclination range, 12 degrees azimuth.
        Calculating ray paths: (30, 50) degrees inclination range, 12 degrees azimuth.
        ...

        Generating ray paths for example_winter-01.met
        ...



Visualize...

    .. code:: none

        stochprop plot prop-model --model-file prop/winter/winter.pgm

    .. figure:: _static/_images/US_RM-az_dev.png
        :width: 500px
        :align: center
        :figclass: align-center

The second plot...

    .. figure:: _static/_images/US_RM-rng_cel.png
        :width: 500px
        :align: center
        :figclass: align-center



Build a transmission loss model...

    .. code:: none

        stochprop prop build-tlm --atmos-dir samples/winter/ --output-path prop/winter/winter --freq 0.2  --cpu-cnt 8

Visualize...

    .. code:: none

        stochprop plot prop-model --model-file prop/winter/winter_0.200Hz.tlm


    .. figure:: _static/_images/winter_0.359_tloss.png
        :width: 400px
        :align: center
        :figclass: align-center


*Note: this visualization method has some weird behavior when using a single TLM and/or a single yield value.  Work is ongoing to debug it.*

**Utilizing Propagation Statistics**

Plotting detection statistics for a single infrasound array...

    .. code:: none

        stochprop plot detection-stats --tlm-files 'prop/US_RM/US_RM-winter_*Hz.tlm' --yield-vals '1, 10, 100' --array-dim 5


    .. figure:: _static/_images/det-stats.png
        :width: 600px
        :align: center
        :figclass: align-center


Also able to apply models to a network of infrasound arrays to quantify detection capabilities...network info...


Write a file that contains latitude and longitude of each network station as well as the number of sensors for each one and a transmission loss model for each...

    .. code:: none

        39.4731, -110.740, 4, prop/US_RM/US_RM-winter_0.500Hz.tlm
        40.6530, -112.119, 4, prop/US_RM/US_RM-winter_0.500Hz.tlm
        38.5337, -113.855, 4, prop/US_RM/US_RM-winter_0.500Hz.tlm
        39.7196, -113.390, 4, prop/US_RM/US_RM-winter_0.500Hz.tlm
        41.6071, -111.564, 4, prop/US_RM/US_RM-winter_0.500Hz.tlm
        37.0109, -113.244, 4, prop/US_RM/US_RM-winter_0.500Hz.tlm
        40.0795, -111.831, 4, prop/US_RM/US_RM-winter_0.500Hz.tlm



Run the analysis and visualize...

    .. code:: none

        stochprop plot network-performance --network-info network_test.dat --lat-min 36 --lat-max 42 --lon-min -117.5 --lon-max -107.5

    .. figure:: _static/_images/network-performance.png
        :width: 600px
        :align: center
        :figclass: align-center


--------------------------------
Perturbing Atmospheric Structure
--------------------------------

Stuff...

    .. code:: none

        stochprop stats perturb --method eof --atmo-file profs/g2stxt_2011010118_39.1026_-84.5123.dat --eofs-path eofs/example_winter --sample-path samples/perturb/test

Visualization...

    .. codee:: none

        stochprop plot atmo-ensemble --atmo-dir samples/perturb/ --atmo-pattern "*.met" 

    .. figure:: _static/_images/perturb1.png
        :width: 400px
        :align: center
        :figclass: align-center

