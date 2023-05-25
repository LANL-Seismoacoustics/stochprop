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

        Usage: stochprop prop season-trends [OPTIONS]

            stochprop prop season-trends
            -----------------------
        
            Example Usage:
                stochprop prop season-trends --atmo-dir profs/ --results-path example

            Options:
              --atmo-dir TEXT        Directory of atmospheric specifications (required)
              --results-path TEXT    Output path and prefix
              --atmo-pattern TEXT    Specification file pattern (default: '*.dat')
              --atmo-format TEXT     Specification format (default: 'zTuvdp')
              --year-selection TEXT  Limit analysis to specific year(s) (default: None)
              --include-NS BOOLEAN   Option to include north/south analysis
              -h, --help             Show this message and exit.


Stuff...

    .. code:: none

        stochprop prop season-trends --atmo-dir profs/

Stuff...

    .. figure:: _static/_images/season-trends.png
        :width: 600px
        :align: center
        :figclass: align-center


Stuff...

    .. code-block:: none

        ##########################$#########
        ##                                ##
        ##           stochprop            ##
        ##      Propagation Methods       ##
        ##   ESS Ratio Seasonal Analysis  ##
        ##                                ##
        ####################################


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

        Usage: stochprop eof build [OPTIONS]

        stochprop eof build
        -----------------------
        
        Example Usage:
            stochprop eof build --atmo-dir profs/ --eofs-path eofs/example
            stochprop eof build --atmo-dir profs/ --eofs-path eofs/example_low_alt --max-alt 80.0 --eof-cnt 50
            stochprop eof build --atmo-dir profs/ --eofs-path eofs/example_winter --month-selection '10:12, 01:03'

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

        stochprop eof build --atmo-dir profs/ --eofs-path eofs/example_winter --week-selection '38:52,1:15'

Stuff...

    .. code:: none


        ##################################
        ##                              ##
        ##           stochprop          ##
        ##          EOF Methods         ##
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

        stochprop eof build --atmo-dir profs/ --eofs-path eofs/example_summer --week-selection '19:35'

        stochprop eof build --atmo-dir profs/ --eofs-path eofs/example_transition --week-selection '16:18,36:37'  --eof-cnt 35

Note that in the spring/fall analysis, there aren't enough atmospheric specifications in the 5 weeks defining the spring and fall transitions and the methods error out if more EOFs are requested than atmospheric specifications provided.  In more general analysis, sampling these weeks across multiple years provide sufficient atmospheric specification samples to produce a full 100 EOFs, but in this example the EOF count needs to be limited to 35.


**Visualize EOFs**

Discussion...

    .. code:: none

        stochprop eof plot --eofs-path eofs/example_winter


Stuff...

    .. figure:: _static/_images/winter_eofs.png
        :width: 600px
        :align: center
        :figclass: align-center


Repeat for the summer season and show the first 10 EOFs (click the image to view it in a larger format).

    .. figure:: _static/_images/summer_eofs.png
        :width: 800px
        :align: center
        :figclass: align-center



**Sampling the Atmospheric Structure**


Stuff...

    .. code:: none 

        stochprop eof coeffs --atmo-dir profs/ --eofs-path eofs/example_winter --coeff-path coeffs/example_winter --week-selection '38:52,1:15'



**Check Seasonal Trends with EOFs (optional) **

Stuff...

    .. code:: none

        stochprop eof build --atmo-dir profs/ --eofs-path eofs/example_all  --max-alt 80.0

Stuff...


    .. code:: none

        stochprop eof coeffs --atmo-dir profs/ --eofs-path eofs/example_all --run-all-weeks True --coeff-path coeffs/example_all --eof-cnt 50

Stuff...

    .. code:: none

        stochprop eof seasonality --eofs-path eofs/example_all --coeff-path coeffs/example_all --eof-cnt 50



-----------------------------------
Constructing Propagation Statistics
-----------------------------------

Stuff...


--------------------------------
Perturbing Atmospheric Structure
--------------------------------

Stuff...

