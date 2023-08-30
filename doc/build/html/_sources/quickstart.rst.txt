.. _quickstart:

==========
Quickstart
==========

Stochastic propagation analysis for infrasound propagation statistics using the *stochprop* library can be separated into a sequence of steps: 1) investigate temporal trends in the atmospheric structure at a given location to identify seasonality and separate historical atmospheric specifications into groups for summer, winter, and the spring/fall transition, 2) apply some statistical means to re-sample the atmospheric structure during each season to produce an ensemble of atmospheric samples that is a reasonable size, 3) run the re-sampled atmosphere ensemble through some propagation modeling methodology and parameterize the computed propagation predictions into a model that can be utilized in localization and characterization analyses.  Beyond these stages, the atmospheric statistics constructed for re-sampling can be used to perturb a reference atmospheric specification for an event-of-interest; though, this last application is an ongoing area of R&D.  

For the following example usage of the software, check that the Anaconda environment is active (:code:`conda activate stochprop_env`) and :code:`cd` into the 'stochprop/examples' directory where the sample data is located. In general, a useful project directory structure for applying stochprop to a given project includes several subdirectories in the main working directory (recommended, but not required): 

    * 'profs' directory with source G2S atmospheric specifications
    * 'eofs' directory to write EOF results into
    * 'coeffs' directory for EOF projection results
    * 'samples' directory for results of EOF coefficient samples, atmospheric fitting, and perturbation results
    * 'prop' directory for propagation statistics model output

The 'profs' directory is included with the *stochprop* git repository and includes a set of G2S specifications for the Quickstart tutorial here.  Also, because a set of example transmission loss models are included, the 'prop' directory is included with the repository.  The other directories will need to be created before continuing the Quickstart:

    .. code-block:: console

        mkdir eofs
        mkdir coeffs
        mkdir samples

---------------
Seasonal Trends
---------------

For mid-latitude regions of the earth, the middle atmosphere winds associated with the circumpolar vortex have persistent orientation (eastward during the winter and westward during the summer in the northern hemisphere) and can be used to identify seasonal trends for a given location.  The temperature maximum at the stratopause is typically cooler than the near-ground temperature so that the middle atmosphere waveguide is dependent on these wind structures and therefore highly directional.  Refraction of acoustic waves back to the ground surface occurs when the sum of the sound speed and along-direction wind at some altitude aloft is greater than the sound speed near the ground surface.  Therefore, one can compute the ratio of the sound speed at various altitudes relative to that at the ground surface to identify when an infrasonic waveguide is present.

The sum of the sound speed and along-direction winds is termed the "effective sound speed" (ESS) and *stochprop* includes a visualization method that analyzes the ratio of the ESS throughout the atmosphere with that near the ground surface for eastward and westward propagation:

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


A year of Ground-to-Space (G2S) atmospheric data is provided with *stochprop* in the 'examples/profs/' directory.  The :code:`ess-ratio` method can ingest these atmospheres to show the seasonal trends in the stratospheric waveguide:

    .. code:: none

        stochprop plot ess-ratio --atmo-dir profs/ --results-path ess-ratio-example

Running this method will output some information to screen and also produce a *matplotlib* figure window showing the below ESS ratio analysis.  In the figure, the upper panel shows the maximum ESS ratio between 40 and 60 km altitude to the east (blue) and west (red).  The colormap in the lower panel shows the ESS ratio throughout the ratio for values greater than unity again to the east and west in blue and red, respectively.  The colormap covers from unity up to a value of 1.1 so easily compare waveguide characteristic at distinct locations.

    .. figure:: _static/_images/example.ess-ratio.png
        :width: 800px
        :align: center
        :figclass: align-center


The information output to screen when running :code:`ess-ratio` methods summarizes the formation and dissipation of the eastward and westward waveguides.  By default, the method only analyzes the eastward and westward components of the ESS.  The :code:`--include-NS` option can be used to also compute the northward and southward components; though, in most cases these components don't provide much additional information on seasonal trends.

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

There is some excess variability for this data set that causes multiple instances of the middle atmosphere waveguide forming and dissipating, but the general result from this analysis is that the eastward waveguide forms in late September and lasts until early April (weeks 38 -- 15) and the westward waveguide forms in early May and lasts until the end of August (weeks 19 - 35).  These seasonal trends will be utilized in constructing atmospheric statistics.

---------------------
Atmosphere Statistics
---------------------

**Build EOFs**

Empirical orthogonal functions (EOFs) are a numerical means of identifying basis functions (similar to eigenvectors) for a vector space.  In the case of atmospheric statistics, a vector can be defined describing the sound speed and wind fields and statistics of the atmospheric structure analyzed using EOF as discussed in the overview of :ref:`analysis`.  The *stochprop* statistics methods include various functions for quantifying statistics related to the atmospheric structure.  In general, the first step in such analysis is to construct a set of EOFs for a given set up atmospheric specifications.  This task is completed using :code:`stochprop stats build-eofs`:

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
        
From the above seasonal trends analysis using :code:`stochprop plot ess-ratio`, it was determined that the winter season corresponds to weeks 38 - 52 and 1 - 15 of the calendar year.  This information can be included in analysis via the :code:`--week-selection` parameter:

    .. code:: none

        stochprop stats build-eofs --atmo-dir profs/ --eofs-path eofs/example_winter --week-selection '38:52,1:15'

Note that the notation *38:52,1:15* is used to denote all weeks between 38 and 52 plus those between 1 and 15.  Running this function reads the atmospheric data from the *profs/* directory, check the file names and/or header info for datetime info to determine whether individual files are within the specified range of weeks, and uses a singular value decomposition (SVD) to  construct the EOFs.  For larger data sets, ingesting large numbers of atmospheric files for analysis can be time consuming; though, for this sample data set it's relatively quick.  

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


This analysis produces a number of output files in the eofs/ directory named according to the :code:`--eofs-path` option.  The contents of these files is summarized in the below table.

+-----------------------------------------+-------------------------------------------------------------------------------------------+
| EOF Output File                         | Description                                                                               |
+=========================================+===========================================================================================+
| eofs/example_winter-mean_atmo.dat       | Mean values of the sound speed, winds, and density                                        |
+-----------------------------------------+-------------------------------------------------------------------------------------------+
| eofs/example_winter-singular_values.dat | Singular values corresponding to each EOF index                                           |
+-----------------------------------------+-------------------------------------------------------------------------------------------+
| eofs/example_winter-snd_spd.eofs        | EOFs for the sound speed, :math:`c = \sqrt{ \gamma \frac{p}{\rho}}`                       |
+-----------------------------------------+-------------------------------------------------------------------------------------------+
| eofs/example_winter-merid_winds.eofs    | EOFs for the meridional (north/south) winds                                               |
+-----------------------------------------+-------------------------------------------------------------------------------------------+
| eofs/example_winter-zonal_winds.eofs    | EOFs for the zonal (east/west) winds                                                      |
+-----------------------------------------+-------------------------------------------------------------------------------------------+

 * NOTE: the current implementation of *stochprop* saves the mean atmospheric structure using only the sound speed, winds, and density.  Other output atmospheric data (samples, perturbations, etc.) are saved in the G2S 'zTuvdp' format with columns containing altitude, temperature, zonal wind, meridional wind, density, and pressure.  The format information needed for the *NCPAprop* package is included in the file header; however, if the mean atmospheric file is used in *infraGA/GeoAc* be sure to specify the column format: 'zcuvd'.  This might be changed in a future update.

Similar analysis can be completed for the summer and spring/fall transition periods when the middle atmosphere waveguide dissipates:

    .. code:: none

        stochprop stats build-eofs --atmo-dir profs/ --eofs-path eofs/example_summer --week-selection '19:35'

        stochprop stats build-eofs --atmo-dir profs/ --eofs-path eofs/example_transition --week-selection '16:18,36:37'  --eof-cnt 35

Note that in the spring/fall analysis, there aren't enough atmospheric specifications in the 5 weeks defining the spring and fall transitions and the methods error out if more EOFs are requested than atmospheric specifications provided.  In more general analysis, sampling these weeks across multiple years provide sufficient atmospheric specification samples to produce a full 100 EOFs, but in this example the EOF count needs to be limited to 35.

**Visualize EOFs**

The EOF analysis results from the EOF construction can be visualized using the EOF function in :code:`stochprop plot`.  All that's requires is to specify the EOFs path from build run:

    .. code:: none

        stochprop plot eofs --eofs-path eofs/example_winter


A figure is generating containing the mean sound speed and wind structures (left-most panels) and the first few EOFs are visualized for comparison.  As noted in the discussion of :ref:`analysis`, the EOFs are computed using stacked sound speed and wind information so the :math:`n^\text{th}` EOF defines perturbations to the sound speed and both wind components in combination.  In the figure, the zonal winds are denoted by the blue lines and meridional by red.

    .. figure:: _static/_images/winter_eofs.png
        :width: 600px
        :align: center
        :figclass: align-center


In addition to visualizing the mean atmosphere and EOF structure, the visualization methods perform an analysis of the singular values to identify the number of EOF terms needed for accurately representing the atmospheric structure.  In the below summary, analysis has been performed to identify the number of singular values needed to describe some percentage of the atmospheric variability by considering the ratio of the :math:`n^\text{th}` EOF's singular value to that of the :math:`n=0` EOF.  In order to include 99% of the atmospheric variability, EOFs with singular values satisfying :math:`\frac{\sigma_n}{\sigma_0} \geq 0.01` must be included.  From the analysis below, this requires inclusion of the first 56 EOFs in analysis.

    .. code:: none

        #######################################
        ##                                   ##
        ##             stochprop             ##
        ##       Visualization Methods       ##
        ##        EOF Analysis Results       ##
        ##                                   ##
        #######################################

        Visualizing EOF results...

        Singular Value Analysis
            90% variability rank: 12
            95% variability rank: 21
            99% variability rank: 56
            99.5% variability rank: 72
            99.9% variability rank: 106 

Note that this visualization also reads in the computed singular values and prints variability ranks.  Such information is useful in determining how many EOFs are needed in continued analysis.  In this case, one finds that the ratio of the singular value of the 56th EOF to that of the rank 0 EOF, :math:`\frac{\sigma_{56}}{\sigma_0}` is less than 0.01 and therefore 99% of the variability in the vector space can be captured by using the first 56 EOFs in analysis.  Similar analysis below using the EOFs to fit a reference atmosphere will further demonstrate how this decision of EOF count can be made visually. 

The EOF visualization methods defaults to show the first 5 EOFs but can be adjusted to show additional contributions using the :code:`--eof-cnt` parameter.  Below is an example of the summer EOF structure showing the first 10 EOFs.  Obviously plotting higher numbers make the plot more difficult to read.  A future update may add an option to plot a specific range or sequence of EOFs (e.g., 10:20 or '1,5,9,13') if such an option is determined to be useful.

    .. code:: none

        stochprop plot eofs --eofs-path eofs/example_summer --eof-cnt 10

    .. figure:: _static/_images/summer_eofs.png
        :width: 800px
        :align: center
        :figclass: align-center



**Analyze Fitting Accuracy**

As noted above in the visualization output, the number of EOFs used in analysis determines how much of the variability in the vector suite is captured.  This can be visualized by fitting a reference atmosphere using a specified number of EOFs.  In the above result, 90% of the variability is captured by the first 12 EOFs.  Using these first 12 EOFs as basis functions to fit an atmospheric state can be done using the :code:`eof-fit` function in the :code:`stochprop plot` methods.

    .. code:: none

        stochprop plot eof-fit --atmo-file profs/g2stxt_2011010118_39.1026_-84.5123.dat --eofs-path eofs/example_winter --eof-cnt 12

    .. figure:: _static/_images/eof-fit_12.png
        :width: 500px
        :align: center
        :figclass: align-center


In this result, the reference atmosphere is projected onto the EOF basis vectors as detailed in :ref:`analysis`.  The resulting estimate is relatively accurate for the sound speed and zonal winds (left and center panels of the figure); however, the meridional winds exhibit a notably poor fit.  Consider instead if the fit uses the first 56 EOFs to capture 99% of the variability in the vector space:

    .. code:: none

        stochprop plot eof-fit --atmo-file profs/g2stxt_2011010118_39.1026_-84.5123.dat --eofs-path eofs/example_winter --eof-cnt 56

    .. figure:: _static/_images/eof-fit_56.png
        :width: 500px
        :align: center
        :figclass: align-center

In this result the structure is much more accurately fit and the finer scale variations in the lower atmosphere are fit relatively accurately (though slightly smoothed).  Repeating this analysis with 72 EOFs to include 99.5% of variability improves the accuracy of these finer scale variations.  

It should be noted that the default behavior of :code:`stochprop stats build-eofs` is to only keep the first 100 EOFs (all singular values are kept) and therefore when the above visualization results are obtained it might be necessary to re-run the EOF construction with an increased EOF count to capture higher accuracy (e.g., 106 EOFs for 99.9% variability captured in the EOF vector space).  In general, the 99.5% accuracy level has been used in evaluations of these methods.

In the case that the EOF fit to the atmospheric structure is needed for continued work, it can be written into a file (again, useful to define a subdirectory in 'samples') using:

    .. code:: none

        mkdir samples/fits/

        stochprop plot eof-fit --atmo-file profs/g2stxt_2011010118_39.1026_-84.5123.dat --eofs-path eofs/example_winter --eof-cnt 56 --output-file 'samples/fits/2011-01-01_fit-56.met'


This produces a file with header information summarizing the EOF fit parameters and the atmospheric data in a format identical to other G2S files:

    .. code:: none

        # Data Source: stochprop v0.1.0
        # Calculated: 2023-08-30 10:54:55.038919
        # Method: Fitting
        # Reference Specification = profs/g2stxt_2011010118_39.1026_-84.5123.dat
        # EOF Set = eofs/example_winter (cwd: /path/to/stochprop/examples)
        # EOF Cnt = 56
        # Fields = [ Z(km), T(K), U(m/s), V(m/s), R(g/cm^3), P(mbar)]
        # The following lines are formatted input for ncpaprop
        #% 0, Z0, km, 0.0
        #% 1, Z, km
        #% 2, T, K
        #% 3, U, m/s
        #% 4, V, m/s
        #% 5, RHO, g/cm3
        #% 6, P, mbar
        0.000000000000000000e+00 2.835923387873799015e+02 5.990006341023006442e+00 2.687147261693745293e-02 1.235309932432432888e-03 1.005431122324828380e+03
        1.000000000000000056e-01 2.828483807135883126e+02 6.076833219597064684e+00 -1.821912076473385100e-02 1.238559089438562863e-03 1.005431122324828380e+03
        ...


**Sampling the Atmospheric Structure**


The primary aim of using EOFs to quantify the statistics of the atmospheric structure over some archived time period is for data reduction.  This can be accomplished as noted in :ref:`analysis` via sampling of coefficient values.  Generating a suite of atmospheric samples representative of a given season requires two steps:  1) compute coefficients for the EOF basis vectors projected onto all elements of the suite, 2) use a kernel density estimate (KDE) of the finite set of coefficient values to produce a set of atmospheric samples representative of the original suite.  The first of these steps can be accomplished using the :code:`eof-coeffs` function in :code:`stochprop stats`.  This analysis requires specification of the EOF basis info, output coefficient info path, and any selection of months or weeks needed to only include the appropriate atmospheric data.  For the example winter analysis:

    .. code:: none 

        stochprop stats eof-coeffs --atmo-dir profs/ --eofs-path eofs/example_winter --coeff-path coeffs/example_winter --week-selection '38:52,1:15'

This analysis produces a single *example_winter-coeffs.npy* file that contains the coefficient values for all EOFs used in analysis.  These coefficients can then be used to sample the atmospheric vector space using :code:`stochprop stats sample-eofs`.  However, it's useful to first create a subdirectory for this set of samples ot keep them separate from other seasonal samples,

    .. code:: none

        mkdir samples/winter 

        stochprop stats sample-eofs --eofs-path eofs/example_winter --coeff-path coeffs/example_winter --sample-path samples/winter/example_winter --sample-cnt 50

The number of samples generated is controlled via :code:`--sample-cnt` and it can be slightly difficult to determine a useful number for continued analysis.  Ideally, the mean and standard deviation computed from individual sample sets should be consistent; therefore, one can repeatedly generate suites and compare these statistics.  That is, if the above method is run using :code:`--sample-cnt 10`, a set of 10 atmospheric specification samples are generated.  If it's run a second time, another set of 10 samples will be generated.  It's unlikely that with only 10 samples the mean and standard deviations of these sets will be consistent; however, if the sample count is increased it's likely that they will converge to the mean and standard deviation of the original set of atmospheres.  In ongoing evaluation of the methods here, a sample count slightly larger than the number of EOFs used in construction has produced sample sets with consistent statistics (e.g., using 75 EOFs to capture ~99.5% of variability and generating 100 atmospheric samples).  A more rigorous means of quantifying the needed number of samples is ongoing R&D.

Once the sampled set of atmospheric specifications has been constructed, they can be visualized using the :code:`atmo-ensemble` method in :code:`stochprop plot`.  In this visualization, the darker/thicker lines denote the mean atmospheric state in the ensemble and the various thinner lines are various samples in the set.  By default, the first 25 samples are visualized, but this can be modified using :code:`--plot-cnt` to include more or fewer samples in the visualization.

    .. code:: none 

        stochprop plot atmo-ensemble --atmo-dir samples/winter/ --atmo-pattern '*.met'

    .. figure:: _static/_images/winter_eof-samples.png
        :width: 500px
        :align: center
        :figclass: align-center



Similar to the fitting output, this produces a file with header information summarizing the EOF sampling parameters and the atmospheric data in a format identical to other G2S files:

    .. code:: none

        # Data Source: stochprop v0.1.0
        # Calculated: 2023-08-30 10:02:07.187834
        # Method: Coefficient KDE Sampling
        # Coeff Label = None
        # EOF Set = eofs/example_winter (cwd: /path/to/stochprop/examples)
        # EOF Cnt = 100
        # Sample: 0/50
        # Fields = [ Z(km), T(K), U(m/s), V(m/s), R(g/cm^3), P(mbar)]
        # The following lines are formatted input for ncpaprop
        #% 0, Z0, km, 0.0
        #% 1, Z, km
        #% 2, T, K
        #% 3, U, m/s
        #% 4, V, m/s
        #% 5, RHO, g/cm3
        #% 6, P, mbar
        0.000000000000000000e+00 2.826111359577614053e+02 9.776727538213997315e-01 1.721851603574645839e+00 1.235309932432432671e-03 1.001952425169503158e+03
        1.000000000000000056e-01 2.819774741146347310e+02 1.301368733916757225e+00 2.318155321169644623e+00 1.238085929951651589e-03 1.001952425169503158e+03
        ...

The 'Coeff Label' field can be set in the CLI call using :code:`--label` and is useful to provide some additional information about the sampling (e.g., :code:`--label 'Sampling of winter atmospheric structure in the western US'`).


**Check Seasonal Trends with EOFs (optional)**

As noted in Blom et al. (2023), for locations away from mid-latitude, the effective sound speed ratio is less useful in identifying seasonal trends and an alternative method is needed.  The EOF coefficient structure across a full year of time can be computed and analyzed to identify those times of the year that have similar structure.  Consider first using the :code:`build-eofs` method to build a full year set of EOFs (without any month or week limitations):

    .. code:: none

        stochprop stats build-eofs --atmo-dir profs/ --eofs-path eofs/example_all  --max-alt 80.0

In this case, the seasonal trends are still expected to be focused in the lower and middle atmosphere, so that the :code:`--max-alt 80.0` parameter limits the EOF construction to these regions of the atmosphere and avoids including the atmospheric tides which are dominated by a 24-hour periodicity due to solar heating.  Given this set of EOFs, one can compute coefficient sets for each week of the year using :code:`eof-coeffs` with the parameter flag :code:`--run-all-weeks` (a similar :code:`--run-all-months` is available for a coarser resolution seasonal analysis, but weekly resolution seems to be more robust).

    .. code:: none

        stochprop stats eof-coeffs --atmo-dir profs/ --eofs-path eofs/example_all --run-all-weeks True --coeff-path coeffs/example_all --eof-cnt 50

The above produces 52 *example_all.week_##-coeffs.py* files in the *coeffs/* directory (one for each week).  In a similar method to the above sampling of coefficient values using KDEs, each weekly set of coefficients can be analyzed and overlap of coefficient values between pairs of weeks computed as discussed in :ref:`analysis`:

    .. code:: none

        stochprop stats coeff-overlap --eofs-path eofs/example_all --coeff-path coeffs/example_all --eof-cnt 50

The resulting overlap can be analyzed using hierarchical clustering to identify those groupings of weeks which exhibit similar EOF projections and therefore similar atmospheric structure.  Visualizing the clustering results:

    .. code:: none

        stochprop plot coeff-overlap --overlap coeffs/example_all-overlap.npy


    .. figure:: _static/_images/example_all-seasonality.png
        :width: 300px
        :align: center
        :figclass: align-center

This identifies seasonal trends such that summer extends from week 20 to 33 and winter covers weeks 38 -- 52 and 1 -- 14.  For comparison, the :code:`stochprop plot ess-ratio` analysis above identified similar seasonal trends though with a slightly longer summer (weeks 19 - 35).  While this additional analysis isn't overtly needed for this mid-latitude location, analysis of seasonal trends near the equatorial and polar regions is often elucidated by this EOF coefficient overlap analysis.



----------------------
Propagation Statistics
----------------------

**Constructing Propagation Statistics**

Once a reduced set of atmospheric specifications representative of a given season and location has been constructed, propagation statistics can be computed using methods in :code:`stochprop prop`.  For localization analysis, the propagation path geometry computed via infrasonic ray tracing (propagation time and arrival range for celerity statistics, back azimuth biases due to cross winds, etc.) can be computed using the `InfraGA/GeoAc <https://github.com/LANL-Seismoacoustics/infraGA>`_ ray tracing methods.  The :code:`build-pgm` methods to build a path geometry model (PGM) requires a directory of atmospheric specifications, an output path, a source location (latitude, longitude, altitude) (typically the reference location of the G2S locations).  As with other functions in *stochprop*, usage information can be displayed with the :code:`--help` flag:

    .. code:: none

        Usage: stochprop prop build-pgm [OPTIONS]

          stochprop prop build-pgm 
          ---------------------
        
          Example Usage:
              stochprop prop build-pgm --atmo-dir samples/winter/ --output-path prop/winter/winter --src-loc '[30.0, -120.0, 0.0]'  --cpu-cnt 8

        Options:
          --atmo-dir TEXT      Directory containing atmospheric specifications
          --atmo-pattern TEXT  Atmosphere file pattern (default: '*.met')
          --output-path TEXT    Path and prefix for PGM output
          --src-loc TEXT        Source location (lat, lon, alt)
          --inclinations TEXT   Inclination min, max, and step (default: [2, 50, 2]
          --azimuths TEXT       Azimuth min, max, and step (default: [0, 360, 6]
          --bounces INTEGER     Number of ground bounces for ray paths (default: 25)
          --z-grnd FLOAT        Ground elevation for simulations (default: 0.0 km)
          --rng-max FLOAT       Maximum range for simulations (default: 1000.0 km)
          --freq FLOAT          Frequency for Sutherland-Bass atten. (default: 0.5 Hz)
          --prof-format TEXT    Profile format (default: 'zTuvdp')
          --infraga-path TEXT   Path to infraGA/Geoac binaries
          --clean-up BOOLEAN    Remove individual results after merge (default: True)
          --cpu-cnt TEXT        Number of CPUs for propagation simulations
          --rng-window FLOAT    Range window in PGM (default: 50 km)
          --rng-step FLOAT      Range resolution in PGM (default: 10 km)
          --az-bin-cnt INTEGER  Number of azimuth bins in PGM (default: 16)
          --az-bin-width FLOAT  Azimuth bin width in PGM (default: 30 deg)
          --verbose BOOLEAN     Output analysis stages as they're done.
          -h, --help            Show this message and exit.

The ray tracing simulation parameters (e.g., inclination and azimuth angle limits and resolutions) can be tuned if necessary.  The default model construction focuses on regional distances (within 1,000 km).  Note that the default configuration assumes the infraGA/GeoAc binaries are on your path, but if that is not the case you can specify where those methods are using the :code:`--infraga-path` parameter (as noted in the :ref:`installation` discussion, a combined PyGS conda environment is not yet implemented for which the *infraga* Python methods can be leveraged without changing environments).

Similar to the sampling analysis, it's useful to define subdirectories for individual seasonal propagation models.  Construction of PGMs for the winter samples generated above can be completed using (update the call with the available number of CPUs; regardless, this will run for a while):

    .. code:: none

        mkdir prop/winter 

        stochprop prop build-pgm --atmo-dir samples/winter/ --output-path prop/winter/winter --src-loc '39.1026, -84.5123, 0.0' --verbose True --cpu-cnt 14 --clean-up False

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
          Source location: 39.1026, -84.5123, 0.0
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
            source location (lat, lon, alt): 39.1026, -84.5123, 0
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


The above will cycle through all atmospheric specifications in the directory and compute ray paths for each one using the inclination and azimuth parameters specified.  Once those computations are complete, arrivals are combined into a single file (*prop/winter/winter.arrivals.dat*) and statistics are constructed for the arrival celerity and back azimuth deviation as a function of range and azimuth.  The width of the sliding range window defaults to 50 km and has a step of 10 km.  These values can be modified via :code:`--rng-window` and :code:`--rng-step`, respectively.  Azimuth bins defuault to a count of 16 with :math:`30^\circ` widths, which can be modified using :code:`--az-bin-cnt` and :code:`--az-bin-width`, respectively.  The resulting model will be saved (*prop/winter/winter.pgm*) using Pickle (output format will likely be updated in the future, considering a Numpy npz file or a JSON file).

Once constructed, the method can be visualized using the :code:`prop-model` method in :code:`stochprop plot`.  For a PGM, two plots are generated showing the statistics at each compass direction (note: only the 8 primary compass directions are shown even if additional azimuth bins are computed).  The first plot shows azimuth deviation statistics with the bias (darker line) and scatter (95% confidence bound in shaded region) while the second shows the probability of an infrasonic arrival with a given celerity.

    .. code:: none

        stochprop plot prop-model --model-file prop/winter/winter.pgm

    .. figure:: _static/_images/US_RM-az_dev.png
        :width: 500px
        :align: center
        :figclass: align-center

    .. figure:: _static/_images/US_RM-rng_cel.png
        :width: 500px
        :align: center
        :figclass: align-center


For source characterization (e.g., yield estimation), transmission loss statistics are needed between the receiver and source.  The normal mode methods in the `NCPAprop <https://github.com/chetzer-ncpa/ncpaprop>`_  software suite can be used to compute the transmission loss using an "N x 2D" framework in discrete azimuths (unlike the fully 3D ray tracing methods used for PGM construction).  Due to recent developments of terrain inclusion for the NCPAprop parabolic equation (PE) methods, an option is planned in the near term to enable TLM construction using either modal methods or PE methods.  Again, usage can be displayed using the :code:`--help` flag.

    .. code:: none

        Usage: stochprop prop build-tlm [OPTIONS]

        stochprop prop build-tlm 
        ---------------------
        
          Example Usage:
            stochprop prop build-tlm --atmo-dir samples/winter/ --output-path prop/winter/winter --freq 0.2  --cpu-cnt 8

        Options:
          --atmo-dir TEXT           Directory containing atmospheric specifications
          --atmo-pattern TEXT       Atmosphere file pattern (default: '*.met')
          --output-path TEXT         Path and prefix for TLM output
          --freq FLOAT               Frequency for simulation (default: 0.5 Hz)
          --azimuths TEXT            Azimuth min, max, and step (default: [0, 360, 6]
          --z-grnd FLOAT             Ground elevation for simulations (default: 0.0)
          --rng-max FLOAT            Maximum range for simulations (default: 1000.0)
          --ncpaprop-path TEXT       Path to NCPAprop binaries
          --clean-up BOOLEAN         Remove individual results after merge (default: True)
          --cpu-cnt TEXT             Number of CPUs for propagation simulations
          --az-bin-cnt INTEGER       Number of azimuth bins in TLM (default: 16)
          --az-bin-width FLOAT       Azimuth bin width in TLM (default: 30 deg)
          --rng-lims TEXT            Range limits in TLM (default: [1, 1000])
          --rng-cnt INTEGER          Range intervals in TLM (default: 100)
          --rng-spacing TEXT         Option for range sampling ('linear' or 'log')
          --use-coherent-tl BOOLEAN  Use coherent transmission loss (default: False
          -h, --help                 Show this message and exit.

As with the PGM construction, a suite of atmospheric specifications is needed as well as an output path and label.  Unlike the ray tracing methods, the azimuthal symmetry assumed in the NCPAprop methods doesn't account for the spherical Earth geometry so that the latitude and longitude are not needed (I need to add the source altitude as an option to build statistics for elevated sources).  Unlike the "infinite frequency" simulation methods using ray tracing, TLM construction requires specification of a frequency for simulation and yield estimation analysis typically requires multiple models covering some frequency range across which signal is observed.  During evaluation of the methods here, frequencies of 0.05, 0.1, 0.2, 0.5, and 1.0 have been used in computation; though, higher resolution models might be useful in more robust analyses.  As with linking infraGA/GeoAc, the path to NCPAprop binaries can be specified via :code:`--ncpaprop-path` if that directory is not on your path. For the winter example being used here, a TLM at 0.2 Hz can be computed using:

    .. code:: none

        stochprop prop build-tlm --atmo-dir samples/winter/ --output-path prop/winter/winter --freq 0.2  --cpu-cnt 8

A note about multi-threading: unlike the multi-threading in infraGA/GeoAc that's built into the software, multi-threading the NCPAprop methods is done via Python's *subprocess* library to simultaneously compute propagation effects for multiple atmospheric specifications at the same time.  Because of this, the maximum number of CPUs useful in such analysis is limited by the number of atmospheric specifications in the suite (for ray tracing, acceleration is limited by the number of unique inclination angles at each azimuth since that's how infraGA/GeoAc is parallelized).  

Once computed, the model is saved into a file named using the specified output path/label and the frequency used in the calculation (*prop/winter/winter_0.200Hz.tlm* in this case).  Visualization of the resulting model can be be done again using the :code:`prop-model` method in the plotting tools:

    .. code:: none

        stochprop plot prop-model --model-file prop/winter/winter_0.200Hz.tlm

    .. figure:: _static/_images/winter_0.359_tloss.png
        :width: 400px
        :align: center
        :figclass: align-center

**Utilizing Propagation Statistics**

The PGM and TLM files constructed above and be ingested for use by `InfraPy <https://github.com/LANL-Seismoacoustics/infrapy>`_ for localization and yield estimation.  Additionally, several more generalized statistical analyses can be completed within *stochprop*.  As discussed in Blom et al. (2023), the IMS infrasound noise background model can be combined with an explosive blastwave model and transmission loss statistics to compute detection probability.  Several TLMs are included in the git repository for *stochprop* describing transmission loss statistics for the western US (in the vicinity of the Rocky Mountains).  These can be used to compute the probability of an arrival with a signal-to-noise ratio (SNR) greater than unity and quantify the likelihood that a signal is detected in such a scenario:

    .. code:: none

        stochprop plot detection-stats --tlm-files 'prop/US_RM/US_RM-winter_*Hz.tlm' --yield-vals '1, 10, 100' --array-dim 5

In this function call, the wildcard in the TLM files specification (note those quotes around it) leads in all available TLMs matching that patter and the set of yield values are defined in tons equivalent TNT.  The resulting visualization is shown below and 

    .. figure:: _static/_images/det-stats.png
        :width: 600px
        :align: center
        :figclass: align-center

In the resulting visualization, each row corresponds to one of the explosive yield values and each column denotes results for a single transmission loss model.  *Note: this visualization method has some weird behavior when using a single TLM and/or a single yield value.  Work is ongoing to debug it.  Also, for other combinations of TLMs (e.g., summer vs. winter), the labelling needs to be corrected.*

In addition to computing detection statistics for a single station, a network of stations can be analyzed to quantify the probability that an explosive source of a given yield will be detected by some limiting number of stations in the network.  Because the information required in such an analysis is a bit more complex, a file containing network info is needed.  The expected format is a comma-separated-value (CSV) file containing the latitude and longitude of each station, the number of sensor elements (to account for signal gain from beamforming analysis), and a TLM.  An example network info file is summarized below for the University of Utah infrasound network.

    .. code:: none

        39.4731, -110.740, 4, prop/US_RM/US_RM-winter_0.500Hz.tlm
        40.6530, -112.119, 4, prop/US_RM/US_RM-winter_0.500Hz.tlm
        38.5337, -113.855, 4, prop/US_RM/US_RM-winter_0.500Hz.tlm
        39.7196, -113.390, 4, prop/US_RM/US_RM-winter_0.500Hz.tlm
        41.6071, -111.564, 4, prop/US_RM/US_RM-winter_0.500Hz.tlm
        37.0109, -113.244, 4, prop/US_RM/US_RM-winter_0.500Hz.tlm
        40.0795, -111.831, 4, prop/US_RM/US_RM-winter_0.500Hz.tlm

It should be noted that all transmission loss models should use the same frequency content, but depending on the spatial extent of the network each station might have a unique model computed for its location.  The analysis is run via :code:`stochprop plot network-performance` and requires specifying this network info as well as a source yield (in tons eq. TNT) and the latitude and longitude bounds of the region of interest.  The default behavior of this analysis is to require at least 3 detecting stations; though, this parameter can be varied with :code:`--min-det-cnt`.

    .. code:: none

        stochprop plot network-performance --network-info network_test.dat --src-yld 10 --lat-min 36 --lat-max 42 --lon-min -117.5 --lon-max -107.5

    .. figure:: _static/_images/network-performance.png
        :width: 600px
        :align: center
        :figclass: align-center

.. 
    COMMENTED OUT SECTION
    --------------------------------
    Perturbing Atmospheric Structure
    --------------------------------

    The current focus of ongoing *stochprop* development is investigation of the EOF methods for atmospheric perturbation studies to quantify propagation uncertainty.  In the prototype function below, a reference atmospheric specification is perturbed using a set of EOFs to produce a suite of samples characterized by a specified standard deviation (10 m/s in this case).

        .. code:: none

            stochprop stats perturb --method eof --atmo-file profs/g2stxt_2011010118_39.1026_-84.5123.dat --eofs-path eofs/example_winter --sample-path samples/perturb/test --std-dev 10.0

    The :code:`atmo-ensemble` visualization method can be used to visualize the resulting atmospheric specification set.

        .. code:: none

            stochprop plot atmo-ensemble --atmo-dir samples/perturb/ --atmo-pattern '*.met' --ref-atmo profs/g2stxt_2011010118_39.1026_-84.5123.dat

        .. figure:: _static/_images/perturb1.png
            :width: 400px
            :align: center
            :figclass: align-center

