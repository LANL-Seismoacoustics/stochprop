.. _propagation:

=====================================
Propagation Statistics
=====================================

* Overview of building propagation statistics and their use in BISL and SpYE

.. figure:: _static/_images/stochprop_fig1.jpg
    :width: 500px
    :align: center
    :alt: alternate text
    :figclass: align-center
    
    Stochastic propagation models are constructing using a suite of possible atmospheric states, propagation modeling applied to each, and a statistical model describing the variability in the resulting set of predicted effects


********************************
Path Geometry Models (PGMs)
********************************

* Stuff...


.. code:: Python

    from stochprop import propagation

	propagation.run_infraga(sample_dirs + "/" + season, results_dir + "/" + season + "/" + season + ".arrivals.dat", cpu_cnt=cpu_cnt, geom="sph", inclinations=[5.0, 45.0, 1.5], azimuths=azimuths, src_loc=src_loc)

* More stuff...

.. code:: Python

        pgm = propagation.PathGeometryModel()
        pgm.build(results_dir + "/" + season + "/" + season + ".arrivals.dat", results_dir + "/" + season + "/" + season + ".pgm", show_fits=False, geom="sph", src_loc=src_loc)


* More stuff...

.. code:: Python

        pgm.load(results_dir + "/" + season + "/" + season + ".pgm", smooth=True)
        pgm.display(file_id=(results_dir + "/" + season + "/" + season), subtitle=season)

* Use in InfraPy BISL

.. code:: Python

        pgm.display(file_id=(results_dir + "/" + season + "/" + season), subtitle=season)



********************************
Transmission Loss Models (TLMs)
********************************
* Stuff...
