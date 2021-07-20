.. _cookbook:

=============================
`GeodeZYX Toolbox`'s Cookbook
=============================

------------
Convert time
------------

The `GeodeZYX Toolbox` can handle a lot of time **scales** and time **representations**. 

What we call `scale` is a physical definition: UTC, GPS Time, TAI...

What we call `representation` is the way the physical value is represented : Julian date, year/day of year, GPS week/day of week, ISO string...

The basic idea behind it is the conversion to and from the ``datetime`` objects of the Python's ``datetime`` module.

The generic structure of the fuctions is ``conv.<input time scale/representation>2<output time scale/representation>``.

Time conversion exemples
~~~~~~~~~~~~~~~~~~~~~~~~

First of all:
::

    import geodezyx.conv as conv                  # Import the conversion module
    import datetime as dt

Get today as GPS Time (week/day of week)
""""""""""""""""""""""""""""""""""""""""
::

    now = dt.datetime.now()    ###  datetime.datetime(2021, 6, 18, 15, 53, 56, 245345)
    gpstime = conv.dt2gpstime(now)

    gpstime
    Out: (2162, 5)


Convert year/day of year to Modified Julian Day
"""""""""""""""""""""""""""""""""""""""""""""""
::

    year,doy = (2021, 169)
    epoch = conv.doy(year, doy)
    mjd = conv.dt2MJD(epoch)

    mjd
    Out: 59383.0

Convert a SP3 name to decimal year
""""""""""""""""""""""""""""""""""
::

    p = "/home/user/igs20726.sp3"

    epoch = conv.sp3name2dt(p)
    yeardec = conv.dt2year_decimal(epoch)

    yeardec
    Out: 2019.739611872146

It handles also RINEX names with ``conv.rinexname2dt(p_rinex)`` (old and new naming)
(Documentation here: :py:func:`geodezyx.conv.conv_time.rinexname2dt`)

Complete reference
~~~~~~~~~~~~~~~~~~

All the module's functionalites can be found here:
:py:mod:`geodezyx.conv.conv_time`

-------------------
Convert coordinates 
-------------------

The `GeodeZYX Toolbox` can easily handle coordinate conversion in Geocentric (X,Y,Z), Geographic (latitude, longitude, height) and topocentric (East, North, Up). 

**Warning**: This is considered as the "low-level" coordinate conversions. It does not deal with the different Reference Frame and their realisations (ITRFxx, ETRFxx...). This is managed by "high-level" functions in the ``reffram`` module.

These functions are optimized for arrays (multiple inputs) but can also handle scalars.

Coordinates conversion exemples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We consider as exemple the coordinates of the Helmertturm given in latitude, longitude, height from Wikipedia and the GNSS station POTS in Geocentric XYZ.
::

    Helmertturm = np.array([52.380278, 13.065278,147.983])
    POTS = np.array([3800689.6341,882077.3857,5028791.3179 ])

We can bring POTS in Geographic and the Helmertturm in Geocentric coordinates
::

    POTS_geo = conv.XYZ2GEO(POTS[0],POTS[1],POTS[2])
    POTS_geo = conv.XYZ2GEO(*POTS)

    POTS_geo
    Out: (52.37929737808202, 13.066091316954145, 144.41769897658378)

::

    Helmertturm_xyz = conv.GEO2XYZ(Helmertturm[0],
                                   Helmertturm[1],
                                   Helmertturm[2])

    Helmertturm_xyz = conv.GEO2XYZ(*Helmertturm)

    Helmertturm_xyz
    Out: (89.92804903383008, 14.005569048843014, -6356604.297220887)


We can also get the vector between POTS and the Helmertturm 
(i.e. the Helmertturm in the Topocentric frame centered on POTS)

::

    conv.XYZ2ENU_2(Helmertturm[0], Helmertturm[1], Helmertturm[2], 
                   POTS[0],POTS[1],POTS[2])
    Out: (array([0.88515353]), array([20735.55848087]), array([-6364723.46820732]))


Complete reference
~~~~~~~~~~~~~~~~~~

All the module's functionalites can be found here:
:py:mod:`geodezyx.conv.conv_coords`

-----------------------
Helmert Transformations
-----------------------

Apply and estimate parameters for an ad-hoc Helmert Transformation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The functions :py:func:`geodezyx.reffram.geometry.helmert_trans_estim` and :py:func:`geodezyx.reffram.geometry.helmert_trans_apply` and  offer an interface to estimate parameters and apply an Helmert transformation respectively.

Transform coordinates between two ITRF/ETRF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The functions :py:func:`geodezyx.reffram.geometry.itrf_helmert_get_parameters` and :py:func:`geodezyx.reffram.geometry.itrf_helmert_trans` offer an interface to get the transformation parameters between two ITREF/ETRF realization, and apply the corresponding Helmert transformation respectively.

Complete reference
~~~~~~~~~~~~~~~~~~

All the module's functionalites can be found here:
:py:mod:`geodezyx.reffram.geometry`

------------------------
Euler pole determination
------------------------

The toolbox proposes tools to manipulate Euler rotation poles:

- to determine the tectonic plate's Euler pole based on some GNSS absolute velocities (:py:func:`geodezyx.geodyn.euler_pole_calc.euler_pole_calc`).
- to analyze the quality of the Euler Pole estimation (:py:func:`geodezyx.geodyn.euler_pole_calc.euler_pole_quality`).
- to convert the estimated Euler pole in a vector form to a latitude/longitude/rate form (:py:func:`geodezyx.geodyn.euler_pole_calc.euler_pole_vector_to_latlongrate`), and also do the reverse conversion (:py:func:`geodezyx.geodyn.euler_pole_calc.euler_pole_vector_from_latlongrate`).
- to substract the plate's velocity to analyze the residual velocities of the stations (located at the plate's boundary for instance) (:py:func:`geodezyx.geodyn.euler_pole_calc.euler_vels_relative_to_ref`).

Complete reference
~~~~~~~~~~~~~~~~~~

All the module's functionalites can be found here:
:py:mod:`geodezyx.geodyn.euler_pole_calc`


---------------------------------
Read and import geodetic products
---------------------------------

Main import functionalities
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The toolbox mainly handles:

- The GNSS products such as clock offsets (.clk files) and orbits (.sp3 files): :py:mod:`geodezyx.files_rw.read_gnss_prods`
- The Earth orientation parameters: :py:mod:`geodezyx.files_rw.read_gnss_prods`
- The Troposphere files: :py:mod:`geodezyx.files_rw.read_athmo`

Complete reference
~~~~~~~~~~~~~~~~~~

All the module's functionalites can be found here:
:py:mod:`geodezyx.files_rw`

------------------------------------
Read and import geodetic time series
------------------------------------

The toolbox is designed to import and pre-process a wide range of geodetic GNSS Time Series.

Read the dedicated Jupyter notebook stored in  ``<...>/geodezyx/000_exemples/timeseries_reader``

Complete reference
~~~~~~~~~~~~~~~~~~

All the module's functionalites can be found here:
:py:mod:`geodezyx.files_rw.read_coords_time_series`

-----------------------------
Read and import GNSS Sitelogs
-----------------------------

You can easily import the content of the IGS's GNSS Sitelogs (a.k.a Logsheets) in dedicated objects

Read the dedicated Jupyter notebook stored in  ``<...>/geodezyx/000_exemples/logsheets_reader``

------------------------------------------
Point and Click to detect offsets manually
------------------------------------------

The toolbox contains a tool to select manually the jumps in the Geodetic Time Series.

Based on `matplotlib`, you can "point and click" the discontinuities you detected visually with your mouse.

Plot first your data (in a theoretical DataFrame DF)

::

    fig,(axn,axe,axu) = plt.subplots(3,1)
    axn.plot(DF["t"],DF["n"])
    axe.plot(DF["t"],DF["e"])
    axu.plot(DF["t"],DF["u"])


Then, create the Point and Click object
::

    PnC = gcls.point_n_click_plot()
    multi , cid = PnC(fig=fig)



The selected jumps/offsets are stored in a list attribute of the Point and Click object
::

    PnC.selectedX


Read the dedicated script stored in  ``<...>/geodezyx/000_exemples/logsheets_reader``.
More details in the class documentation : :py:class:`geodezyx.time_series.ts_class.point_n_click_plot`

----------------------------------------------------
Statistics and plots for orbit and clock comparisons
----------------------------------------------------

TBC






