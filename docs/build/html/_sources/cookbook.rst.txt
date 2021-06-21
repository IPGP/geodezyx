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

Complete reference
==================

All the module's functionalites can be found here:
:py:mod:`geodezyx.conv.conv_time`


Time conversion exemples
========================

First of all:
::

    import geodezyx.conv as conv                  # Import the conversion module
    import datetime as dt

Get today as GPS Time (week/day of week)
----------------------------------------
::

    now = dt.datetime.now()    ###  datetime.datetime(2021, 6, 18, 15, 53, 56, 245345)
    gpstime = conv.dt2gpstime(now)

    gpstime
    Out: (2162, 5)


Convert year/day of year to Modified Julian Day
-----------------------------------------------
::

    year,doy = (2021, 169)
    epoch = conv.doy(year, doy)
    mjd = conv.dt2MJD(epoch)

    mjd
    Out: 59383.0

Convert a SP3 name to decimal year
----------------------------------
::

    p = "/home/user/igs20726.sp3"

    epoch = conv.sp3name2dt(p)
    yeardec = conv.dt2year_decimal(epoch)

    yeardec
    Out: 2019.739611872146

It handles also RINEX names with ``conv.rinexname2dt(p_rinex)`` (old and new naming)
(Documentation here: :py:func:`geodezyx.conv.conv_time.rinexname2dt`)

-------------------
Convert coordinates 
-------------------

The `GeodeZYX Toolbox` can easily handle coordinate conversion in Geocentric (X,Y,Z), Geographic (latitude, longitude, height) and topocentric (East, North, Up). 

**Warning**: This is considered as the "low-level" coordinate conversions. It does not deal with the different Reference Frame and their realisations (ITRFxx, ETRFxx...). This is managed by "high-level" functions in the ``reffram``module.

These functions are optimized for arrays (multiple inputs) but can also handle scalars.

Coordinates conversion exemples
===============================

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
==================

All the module's functionalites can be found here:
:py:mod:`geodezyx.conv.conv_coords`

------------------------------
Read and import geodetic files 
------------------------------

------------------------------------------
Point and Click to detect offsets manually
------------------------------------------




