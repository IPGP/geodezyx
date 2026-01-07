
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 30/12/2025 21:56:47

@author: psakic
"""

import conv_time
import datetime as dt

tm = dt.datetime(2020,1,5,12,0,0)
tl = conv_time.dt_range(dt.datetime(2020,1,1), dt.datetime(2020,1,10))

a = conv_time.dt2gpstime(tl)
