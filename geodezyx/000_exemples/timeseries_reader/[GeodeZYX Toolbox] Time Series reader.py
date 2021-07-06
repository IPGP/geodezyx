#!/usr/bin/env python
# coding: utf-8

# ### Import the librairies

# In[1]:


#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names


# ### Read TimeSeries

# We take Nevada Geodetic Laboratory's Time Series as exemple. But the GeodeZYX Toolbox can handle a lot of formats. See the documentation for more details.

# In[94]:


#### file path of the ENU coordinates
exemple_dir = "/home/psakicki/CODES/GeodeZYX-Toolbox_v4/geodezyx/000_exemples/timeseries_reader/exemples"
p_nevada_xyz = exemple_dir + "/POTS.txyz2"
#### Import the time series
TS = files_rw.read_nevada(p_nevada_xyz,"xyz")


# Data are imported in a TimeSeries object. Let's see what we can do with it.

# ### Access the coordinates

# #### Get the number of points

# In[59]:


TS.nbpts


# #### Get the mean position

# In[63]:


TS.mean_posi()


# #### Get a random point

# In[64]:


TS.aleapt()


# #### As lists: X,Y,Z,T in posix time, sigmaX, sigmaY, sigmaZ

# In[54]:


TS.to_list()


# #### As a DataFrame (recommended)

# In[65]:


TS.to_dataframe()


# #### As DataFrame, but with geographic coordinates

# In[66]:


TS.to_dataframe('FLH') ### F: Phi, L:Lambda, H: Height


# ### Determine the East North Up coordinates

# #### w.r.t. mean position

# In[67]:


TS.ENUcalc_from_mean_posi()


# #### w.r.t. a ref position

# In[68]:


pt_ref = time_series.Point(3800689.6341,882077.3857,5028791.3179 ,2008.0,"XYZ")
TS.ENUcalc(pt_ref)


# #### Now we can access to the ENU coordinates

# In[69]:


TS.to_dataframe('ENU') ### F: Phi, L:Lambda, H: Height


# ### Plot The TimeSeries

# In[37]:


TS.plot('ENU')


# ### Clean the outliers from the TimeSeries

# #### Based on high sigmas

# In[52]:


std_val_lim = 0.005 ### exclude point with std higher than 5mm (exemple)
coords_type = 'ENU'
TSclean = time_series.std_dev_cleaner(TS,std_val_lim,coords_type,verbose=True)


# #### Based on MAD detection

# In[71]:


TSclean = time_series.mad_cleaner(TSclean,coortype=coords_type,detrend_first=True,verbose=True)


# ### Shorten a TimeSeries

# We keep 2 windows: 2000-2005 and 2010 until the last epoch

# In[78]:


TSshort = time_series.time_win(TSclean,[(dt.datetime(2000,1,1),dt.datetime(2005,12,31)),
                                        (dt.datetime(2010,1,1),dt.datetime(2099,12,1))]) 


# ### Add discontinuities

# #### Extract the discontinuities from the Sitelogs

# In[92]:


Ant_dico = dict()
Rec_dico = dict()

ls = exemple_dir + '/pots_20210215.log'

#### Read the Site name, Receiver blocks, and Antenna Blocks
Site     = gfc.read_blocks_logsheet(ls,1)
Rec_list = gfc.read_blocks_logsheet(ls,3)
Ant_list = gfc.read_blocks_logsheet(ls,4)
    
Rec_dico[Site[0].Four_Character_ID] = Rec_list
Ant_dico[Site[0].Four_Character_ID] = Ant_list 


# In[93]:


#### Extract the discontinuities epochs
Discont_Rec = [Rec.Date_Installed for Rec in Rec_dico[TSshort.stat]]
Discont_Ant = [Ant.Date_Installed for Ant in Ant_dico[TSshort.stat]]
        
Discont = sorted(Discont_Rec + Discont_Ant)


# #### Apply the discontinuities to the TimeSeries

# In[85]:


TSshort.set_discont(Discont)


# ### Plot the result

# In[89]:


TSshort.plot("ENU")
TSshort.plot_discont()


# ### Export

# #### Export content in an HECTOR's compatible NEU file

# In[101]:


export_path_hector = exemple_dir
time_series.export_ts_as_neu(TSshort, export_path_hector , "",coordtype=coords_type)


# #### Export content in an MIDAS's compatible TNEU file 

# In[102]:


export_path_midas = exemple_dir
time_series.export_ts_as_midas_tenu(TSshort, export_path_midas , "",coordtype=coords_type)


# In[ ]:




