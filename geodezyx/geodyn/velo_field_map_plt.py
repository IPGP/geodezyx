# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.geodyn contains functions to plot velocity 
vectors on a map.

it can be imported directly with:
from geodezyx import geodyn

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""



########## BEGIN IMPORT ##########
#### External modules
import math

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

try:
    import netCDF4
    from mpl_toolkits.basemap import Basemap
except:
    pass

from matplotlib.patches import Ellipse

#### geodeZYX modules
from geodezyx.utils_xtra import plot_utils

#### Import star style
# from geodezyx import *                   # Import the GeodeZYX modules
# from geodezyx.externlib import *         # Import the external modules
# from geodezyx.megalib.megalib import *   # Import the legacy modules names


#### Import the logger
import logging
log = logging.getLogger('geodezyx')

##########  END IMPORT  ##########





# ------------------
# Plot land mask
# ------------------

def landmask(M, color='0.8'):
   # Make a constant colormap, default = grey
   constmap = matplotlib.colors.ListedColormap([color])

   jmax, imax = M.shape
   # X and Y give the grid cell boundaries,
   # one more than number of grid cells + 1
   # half integers (grid cell centers are integers)
   X = -0.5 + np.arange(imax+1)
   Y = -0.5 + np.arange(jmax+1)

   # Draw the mask by pcolor
   M = np.ma.masked_where(M > 0, M)
   plt.pcolor(X, Y, M, shading='flat', cmap=constmap)

# -------------
# Colormap
# -------------

# Colormap, smlgn. med Rob Hetland

def LevelColormap(levels, cmap=None):
    """
    Make a colormap based on an increasing sequence of levels
    """
    
    # Start with an existing colormap
    if cmap == None:
        cmap = plt.get_cmap()

    # Spread the colours maximally
    nlev = len(levels)
    S = np.arange(nlev, dtype='float')/(nlev-1)
    A = cmap(S)

    # Normalize the levels to interval [0,1]
    levels = np.array(levels, dtype='float')
    L = (levels-levels[0])/(levels[-1]-levels[0])

    # Make the colour dictionary
    R = [(L[i], A[i,0], A[i,0]) for i in range(nlev)]
    G = [(L[i], A[i,1], A[i,1]) for i in range(nlev)]
    B = [(L[i], A[i,2], A[i,2]) for i in range(nlev)]
    cdict = dict(red=tuple(R),green=tuple(G),blue=tuple(B))

    # Use 
    return matplotlib.colors.LinearSegmentedColormap(
        '%s_levels' % cmap.name, cdict, 256)


def draw_map(station_etude,latm,latM,lonm,lonM,path,
             all_pos,hw,vn_ITRF,ve_ITRF,plot_vertical_ITRF,
             incvn_ITRF,incve_ITRF,incplot_vertical_ITRF,
             plot_GEBCO=False,
             plot_vertical=False,
             plot_topo=True,
             plot_ellipses=True,
             coarse_lines=False,
             legend_arrow_length=("1 cm/yr",0.01),
             legend_ellipse_size=("2 mm/yr",0.002),
             legend_position=(0.5,0.9),
             scale_arrow=15000000,
             scale_ellipse=10000000,
             name_stats=True,
             name_stats_font_size=8,
             name_stats_offset=(0.005,0.01),
             shorten_oversized_arrows=True,
             exclude_points_out_of_range = True,
             adjust_text=False,
             pixels_hires_backgrnd=2000,
             draw_borders=True,
             full_return=False,
             draw_latlon_lines=True):
    
    nstation = len(station_etude)
    fig,ax=plt.subplots(figsize=(7,8))  
     
    if incvn_ITRF is None:
        incvn_ITRF = np.zeros(len(station_etude))
    if incve_ITRF is None:
        incve_ITRF = np.zeros(len(station_etude))
    if incplot_vertical_ITRF is None:
        incplot_vertical_ITRF = np.zeros(len(station_etude))

    if coarse_lines:
        resolution="l"
    else:
        resolution="h"
    m= Basemap(projection='merc', llcrnrlat=latm,urcrnrlat = latM, 
               llcrnrlon = lonm, urcrnrlon = lonM, lat_ts=-20, 
               resolution=resolution,
               area_thresh = 1,epsg=3395)
    
    if draw_borders:
        m.drawcoastlines(linewidth=0.25)    
        m.drawcountries()
        m.drawstates()

    if plot_topo:
        m.arcgisimage(service='World_Shaded_Relief', 
                      xpixels = pixels_hires_backgrnd, verbose= True)
    else:
        m.drawmapboundary(fill_color='#97B6E1')
        m.fillcontinents(color='#EFEFDB', lake_color='#97B6E1',zorder=1)
        m.shadedrelief()
    
    if not plot_GEBCO:
        color1='black'
    else:
        color1='white'
        levels = [-8000,-7000,-6000,-5500,-5000,-4500,-4000] #Niveaux pour l'echelle de couleur
        gebconc = netCDF4.Dataset(path+'GEBCO\\GRIDONE_2D.nc')    
        #Vecteur longitudes/latitudes
        long_ini = gebconc.variables['lon'][:]
        lats_ini = gebconc.variables['lat'][:]
        (lats_area,longs_area,gebco_area) = area(long_ini,lats_ini,gebconc,latm,latM,lonm,lonM,'0-180') #lat bas,lat haut, long
        (longs,lats,gebco) = split_grid(longs_area,lats_area,gebco_area,2)
        (longs_gr,lats_gr)=np.meshgrid(longs,lats)  
        fond = m.pcolormesh(longs_gr,lats_gr,
                            gebco,shading='flat',
                            cmap=LevelColormap(list(np.asarray(levels)*(1)),
                                               cmap='deep_r'), # cmap=cm.deep_r
                            latlon=True) #ice,deep
        cbar=m.colorbar(location='top',pad=0.5)
        cbar.set_label('Depth [m] ', rotation=0)
        
    if draw_latlon_lines:
        m.drawparallels(np.arange(latm,latM,1.),fontsize=10,color=color1,labels=[True,False,False,False],linewidth=0.25, dashes=[10000,1]) #de -180 a 360 par pas de 5° ;   
        m.drawmeridians(np.arange(lonm,lonM,1.),fontsize=10,color=color1,labels=[False,False,False,True],linewidth=0.25, dashes=[10000,1]) # Haut / G/D/BAS 
        ### Labels
        plt.xlabel('Longitude (°) ',labelpad=25,fontsize=10)
        plt.ylabel('Latitude (°) ',labelpad=40,fontsize=10)
        ax.yaxis.set_label_position("left")
       
    all_posx_proj=[0]*nstation
    all_posy_proj=[0]*nstation
    for i in range(nstation):
          all_posx_proj[i],all_posy_proj[i] = m(all_pos[i][1],all_pos[i][0]) #positions xy des stations converties a la projection de la carte  [m]         
          ## print(i,all_posx_proj[i],all_posy_proj[i])
          m.plot(all_posx_proj[i],all_posy_proj[i],'xkcd:black',marker='.',linestyle='none',markersize=4,lw=4,zorder=20) #point rouge à la station, zorder indique l'ordre vertical de dessin sur la carte, les nb elevés sont dessinés les plus hauts
#          plt.text(all_posx_proj[i], all_posy_proj[i], station_etude[i], size=10,ha="center", va="center",bbox=dict(boxstyle="round",
#                           ec=(1., 0.5, 0.5),
#                           fc=(1., 0.8, 0.8),
#                           )
#                 )


    ############### DICO POUR FLECHES DEFINITION
    arrow_prop_dict = dict(arrowstyle='->,head_length=0.8,head_width=0.3',
                           ec='xkcd:black',fc='xkcd:black',linewidth=2)

    arrow_prop_dict_V_up = dict(arrowstyle='->,head_length=0.8,head_width=0.3',
                           ec='xkcd:green',fc='xkcd:green',linewidth=2)

    arrow_prop_dict_V_down = dict(arrowstyle='->,head_length=0.8,head_width=0.3',
                           ec='xkcd:red orange',fc='xkcd:red orange',linewidth=2)

    arrow_prop_dict_V_up_short = dict(arrowstyle='->,head_length=0.8,head_width=0.3',
                           ec='xkcd:green',fc='xkcd:green',linewidth=2,
                           linestyle="--")

    arrow_prop_dict_V_down_short = dict(arrowstyle='->,head_length=0.8,head_width=0.3',
                           ec='xkcd:red orange',fc='xkcd:red orange',linewidth=2,
                           linestyle="--")

    if not plot_vertical: #CHAMP VITESSES HORIZONTALES        
        for i in range(nstation):
            
            ######### Exclude point if not in range
            #latm,latM,lonm,lonM
            bool_good_lat = latm < all_pos[i][0] < latM
            bool_good_lon = lonm < all_pos[i][1] < lonM
            if exclude_points_out_of_range and not (bool_good_lat and bool_good_lon):
                log.info("INFO : exclude point because out of range")
                continue
        
            
            x_end_arrow_ok = all_posx_proj[i]+np.multiply(ve_ITRF[i],scale_arrow)
            y_end_arrow_ok = all_posy_proj[i]+np.multiply(vn_ITRF[i],scale_arrow)
            
            ax.annotate(text='', 
                         xy=(x_end_arrow_ok,
                             y_end_arrow_ok), 
                         xytext=(all_posx_proj[i], all_posy_proj[i]) , 
                         arrowprops=arrow_prop_dict,
                         annotation_clip=True) #xy point arrivee de la fleche, xytext l'origine de la fleche
            
            a=math.atan((all_posx_proj[i]+np.multiply(ve_ITRF[i],scale_arrow)-all_posy_proj[i])/(all_posy_proj[i]+np.multiply(ve_ITRF[i],scale_arrow)-all_posy_proj[i])) # l'angle d'orientation de l'ellipsoide, avec 0 au nord et +90 a l'ouest
                
            e = Ellipse(xy=(x_end_arrow_ok,
                            y_end_arrow_ok), 
                            width=np.multiply(2,incve_ITRF[i])*scale_ellipse,
                            height=np.multiply(2,incvn_ITRF[i])*scale_ellipse,
                            angle=a) #multiplication par 2 pour obtenir la largeur et la hauteur de l'ellipse 
                
            # STATION NAME
            if name_stats:
                offset_x_ok , offset_y_ok = plot_utils.axis_data_coords_sys_transform(ax,
                                                                              name_stats_offset[0],
                                                                              name_stats_offset[1])
                plt.text(all_posx_proj[i] + offset_x_ok,
                         all_posy_proj[i] + offset_y_ok,
                         station_etude[i], fontsize=name_stats_font_size)           
            
    
            # ELIPSES
            if plot_ellipses:
                ax.add_artist(e)
                e.set_clip_box(ax.bbox)
                e.set_edgecolor('none')
                e.set_facecolor('black')
                e.set_alpha(0.3)
                e.set_zorder(1) 
            
    else:# Champ de vitesses verticales ITRF	   
        ############### PLOT FLECHES
        Text = []
        for i in range(nstation): 
            
            ######### Exclude point if not in range
            #latm,latM,lonm,lonM
            bool_good_lat = latm < all_pos[i][0] < latM
            bool_good_lon = lonm < all_pos[i][1] < lonM
            if exclude_points_out_of_range and not (bool_good_lat and bool_good_lon):
                log.info("INFO : exclude point because out of range")
                continue
        
            ######### AJUSTEMENT SI FLECHES TROP GRANDES
            x_end_arrow , y_end_arrow = all_posx_proj[i],all_posy_proj[i]+np.multiply(plot_vertical_ITRF[i],scale_arrow)
            
            x_end_axis_ref , y_end_axis_ref = plot_utils.axis_data_coords_sys_transform(ax,
                                                              x_end_arrow,
                                                              y_end_arrow,
                                                              inverse=True) 
            if (y_end_axis_ref < 0.) and shorten_oversized_arrows:
                shortened_arrow = True
                x_end_arrow_ok = x_end_arrow
                _ , y_end_arrow_ok = plot_utils.axis_data_coords_sys_transform(ax,
                                                              x_end_axis_ref,
                                                              0,
                                                              inverse=False)
            elif (y_end_axis_ref > 1.) and shorten_oversized_arrows:
                shortened_arrow = True
                x_end_arrow_ok = x_end_arrow
                _ , y_end_arrow_ok = plot_utils.axis_data_coords_sys_transform(ax,
                                                              x_end_axis_ref,
                                                              1.,
                                                              inverse=False) 
            else:
                shortened_arrow = False
                x_end_arrow_ok = x_end_arrow
                y_end_arrow_ok = y_end_arrow
                
            ######### FIN AJUSTEMENT SI FLECHES TROP GRANDES

            ### SELECT ARROW PROPERTY DEPENDING ON THE TYPE
            if plot_vertical_ITRF[i]>0 and not shortened_arrow:
                arrow_prop_dict_V = arrow_prop_dict_V_up
            elif plot_vertical_ITRF[i]>0 and shortened_arrow:
                arrow_prop_dict_V = arrow_prop_dict_V_up_short
            elif plot_vertical_ITRF[i]<=0 and not shortened_arrow:
                arrow_prop_dict_V = arrow_prop_dict_V_down
            elif plot_vertical_ITRF[i]<=0 and shortened_arrow:
                arrow_prop_dict_V = arrow_prop_dict_V_down_short
                                
            ### PLOT      
            ax.annotate(text='',
                      xy=(x_end_arrow_ok , y_end_arrow_ok),
                      xytext=(all_posx_proj[i], all_posy_proj[i]) , 
                      arrowprops=arrow_prop_dict_V,
                      annotation_clip=False)
            #xy point arrivee de la fleche, xytext l'origine de la fleche
                
            ### STATION NAME
            if name_stats:
                offset_x_ok , offset_y_ok = plot_utils.axis_data_coords_sys_transform(ax,
                                                                              name_stats_offset[0],
                                                                              name_stats_offset[1])
                Text.append(plt.text(all_posx_proj[i] + offset_x_ok,
                         all_posy_proj[i] + offset_y_ok,
                         station_etude[i], fontsize=name_stats_font_size))
                


            ############### PLOT ELLIPSES
            # angle de l'ellipse dépend de la direction du vecteur vitesse            
            a=math.atan((all_posx_proj[i]+np.multiply(plot_vertical_ITRF[i],scale_arrow)-all_posx_proj[i])/(all_posy_proj[i]+np.multiply(plot_vertical_ITRF[i],scale_arrow)-all_posy_proj[i])) # l'angle d'orientation de l'ellipsoide, avec 0 au nord et +90 a l'ouest
            
            # weird previous xy definition
            # [all_posx_proj[i]+np.multiply(plot_vertical_ITRF[i],0),all_posy_proj[i]+0.85*np.multiply(plot_vertical_ITRF[i],scale_arrow)]
            
            e = Ellipse(xy=(x_end_arrow_ok , y_end_arrow_ok), 
                        width=np.multiply(2,incplot_vertical_ITRF[i])*scale_ellipse,
                        height=np.multiply(2,incplot_vertical_ITRF[i])*scale_ellipse, angle=a) #multiplication par 2 pour obtenir la largeur et la hauteur de l'ellipse 
            # ax.plot(all_posx_proj[i]+np.multiply(plot_vertical_ITRF[i],0),all_posy_proj[i]+np.multiply(plot_vertical_ITRF[i],scale_arrow),'xkcd:green',marker='.',linestyle='none',markersize=4,lw=4,zorder=20) #point rouge à la station, zorder indique l'ordre vertical de dessin sur la carte, les nb elevés sont dessinés les plus hauts
      
            ax.add_artist(e)
            e.set_clip_box(ax.bbox)
            e.set_edgecolor('none')
            e.set_facecolor('black')
            e.set_alpha(0.3)
            e.set_zorder(1)       
            
    ############### LEGENDE
    legend_arrow_length_label  = legend_arrow_length[0]
    legend_ellipse_size_label  = legend_ellipse_size[0]
    legend_arrow_length_metric = legend_arrow_length[1]
    legend_ellipse_size_metric = legend_ellipse_size[1]
    legend_position_x = legend_position[0]
    legend_position_y = legend_position[1]
    
    ######## Fleche de Legende    
    ### origine de la fleche de legende        
    xl,yl = plot_utils.axis_data_coords_sys_transform(ax,legend_position_x,
                                              legend_position_y + 0.02,
                                              inverse=False)
    ### fin de la fleche de legende                
    xe,ye = xl+np.multiply(legend_arrow_length_metric,scale_arrow),yl
    
    plt.annotate(text='', xy=(xe,ye), xytext=(xl, yl) , 
                 arrowprops=arrow_prop_dict) #xy point arrivee de la fleche, xytext l'origine de la fleche

    props = dict(boxstyle='round', facecolor='black', alpha=0)
    ax.text(legend_position_x, legend_position_y + 0.035, 
            'Velocity : ' + legend_arrow_length_label,
            transform=ax.transAxes, fontsize=11,color='black',
        verticalalignment='bottom', bbox=props)# place a text box in upper left in axes coords

    ######## Légende de l'ellipse d'incertitude
    ax.text(legend_position_x, legend_position_y,
            '+/- ' + legend_ellipse_size_label + ' (radius)',
            transform=ax.transAxes, fontsize=11,color='black',
    verticalalignment='top', bbox=props)
    ells_legend = Ellipse(xy=[xe,ye],
                          width=np.multiply(2,np.multiply(legend_ellipse_size_metric,scale_ellipse)), 
                          height=np.multiply(2,np.multiply(legend_ellipse_size_metric,scale_ellipse)))
                          #angle=a)
        
    if plot_ellipses:
        ax.add_artist(ells_legend)
        ells_legend.set_clip_box(ax.bbox)
        ells_legend.set_edgecolor('none')
        ells_legend.set_facecolor('black')
        ells_legend.set_alpha(0.3)
        ells_legend.set_zorder(1)  

    ############### FINITION
    
    for T in Text:   #Put label in front
        T.set_zorder(100)

    if adjust_text:
        from adjustText import adjust_text
        adjust_text(Text,on_basemap=True) #arrowprops=dict(arrowstyle='->', color='red')
        

    plt.tight_layout() 
        
    if full_return:
        return fig , ax , m , Text
    else:
        return fig , ax , m , Text


#
#
#
#    
#def draw_map(path,stat,lon,lat,statout,latm,latM,lonm,lonM,stat_res,GEBCO,OUTSTAT,RES): 
#    mscamp= 8#taille des marqueurs stations de campagne
#    ml=0.4 #taille des bordures noires des marqueurs, stations campagne
#    msper=12 #taille marqueurs stations permanentes
#    mlp=1.7 #taille des bordures noires des marqueurs, stations campagne
#    
#    fig,ax=plt.subplots(figsize=(12,8))       
#    m= Basemap(projection='merc', llcrnrlat=latm,urcrnrlat = latM, llcrnrlon = lonm, urcrnrlon = lonM, lat_ts=-20, resolution='h', area_thresh = 1)
#    m.drawcoastlines(linewidth=0.25)
##    m.drawcountries(linewidth=0.25)
#    m.drawmapboundary()
#    
#    if GEBCO==0:
#        m.fillcontinents(color='#bfa367', lake_color='white',zorder=1)
#        color1='black'
#        m.drawmapboundary(fill_color='#edf1f7')
#
#    if GEBCO==1:
#        color1='white'
#        m.fillcontinents(color='#bfa367', lake_color='white',zorder=1)
#
#        levels = [-8000,-7000,-6000,-5500,-5000,-4500,-4000]#Niveaux pour l'echelle de couleur
#        gebconc = Dataset(path+'GEBCO\\GRIDONE_2D.nc')    
#        #Vecteur longitudes/latitudes
#        long_ini = gebconc.variables['lon'][:]
#        lats_ini = gebconc.variables['lat'][:]
#        (lats_area,longs_area,gebco_area) = area(long_ini,lats_ini,gebconc,latm,latM,lonm,lonM,'0-180') #lat bas,lat haut, long
#        (longs,lats,gebco) = split_grid(longs_area,lats_area,gebco_area,2)
#        (longs_gr,lats_gr)=np.meshgrid(longs,lats)  
#        fond = m.pcolormesh(longs_gr,lats_gr,gebco,shading=geodyn.LevelColormap(list(np.asarray(levels)*(1)),cmap=cm.deep_r),latlon=True) #ice,deep
#        cbar=m.colorbar(location='top',pad=0.5)
#        cbar.set_label('Depth [m] ', rotation=0)
#    m.drawparallels(np.arange(latm,latM,10.),fontsize=10,color=color1,labels=[True,False,False,False],linewidth=0.25, dashes=[10000,1]) #de -180 a 360 par pas de 5° ;   
#    m.drawmeridians(np.arange(lonm,lonM,40.),fontsize=10,color=color1,labels=[False,False,False,True],linewidth=0.25, dashes=[10000,1]) # Haut / G/D/BAS 
#    
#    ### Labels
#    plt.xlabel('Longitude (°) ',labelpad=25,fontsize=10)
#    plt.ylabel('Latitude (°) ',labelpad=40,fontsize=10)
#    ax.yaxis.set_label_position("left")
#    all_posx_proj,all_posy_proj = m(lon,lat) #positions xy des stations converties a la projection de la carte  [m]         
#    
#    if RES==0:
#        for i in range(len(all_posx_proj)-1,len(all_posx_proj)):
#              if stat[i] in ['KOUC',  'NOUM', 'LPIL', 'VILA' , 'TGWS' , 'SANC','ESPI', 'SSAN', 'Aplot_verticalN', 'AMBA', 'VATU', 'MAWO', 'VANU', 'DVIL', 'HGHN' , 'THIO']:
#                  m.plot(all_posx_proj[i],all_posy_proj[i],'xkcd:red orange',marker='^',markeredgecolor='black',markeredgewidth=mlp,linestyle='none',markersize=msper,lw=6,zorder=3,label='Permanent stations') #point rouge à la station, zorder indique l'ordre vertical de dessin sur la carte, les nb elevés sont dessinés les plus hauts   
#              m.plot(all_posx_proj[i],all_posy_proj[i],'xkcd:wheat',marker='^',markeredgecolor='black',markeredgewidth=ml,linestyle='none',markersize=mscamp,lw=6,zorder=2,label='Campaign stations') #point rouge à la station, zorder indique l'ordre vertical de dessin sur la carte, les nb elevés sont dessinés les plus hauts   
#              
#              if OUTSTAT==1:    
#                  if stat[i] in [statout]:
#                        m.plot(all_posx_proj[i],all_posy_proj[i],'xkcd:hotpurple',marker='^',markeredgecolor='black',markeredgewidth=ml,linestyle='none',markersize=msper,lw=6,zorder=4,label=statout) #point rouge à la station, zorder indique l'ordre vertical de dessin sur la carte, les nb elevés sont dessinés les plus hauts 
#                                  
#                  
#        for i in range(len(all_posx_proj)-1):
#              m.plot(all_posx_proj[i],all_posy_proj[i],'xkcd:wheat',marker='^',markeredgecolor='black',markeredgewidth=ml,linestyle='none',markersize=mscamp,lw=6,zorder=2) #point rouge à la station, zorder indique l'ordre vertical de dessin sur la carte, les nb elevés sont dessinés les plus hauts   
#              if stat[i] in ['KOUC',  'NOUM', 'LPIL', 'VILA' , 'TGWS' , 'SANC','ESPI', 'SSAN', 'Aplot_verticalN', 'AMBA', 'VATU', 'MAWO', 'VANU', 'DVIL', 'HGHN' , 'THIO']:
#                  m.plot(all_posx_proj[i],all_posy_proj[i],'xkcd:red orange',marker='^',markeredgecolor='black',markeredgewidth=mlp,linestyle='none',markersize=msper,lw=6,zorder=3) #point rouge à la station, zorder indique l'ordre vertical de dessin sur la carte, les nb elevés sont dessinés les plus hauts   
#              if OUTSTAT==1:    
#                  if stat[i] in [statout]:
#                        m.plot(all_posx_proj[i],all_posy_proj[i],'xkcd:hot purple',marker='^',markeredgecolor='black',markeredgewidth=mlp,linestyle='none',markersize=msper,lw=6,zorder=4,label=statout) #point rouge à la station, zorder indique l'ordre vertical de dessin sur la carte, les nb elevés sont dessinés les plus hauts 
##              plt.text(all_posx_proj[i], all_posy_proj[i], stat[i], size=10,ha="center", va="center",bbox=dict(boxstyle="round",
##                       ec=(1., 0.5, 0.5),
##                       fc=(1., 0.8, 0.8),
##                       )
##             )
#    
#    if RES==1:
#         for i in range(len(all_posx_proj)-1,len(all_posx_proj)):
##              if stat[i] in ['KOUC',  'NOUM', 'LPIL', 'VILA' , 'TGWS' , 'SANC','ESPI', 'SSAN', 'Aplot_verticalN', 'AMBA', 'VATU', 'MAWO', 'VANU', 'DVIL', 'HGHN' , 'THIO']:
#              if stat[i] in stat_res:
#                  m.plot(all_posx_proj[i],all_posy_proj[i],'xkcd:red orange',marker='^',markeredgecolor='black',markeredgewidth=mlp,linestyle='none',markersize=msper,lw=6,zorder=3,label='ulr6 Network') #point rouge à la station, zorder indique l'ordre vertical de dessin sur la carte, les nb elevés sont dessinés les plus hauts   
##              m.plot(all_posx_proj[i],all_posy_proj[i],'xkcd:wheat',marker='^',markeredgecolor='black',markeredgewidth=ml,linestyle='none',markersize=mscamp,lw=6,zorder=2,label='ulr6') #point rouge à la station, zorder indique l'ordre vertical de dessin sur la carte, les nb elevés sont dessinés les plus hauts   
#                  
#         for i in range(len(all_posx_proj)-1):
##              m.plot(all_posx_proj[i],all_posy_proj[i],'xkcd:wheat',marker='^',markeredgecolor='black',markeredgewidth=ml,linestyle='none',markersize=mscamp,lw=6,zorder=2) #point rouge à la station, zorder indique l'ordre vertical de dessin sur la carte, les nb elevés sont dessinés les plus hauts   
#              if stat[i] in stat_res:
#                  m.plot(all_posx_proj[i],all_posy_proj[i],'xkcd:red orange',marker='^',markeredgecolor='black',markeredgewidth=mlp,linestyle='none',markersize=msper,lw=6,zorder=3) #point rouge à la station, zorder indique l'ordre vertical de dessin sur la carte, les nb elevés sont dessinés les plus hauts   
#             
#                
#    pylab.legend(loc='lower left')               
##    plt.title('Mercator projection, bathymetry from GEBCO')
#    fig.tight_layout()           
#    return fig,m


# fonctions pour plot gebco
def area(lo,la,grid,minla,maxla,minlo,maxlo,type_long='0-180'):
    if type_long == '0-360':

        #index des bornes - grille 1 
        ilo_min = lo.tolist().index(minlo)
        ilo_max = lo.tolist().index(180)
        ila_min = la.tolist().index(minla)
        ila_max = la.tolist().index(maxla)
    
        #Nouvelle grille 1 et vecteurs de coordonnees
        newgrid1 = grid.variables['elevation'][ila_min:ila_max,ilo_min:ilo_max]
        
        #index des bornes - grille 2 
        ilo_min = lo.tolist().index(-180)
        ilo_max = lo.tolist().index(maxlo - 360)
        ila_min = la.tolist().index(minla)
        ila_max = la.tolist().index(maxla)
    
        #Nouvelle grille 2 et vecteurs de coordonnees
        newgrid2 = grid.variables['elevation'][ila_min:ila_max,ilo_min:ilo_max]
        
        #grille finale
        newgrid = np.concatenate((newgrid1,newgrid2),axis=1)
        lo += 180
        ilo_min = lo.tolist().index(minlo)
        ilo_max = lo.tolist().index(maxlo)
        
        longs = lo[ilo_min:ilo_max]
        lats = grid.variables['lat'][ila_min:ila_max]
        
    else:    
        #index des bornes
        ilo_min = lo.tolist().index(minlo)
        ilo_max = lo.tolist().index(maxlo)
        ila_min = la.tolist().index(minla)
        ila_max = la.tolist().index(maxla)
    
        #Nouvelle grille et vecteurs de coordonnees
        newgrid = grid.variables['elevation'][ila_min:ila_max,ilo_min:ilo_max]
        longs = grid.variables['lon'][ilo_min:ilo_max]
        lats = grid.variables['lat'][ila_min:ila_max]
    
    return (lats,longs,newgrid)
    
def split_grid(lo,la,gr,fact):
    
    nblo = len(lo)
    nbla = len(la)
    nblo_new = int(np.floor(nblo / float(fact)))
    nbla_new = int(np.floor(nbla / float(fact)))
    
    new_grid = np.zeros((nbla_new, nblo_new))
    new_lo = []
    new_la = []
    

    for i in range(0,nbla_new):
            new_la.append(la[i*fact])
            for j in range(0,nblo_new): 
                new_grid[i,j] = gr[i*fact,j*fact]
    
    for k in range(0,nblo_new):
        new_lo.append(lo[k*fact])
        
    return (new_lo,new_la,new_grid)



########################################################################

