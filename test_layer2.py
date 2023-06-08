from astropy.coordinates import Angle, SkyCoord, SkyOffsetFrame
from regions import RectangleSkyRegion, PixCoord, RectanglePixelRegion
from astropy.wcs import WCS
from astropy.io import fits
import matplotlib.pyplot as plt
from astroquery.skyview import SkyView
import astropy.units as u
from astropy.nddata import Cutout2D
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QMessageBox
import sys
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout, QHeaderView
from PyQt5.QtCore import pyqtSignal, pyqtSlot, QObject, QThread, QTimer
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as Navi
import sip
import pandas as pd
import numpy as np
import threading
import aplpy
from scipy import ndimage, misc
from astropy.visualization.wcsaxes import WCSAxes
from matplotlib.transforms import Affine2D


def rotate(p, origin=(0, 0), degrees=0):
    angle = np.deg2rad(degrees)
    R = np.array([[np.cos(angle), -np.sin(angle)],
                  [np.sin(angle),  np.cos(angle)]])
    o = np.atleast_2d(origin)
    p = np.atleast_2d(p)
    return np.squeeze((R @ (p.T-o.T) + o.T).T)

source = "M 51"
paths = SkyView.get_images(position=source, coordinates='J2000',
                           radius=10*u.arcmin, height=10*u.arcmin,
                        survey=['2MASS-K'])
coord = SkyCoord.from_name(source, frame = 'icrs')
hdu = paths[0][0]
wcs = WCS(hdu.header)



image_rotation = 20    #in degrees

rotated_array = ndimage.rotate(hdu.data, image_rotation, reshape=False)


#ax = fig.add_subplot(111,)
fig = plt.figure()


ax = fig.add_subplot(111,)
ax.imshow(rotated_array, origin='lower', vmin=472, vmax=534,
          cmap="gist_earth")



# Slit configuration
pixel_scale = 1/(hdu.header['CDELT2']*3600) # arcsec/pixel

fov_width = 3.1*60  # in arcsec
delta_fov = 9.1*60  # in arcsec


coord_fov = SkyCoord(coord.ra+(0/3600)*u.deg,
                     coord.dec+(0)*u.arcsec)

w_rotate_angle=60
sky_region = RectangleSkyRegion(coord_fov,
                                width=fov_width*u.arcsec,
                                height=delta_fov*u.arcsec,
                                angle=w_rotate_angle*u.deg)
pixel_region = sky_region.to_pixel(wcs)
artist = pixel_region.as_artist(color='w', lw=5)
#ax.add_artist(artist)

pix_slit_angle=0
reg = RectanglePixelRegion(PixCoord(x=150, y=150), width=fov_width*pixel_scale ,
                           height=delta_fov*pixel_scale , angle=Angle(pix_slit_angle, 'deg'))
patch = reg.plot(ax=ax, facecolor='none', edgecolor='red', lw=2,
                 label='Rectangle')

rot_matrix = np.array([[np.cos(pix_slit_angle*np.pi/180), -np.sin(pix_slit_angle*np.pi/180)],
                       [np.sin(pix_slit_angle*np.pi/180), np.cos(pix_slit_angle*np.pi/180)]])



c = ['r', 'tab:orange', 'yellow', 'lime', 'pink']
slit_fov_check = []




for i in range(5):
    slit0_width = 10
    delta_x0 = 5 * 1e3 * 1/18 ##in microns
    
    coord_slit0 = SkyCoord(coord.ra+(delta_x0/3600)*u.deg,
                           coord.dec+9.1*(2-i)/5*u.arcmin)
        
    ds = coord_slit0.transform_to(SkyOffsetFrame(origin=coord, 
                                            rotation=w_rotate_angle*u.deg))
    
    
    coord_slit0 = SkyCoord(coord_slit0.ra+ds.lon,
                           coord_slit0.dec+ds.lat)
    sky_region = RectangleSkyRegion(coord_slit0,
                                    width=slit0_width*u.arcsec,
                                    height=9.1/5*u.arcmin,
                                    angle=w_rotate_angle*u.deg)
    pixel_region = sky_region.to_pixel(wcs)
    artist = pixel_region.as_artist(color=c[i], lw=1)
    #ax.add_artist(artist)
    ax.text(0.1, 0.9-i*0.05, f'Slit {i}', transform=ax.transAxes,
            c=c[i],  weight="bold")
    
    delta_x0 = 5 * 1e3 * 1/18 ## converted to pixels from mm
    
    xold = 150+delta_x0
    yold = 150+9.1*60/5*pixel_scale*(2-i)
    points = [(xold, yold)]
    origins = (150, 150)
    
    new_points = rotate(points, origin=origins, degrees=pix_slit_angle)
    xnew, ynew = new_points[0], new_points[1]
    
    
    reg = RectanglePixelRegion(PixCoord(x=xnew, y=ynew),
                               width=slit0_width*pixel_scale,
                               height=9.1*60/5*pixel_scale , angle=Angle(pix_slit_angle, 'deg'))
    patch = reg.plot(ax=ax, facecolor='none', edgecolor=c[i], lw=2,
                     label='Rectangle')
    if abs(delta_x0) > fov_width/2:
        slit_fov_check.append(delta_x0)

