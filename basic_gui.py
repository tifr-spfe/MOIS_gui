#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 19:43:16 2022

@author: black_neko
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 00:09:39 2022

@author: Himanshu Tyagi
"""
# Form implementation generated from reading ui file 'sign_up.ui'
#
# Created by: PyQt5 UI code generator 5.15.7
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.




from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle, SkyCoord
from regions import RectangleSkyRegion
from astropy.wcs import WCS
from astropy.io import fits
#from astropy.utils.data import get_pkg_data_filename
import matplotlib.pyplot as plt
from astroquery.skyview import SkyView
import astropy.units as u
from astropy.nddata import Cutout2D
from PyQt5 import QtCore, QtGui, QtWidgets
import matplotlib.pyplot as plt
import sys
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout, QHeaderView
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as Navi
from matplotlib.figure import Figure
import platform

from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc

from PySide6.QtCore import Qt, Slot
from PySide6.QtWidgets import (
    QApplication,
    QWidget,
    QDoubleSpinBox,
    QVBoxLayout,
    QHBoxLayout,
)

import sip

class MatplotlibCanvas(FigureCanvasQTAgg):
	def __init__(self,parent=None, dpi = 120):
		self.fig = Figure(figsize=[10,10], dpi = dpi)
		#self.axes = fig.add_subplot(111)

		super(MatplotlibCanvas,self).__init__(self.fig)
		#self.fig.set_positions(0.2, 0.2, 0.8, 0.8) # left,bottom,width,height
        #fig.tight_layout()
class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1535, 1266)
        MainWindow.setFocusPolicy(QtCore.Qt.TabFocus)
        MainWindow.setStyleSheet("\n"
"\n"
"background-color: rgb(255,255,255);\n"
"font: 14pt \"Times\";")
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        """
        self.pushButton = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton.setGeometry(QtCore.QRect(280, 300, 106, 30))
        self.pushButton.setObjectName("pushButton")
        """ 
        self.textEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.textEdit.setGeometry(QtCore.QRect(280, 110, 371, 41))
        self.textEdit.setObjectName("textEdit")
        
        self.textEdit_filename = QtWidgets.QLineEdit(self.centralwidget)
        self.textEdit_filename.setGeometry(QtCore.QRect(280, 180, 371, 41))
        self.textEdit_filename.setObjectName("textEdit_filename")
        
        """
        self.textEdit_2 = QtWidgets.QTextEdit(self.centralwidget)
        self.textEdit_2.setGeometry(QtCore.QRect(280, 170, 371, 41))
        self.textEdit_2.setObjectName("textEdit_2")
        """
        """
        self.plainTextEdit = QtWidgets.QPlainTextEdit(self.centralwidget)
        self.plainTextEdit.setGeometry(QtCore.QRect(280, 230, 371, 41))
        self.plainTextEdit.setObjectName("plainTextEdit")
        """
        
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(120, 110, 131, 31))
        self.label.setObjectName("label")
        
        self.label_filename = QtWidgets.QLabel(self.centralwidget)
        self.label_filename.setGeometry(QtCore.QRect(120, 180, 131, 31))
        self.label_filename.setObjectName("label")
        
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(120, 250, 180, 120))
        self.label_2.setObjectName("label_2")
        """
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(120, 230, 131, 31))
        self.label_3.setObjectName("label_3")
        
        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(190, 510, 131, 31))
        self.label_4.setObjectName("label_4")
        
        self.plainTextEdit_2 = QtWidgets.QPlainTextEdit(self.centralwidget)
        self.plainTextEdit_2.setGeometry(QtCore.QRect(340, 500, 371, 41))
        self.plainTextEdit_2.setObjectName("plainTextEdit_2")
        self.plainTextEdit_3 = QtWidgets.QPlainTextEdit(self.centralwidget)
        self.plainTextEdit_3.setGeometry(QtCore.QRect(340, 550, 371, 41))
        self.plainTextEdit_3.setObjectName("plainTextEdit_3")
        """

        self.label_Temp = QtWidgets.QLabel(self.centralwidget)
        self.label_Temp.setGeometry(QtCore.QRect(120, 565, 131, 31))
        self.label_Temp.setObjectName("Temperature")
        

        
        self.lcdNumber = QtWidgets.QLCDNumber(self.centralwidget)
        self.lcdNumber.setGeometry(QtCore.QRect(300, 560, 161, 61))
        self.lcdNumber.setFocusPolicy(QtCore.Qt.NoFocus)
        self.lcdNumber.setStyleSheet("color: rgb(255, 0, 0);")
        self.lcdNumber.setInputMethodHints(QtCore.Qt.ImhNoEditMenu)
        self.lcdNumber.setSegmentStyle(QtWidgets.QLCDNumber.Filled)
        self.lcdNumber.setProperty("intValue", 0)
        self.lcdNumber.setObjectName("lcdNumber")
        """
        self.pushButton_2 = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_2.setGeometry(QtCore.QRect(250, 620, 106, 30))
        self.pushButton_2.setLocale(QtCore.QLocale(QtCore.QLocale.Erzya, QtCore.QLocale.Russia))
        self.pushButton_2.setAutoDefault(False)
        self.pushButton_2.setDefault(False)
        self.pushButton_2.setFlat(False)
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton_3 = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_3.setGeometry(QtCore.QRect(410, 620, 106, 30))
        self.pushButton_3.setObjectName("pushButton_3")
        
        self.label_6 = QtWidgets.QLabel(self.centralwidget)
        self.label_6.setGeometry(QtCore.QRect(750, 30, 661, 551))
        self.label_6.setText("")
        self.label_6.setPixmap(QtGui.QPixmap("../../../../../../Downloads/IRAS 16253 (1).png"))
        self.label_6.setScaledContents(False)
        self.label_6.setObjectName("label_6")
        
        self.pushButton_4 = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_4.setGeometry(QtCore.QRect(560, 620, 91, 31))
        self.pushButton_4.setObjectName("pushButton_4")
        """
        self.tableWidget = QtWidgets.QTableWidget(self.centralwidget)
        self.tableWidget.setGeometry(QtCore.QRect(377, 250, 275, 185))
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(2)
        self.tableWidget.setRowCount(3)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.horizontalHeader().setStretchLastSection(True)
        self.tableWidget.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        
        self.tableWidget.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(1, item)
        self.pushButton_5 = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_5.setGeometry(QtCore.QRect(280, 410, 371, 41))
        self.pushButton_5.setObjectName("pushButton_5")
        self.gridLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(750,50,661,661))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        #self.label_6.raise_()
   #     self.pushButton.raise_()
        self.textEdit.raise_()
#        self.textEdit_2.raise_()
 #       self.plainTextEdit.raise_()
        self.label.raise_()
        self.label_filename.raise_()
        self.label_2.raise_()
     #   self.label_3.raise_()
        #self.label_4.raise_()
        self.label_Temp.raise_()
        #self.plainTextEdit_2.raise_()
        #self.plainTextEdit_3.raise_()
        self.lcdNumber.raise_()
        #self.pushButton_2.raise_()
        #self.pushButton_3.raise_()
        #self.pushButton_4.raise_()
        self.tableWidget.raise_()
        self.pushButton_5.raise_()
        self.gridLayoutWidget.raise_()
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        
###################### 
             
        #self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        #self.gridLayout = self.gridLayout
        self.gridLayout.setObjectName("gridLayout")
        
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.spacerItem = QtWidgets.QSpacerItem(40, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(self.spacerItem)
        self.gridLayout.addLayout(self.horizontalLayout, 0, 0, 1, 1)
        
        
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.spacerItem1 = QtWidgets.QSpacerItem(40, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        self.verticalLayout.addItem(self.spacerItem1)
        self.gridLayout.addLayout(self.verticalLayout, 1, 0, 1, 1)
        #self.gridLayout.addLayout(self.horizontalHeaderItem, 0.5, 0.5, 0.25, 0.25)
        
        
        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
        
        self.canv = MatplotlibCanvas(self)
        #self.df = []
		
        self.toolbar = Navi(self.canv,self.centralwidget)
        self.horizontalLayout.addWidget(self.toolbar)
        
##########################
        
        #self.pushButton_2.clicked.connect(self.sum)
        #self.pushButton_3.clicked.connect(self.multiply)
        #self.pushButton_4.clicked.connect(self.plainTextEdit_2.clear)
        #self.pushButton_4.clicked.connect(self.plainTextEdit_3.clear)
        self.pushButton_5.clicked.connect(self.plot_table)
                
        


    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
#        self.pushButton.setText(_translate("MainWindow", "Sign Up"))
        self.label.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" color:#000000;\">Target name</span></p></body></html>"))
        self.label_filename.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" color:#000000;\">File name</span></p></body></html>"))
        self.label_2.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" color:#000000;\">Slit Configurations <br> (Enter in arcsec) </span></p></body></html>"))
  #      self.label_3.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" color:#000000;\">Confirm Password</span></p></body></html>"))
        #self.label_4.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" color:#000000;\">A</span></p><p><br/></p></body></html>"))
        self.label_Temp.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" color:black;\">Temperature</span></p><p><br/></p></body></html>"))
        #self.pushButton_2.setText(_translate("MainWindow", "SUM"))
        #self.pushButton_3.setText(_translate("MainWindow", "Multiply"))
        #self.pushButton_4.setText(_translate("MainWindow", "Clear"))
        item = self.tableWidget.verticalHeaderItem(0)
        item.setText(_translate("MainWindow", "0"))
        item = self.tableWidget.verticalHeaderItem(1)
        item.setText(_translate("MainWindow", "1"))
        item = self.tableWidget.verticalHeaderItem(2)
        item.setText(_translate("MainWindow", "2"))
        
        item = self.tableWidget.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "Slit Width"))
        item = self.tableWidget.horizontalHeaderItem(1)
        item.setText(_translate("MainWindow", "Offset"))
        self.pushButton_5.setText(_translate("MainWindow", "Load Image"))
    
    def sum(self):
    
        a = self.plainTextEdit_2.toPlainText()
        b = self.plainTextEdit_3.toPlainText()

        try:
            c = int(a)+int(b)
            self.lcdNumber.display(c)
        except:
            print("Invalid Inputs")

    def multiply(self):
        a = self.plainTextEdit_2.toPlainText()
        b = self.plainTextEdit_3.toPlainText()

        try:
            c = int(a)*int(b)
            self.lcdNumber.display(c)
        except:
            print("Invalid Inputs")
            
    def plot_table(self):
        try:
            source = self.textEdit.text()
            file_name = self.textEdit_filename.text()
            
            coord = SkyCoord.from_name(source, frame = 'icrs')
            print('Coordinates of the requested objects are: ',coord.ra, coord.dec)
    
            if file_name[-5:]==".fits":
                hdu = fits.open(file_name)[0]
                wcs = WCS(hdu.header)
                
            else:    
                paths = SkyView.get_images(position=source, height=10*u.arcmin,  width=10*u.arcmin,
                                        survey=['2MASS-K'])
                hdu = paths[0][0]
                wcs = WCS(hdu.header)
    
            size = u.Quantity((20, 20), u.arcmin)
            stamp = Cutout2D(hdu.data, coord, size, wcs=wcs)
            #fig, ax = plt.subplots( subplot_kw={'projection': stamp.wcs}, dpi = 300)
            projection = stamp.wcs
            """
            x = []
            y = []
            for i in range(4):    
                x.append(float(self.tableWidget.item(i, 0).text()))
                y.append(float(self.tableWidget.item(i, 1).text()))
            
            print(x,y)
            fg = plt.plot(x,y,'.')
            """
            try:
                self.horizontalLayout.removeWidget(self.toolbar)
                self.verticalLayout.removeWidget(self.canv)
                sip.delete(self.toolbar)
                sip.delete(self.canv)
                self.toolbar = None         
                self.canv = None
                self.verticalLayout.removeItem(self.spacerItem1)
                
            except Exception as e:
                print(e)
                pass
            #self.label_6.pixmap(QtGui.QPixmap(fg))
            self.canv = MatplotlibCanvas(self)
            self.toolbar = Navi(self.canv,self.centralwidget)
            self.horizontalLayout.addWidget(self.toolbar)
            self.verticalLayout.addWidget(self.canv)
            #self.canv.axes.cla()
            self.canv.fig
            
            
            ax = self.canv.fig.add_subplot(111, projection=stamp.wcs)
            ax.imshow(stamp.data, origin='lower', cmap="jet", norm=LogNorm()) 
            
            ## Slit configuration

            fov_width = 3.1*60 ## in arcsec
            delta_fov = 9.1*60 ## in arcsec
            coord_fov = SkyCoord(coord.ra+(0/3600)*u.deg, coord.dec+(0)*u.arcsec)                     
            sky_region = RectangleSkyRegion(coord_fov, width=fov_width *u.arcsec, height=9.1*u.arcmin)
            pixel_region = sky_region.to_pixel(stamp.wcs)
            artist = pixel_region.as_artist(color='gray', lw=2)            
            ax.add_artist(artist)
#            ax.text(0.1, 0.9, '', transform=ax.transAxes, c='lime', weight="bold")
            

            
            slit0_width = float(self.tableWidget.item(0, 0).text())
            delta_x0 = float(self.tableWidget.item(0, 1).text())
            coord_slit0 = SkyCoord(coord.ra+(delta_x0/3600)*u.deg, coord.dec+9.1/3*u.arcmin)                     
            sky_region = RectangleSkyRegion(coord_slit0, width=slit0_width *u.arcsec, height=9.1/3*u.arcmin)
            pixel_region = sky_region.to_pixel(stamp.wcs)
            artist = pixel_region.as_artist(color='lime', lw=2)            
            ax.add_artist(artist)
            ax.text(0.1, 0.9, 'Slit 0', transform=ax.transAxes, c='lime', weight="bold")
            
            slit1_width = float(self.tableWidget.item(1, 0).text())
            delta_x1 = float(self.tableWidget.item(1, 1).text())
            coord_slit1 = SkyCoord(coord.ra+(delta_x1/3600)*u.deg, coord.dec)                     
            sky_region = RectangleSkyRegion(coord_slit1, width=slit1_width *u.arcsec, height=9.1/3*u.arcmin)
            pixel_region = sky_region.to_pixel(stamp.wcs)
            artist = pixel_region.as_artist(color='red', lw=2)            
            ax.add_artist(artist)
            ax.text(0.1,0.85, 'Slit 1', transform=ax.transAxes, c='red', weight="bold")
            #ax.scatter(coord.ra.value, coord.dec.value, transform=ax.get_transform('world'))

            slit2_width = float(self.tableWidget.item(2, 0).text())
            delta_x2 = float(self.tableWidget.item(2, 1).text())
            coord_slit2 = SkyCoord(coord.ra+(delta_x2/3600)*u.deg, coord.dec-9.1/3*u.arcmin)                     
            sky_region = RectangleSkyRegion(coord_slit2, width=slit2_width *u.arcsec, height=9.1/3*u.arcmin)
            pixel_region = sky_region.to_pixel(stamp.wcs)
            artist = pixel_region.as_artist(color='yellow', lw=2)            
            ax.add_artist(artist)
            ax.text(0.1, 0.8, 'Slit 2', transform=ax.transAxes, c='yellow', weight="bold")
            
            ax.plot_coord(coord)
            ax.set_xlabel('RA')
            ax.set_ylabel('Dec')
            ax.set_title("2MASS Image (band Ks)")            
            """
            ax.plot(x,y,'o')
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_title("Matplotlib integration")
            """
            
            #plt.legend(loc='best')
            self.canv.draw()
        except:
            print("Object not found")
            from PyQt5.QtWidgets import QMessageBox

            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Error")
            msg.setInformativeText('Object could not be found or slit configuration invalid')
            msg.setWindowTitle("Error")
            msg.exec_()



if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
