# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 00:09:39 2022

@author: Himanshu Tyagi

WARNING: Do not edit this file unless you know what you are doing.
"""


import serial
import os
from astropy import log
from astropy.coordinates import Angle, SkyCoord
from astropy.visualization import AsymmetricPercentileInterval
from regions import RectangleSkyRegion
from astropy.wcs import WCS
from astropy.io import fits
import matplotlib.pyplot as plt
from astroquery.skyview import SkyView
import astropy.units as u
from astropy.nddata import Cutout2D
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QMessageBox
import sys
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout, QHeaderView, QTabWidget
from PyQt5.QtCore import pyqtSignal, pyqtSlot, QObject, QThread, QTimer
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as Navi
import sip
import pandas as pd
from scipy import ndimage, misc
from astropy.visualization.interval import (PercentileInterval,
                                            AsymmetricPercentileInterval,
                                            ManualInterval, MinMaxInterval)

from astropy.visualization.stretch import (LinearStretch, SqrtStretch,
                                           PowerStretch, LogStretch,
                                           AsinhStretch)

from astropy.visualization.mpl_normalize import ImageNormalize
import threading
import numpy as np
#%%
y = 70


class SerialReaderThread(QThread):
    received_data = pyqtSignal(str)
    def __init__(self, serial_port):
        super().__init__()
        self.serial_port = serial_port
    def run(self):
        while True:
            if self.serial_port.in_waiting > 0:
                received_data = self.serial_port.readline().decode().strip()
                self.received_data.emit(received_data)

class MatplotlibCanvas(FigureCanvasQTAgg):
	def __init__(self,parent=None, dpi = 120):
		self.fig = Figure(figsize=[10,10], dpi = dpi)
		super(MatplotlibCanvas,self).__init__(self.fig)


def simple_norm(data, stretch='linear', power=1.0, asinh_a=0.1, log_a=1000,
                min_cut=None, max_cut=None, min_percent=None, max_percent=None,
                percent=None, clip=True):

    if percent is not None:
        interval = PercentileInterval(percent)
    elif min_percent is not None or max_percent is not None:
        interval = AsymmetricPercentileInterval(min_percent or 0.,
                                                max_percent or 100.)
    elif min_cut is not None or max_cut is not None:
        interval = ManualInterval(min_cut, max_cut)
    else:
        interval = MinMaxInterval()

    if stretch == 'linear':
        stretch = LinearStretch()
    elif stretch == 'sqrt':
        stretch = SqrtStretch()
    elif stretch == 'power':
        stretch = PowerStretch(power)
    elif stretch == 'log':
        stretch = LogStretch(log_a)
    elif stretch == 'asinh':
        stretch = AsinhStretch(asinh_a)
    else:
        raise ValueError('Unknown stretch: {0}.'.format(stretch))

    vmin, vmax = interval.get_limits(data)

    return ImageNormalize(vmin=vmin, vmax=vmax, stretch=stretch, clip=clip)


def rotate(p, origin=(0, 0), degrees=0):
    angle = np.deg2rad(degrees)
    R = np.array([[np.cos(angle), -np.sin(angle)],
                  [np.sin(angle),  np.cos(angle)]])
    o = np.atleast_2d(origin)
    p = np.atleast_2d(p)
    return np.squeeze((R @ (p.T-o.T) + o.T).T)


class SkyViewWorker(QtCore.QThread):
    finished = QtCore.pyqtSignal(object, object)  # (hdu, wcs) or (None, None) on error

    def __init__(self, position, height, width, survey, parent=None):
        super().__init__(parent)
        self.position = position
        self.height = height
        self.width = width
        self.survey = survey

    def run(self):
        try:
            paths = SkyView.get_images(position=self.position, height=self.height, width=self.width, survey=self.survey)
            hdu = paths[0][0]
            wcs = WCS(hdu.header).celestial
            self.finished.emit(hdu, wcs)
        except Exception as e:
            print("SkyView download error:", e)
            self.finished.emit(None, None)


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1480, 1000)
        MainWindow.setFocusPolicy(QtCore.Qt.TabFocus)
        MainWindow.setStyleSheet("\n"
                                 "\n"
                                 "background-color: rgb(255,255,255);\n"
                                 "font: 14pt \"Times\";")
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")

        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(10, 10, 1480, 1000))
        self.tabWidget.setObjectName("tabWidget")
        self.tab_3 = QtWidgets.QWidget()
        self.tab_3.setObjectName("tab_3")

        self.label = QtWidgets.QLabel(self.tab_3)
        self.label.setGeometry(QtCore.QRect(120, y+70*0, 131, 31))
        self.label.setObjectName("label")
        

        self.textEdit = QtWidgets.QLineEdit(self.tab_3)
        self.textEdit.setGeometry(QtCore.QRect(280, y+ 0*70, 371, 41))
        self.textEdit.setObjectName("textEdit")
        self.textEdit.setToolTip("Enter simbad searchable name")
        
        
        self.coords_label = QtWidgets.QLabel(self.tab_3)
        self.coords_label.setGeometry(QtCore.QRect(120, y+70*1, 131, 31))
        self.coords_label.setObjectName("coords_label")
        
        
        self.coords = QtWidgets.QLineEdit(self.tab_3)
        self.coords.setGeometry(QtCore.QRect(280, y+70*1, 371, 41))
        self.coords.setObjectName("coords")
        self.coords.setToolTip("Enter Coordinates (Optional, if Object name is given.)")
        
        
        self.field_rotation_value = QtWidgets.QLineEdit(self.tab_3)
        self.field_rotation_value.setGeometry(QtCore.QRect(280, y+70*2, 371, 41))
        self.field_rotation_value.setObjectName("field_rotation")
        self.field_rotation_value.setToolTip("Enter field rotation angle in degrees")
        
        
        self.field_rotation = QtWidgets.QLabel(self.tab_3)
        self.field_rotation.setGeometry(QtCore.QRect(120, y+70*2, 131, 31))
        self.field_rotation.setObjectName("field_rotation_value")
        
        
        self.label_filename = QtWidgets.QLabel(self.tab_3)
        self.label_filename.setGeometry(QtCore.QRect(120, y+70*3, 131, 31))
        self.label_filename.setObjectName("label")
        
        self.textEdit_filename = QtWidgets.QLineEdit(self.tab_3)
        self.textEdit_filename.setGeometry(QtCore.QRect(280, y+70*3, 250, 41))
        self.textEdit_filename.setObjectName("textEdit_filename")
        self.textEdit_filename.setToolTip("Optional (overrides Target Name)")

        self.file_upload = QtWidgets.QPushButton(self.tab_3)
        self.file_upload.setGeometry(QtCore.QRect(531, y+70*3, 120, 41))
        self.file_upload.setIcon( QtGui.QIcon("upload2.png"))
        self.file_upload.setObjectName("file_upload")
        self.file_upload.setToolTip("Upload a .fits file")
        
        
        self.target_list_fname = QtWidgets.QLabel(self.tab_3)
        self.target_list_fname.setGeometry(QtCore.QRect(120, y+70*4, 131, 31))
        self.target_list_fname.setObjectName("target_list_fname")
        
        self.target_list = QtWidgets.QLineEdit(self.tab_3)
        self.target_list.setGeometry(QtCore.QRect(280, y+70*4, 250, 41))
        self.target_list.setObjectName("target_list")
        self.target_list.setToolTip("Displays Targets in FOV (optional)")
        
        self.load_targets = QtWidgets.QPushButton(self.tab_3)
        self.load_targets.setGeometry(QtCore.QRect(531, y+70*4, 120, 41))
        self.load_targets.setIcon( QtGui.QIcon("csv.png"))
        self.load_targets.setObjectName("load_targets")
        self.load_targets.setToolTip("Load target list")  
        

        self.label_2 = QtWidgets.QLabel(self.tab_3)
        self.label_2.setGeometry(QtCore.QRect(120, y+70*5, 180, 120))
        self.label_2.setObjectName("label_2")

        
        self.tableWidget = QtWidgets.QTableWidget(self.tab_3)
        self.tableWidget.setGeometry(QtCore.QRect(377, y+70*5, 275, 225))
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(2)
        self.tableWidget.setRowCount(5)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(4, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.horizontalHeader().setStretchLastSection(True)
        self.tableWidget.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.tableWidget.verticalHeader().setStretchLastSection(True)
        self.tableWidget.verticalHeader().setSectionResizeMode(QHeaderView.Stretch)
        
        self.tableWidget.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(1, item)

        self.pushButton_5 = QtWidgets.QPushButton(self.tab_3)
        self.pushButton_5.setGeometry(QtCore.QRect(377, y+580, 275, 41))
        self.pushButton_5.setObjectName("pushButton_5")
        self.pushButton_5.setToolTip("Plot slit configuration over sky image")
        
        
        self.label_Temp = QtWidgets.QLabel(self.tab_3)
        self.label_Temp.setGeometry(QtCore.QRect(120, y+70*10, 131, 31))
        self.label_Temp.setObjectName("Temperature")
        
        self.lcdNumber = QtWidgets.QLCDNumber(self.tab_3)
        self.lcdNumber.setGeometry(QtCore.QRect(300, y+70*10, 161, 61))
        self.lcdNumber.setFocusPolicy(QtCore.Qt.NoFocus)
        self.lcdNumber.setStyleSheet("color: rgb(255, 0, 0);")
        self.lcdNumber.setInputMethodHints(QtCore.Qt.ImhNoEditMenu)
        self.lcdNumber.setSegmentStyle(QtWidgets.QLCDNumber.Filled)
        self.lcdNumber.setProperty("intValue", 0)
        self.lcdNumber.setObjectName("lcdNumber")
        
        
        self.gridLayoutWidget = QtWidgets.QWidget(self.tab_3)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(750,y,700,700))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")

        
        MainWindow.setCentralWidget(self.tab_3)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

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

        self.canv = MatplotlibCanvas(self)
        self.toolbar = Navi(self.canv,self.tab_3)
        self.horizontalLayout.addWidget(self.toolbar)
        
        self.pushButton_5.clicked.connect(self.plot_table)
        
        self.file_upload.clicked.connect(self.open_file)    
        self.load_targets.clicked.connect(self.open_target)    

        self.tabWidget.addTab(self.tab_3, "")

        ########   TAB ENG  #########################################################
        self.tab_4 = QtWidgets.QWidget()
        
        self.eng_ser_portno_fname = QtWidgets.QLabel(self.tab_4)
        self.eng_ser_portno_fname.setGeometry(QtCore.QRect(120, y+70*2, 151, 31))
        self.eng_ser_portno_fname.setObjectName("Serial port ID")
        
        self.eng_ser_port_no = QtWidgets.QLineEdit(self.tab_4)
        self.eng_ser_port_no.setGeometry(QtCore.QRect(377, y+70*2, 250, 41))
        self.eng_ser_port_no.setToolTip("Serial port ID input")
        
        
        self.eng_ser_port_fname = QtWidgets.QLabel(self.tab_4)
        self.eng_ser_port_fname.setGeometry(QtCore.QRect(120, y+70*3, 151, 31))
        self.eng_ser_port_fname.setObjectName("Serial port input")
        
        self.eng_ser_port = QtWidgets.QLineEdit(self.tab_4)
        self.eng_ser_port.setGeometry(QtCore.QRect(377, y+70*3, 250, 41))
        self.eng_ser_port.setObjectName("target_list")
        self.eng_ser_port.setToolTip("Serial port input")
        
        self.eng_ser_port_out = QtWidgets.QLabel(self.tab_4)
        self.eng_ser_port_out.setGeometry(QtCore.QRect(120, y+70*10, 200, 41))
        self.eng_ser_port_out.setObjectName("Serial port output")

        self.eng_ser_port_out_edit = QtWidgets.QTextEdit(self.tab_4, readOnly=True)
        self.eng_ser_port_out_edit.setGeometry(QtCore.QRect(377, y+70*10, 281, 41))
        self.eng_ser_port_out.setObjectName("Serial port outputbox")

        self.eng_label_2 = QtWidgets.QLabel(self.tab_4)
        self.eng_label_2.setGeometry(QtCore.QRect(120, y+70*4, 180, 120))
        self.eng_label_2.setObjectName("label_2")

        
        self.eng_tableWidget = QtWidgets.QTableWidget(self.tab_4)
        self.eng_tableWidget.setGeometry(QtCore.QRect(377, y+70*4, 275, 225))
        self.eng_tableWidget.setObjectName("tableWidget")
        self.eng_tableWidget.setColumnCount(2)
        self.eng_tableWidget.setRowCount(5)


        item = QtWidgets.QTableWidgetItem()
        self.eng_tableWidget.setVerticalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.eng_tableWidget.setVerticalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.eng_tableWidget.setVerticalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.eng_tableWidget.setVerticalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        self.eng_tableWidget.setVerticalHeaderItem(4, item)
        item = QtWidgets.QTableWidgetItem()

        self.eng_tableWidget.horizontalHeader().setStretchLastSection(True)
        self.eng_tableWidget.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.eng_tableWidget.verticalHeader().setStretchLastSection(True)
        self.eng_tableWidget.verticalHeader().setSectionResizeMode(QHeaderView.Stretch)
        
        self.eng_tableWidget.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.eng_tableWidget.setHorizontalHeaderItem(1, item)



        self.eng_pushButton_5 = QtWidgets.QPushButton(self.tab_4)
        self.eng_pushButton_5.setGeometry(QtCore.QRect(377, y+70*8, 275, 41))
        self.eng_pushButton_5.setObjectName("pushButton_5")
        self.eng_pushButton_5.setToolTip("Plot slit configuration over sky image")

        
        self.eng_gridLayoutWidget = QtWidgets.QWidget(self.tab_4)
        self.eng_gridLayoutWidget.setGeometry(QtCore.QRect(750,y,661,661))
        self.eng_gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.eng_gridLayout = QtWidgets.QGridLayout(self.eng_gridLayoutWidget)
        self.eng_gridLayout.setContentsMargins(0, 0, 0, 0)
        self.eng_gridLayout.setObjectName("gridLayout")

        
        MainWindow.setCentralWidget(self.tab_4)
        self.eng_statusbar = QtWidgets.QStatusBar(MainWindow)
        self.eng_statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.eng_statusbar)

        self.eng_gridLayout.setObjectName("gridLayout")        
        self.eng_horizontalLayout = QtWidgets.QHBoxLayout()
        self.eng_horizontalLayout.setObjectName("horizontalLayout")
        self.eng_spacerItem = QtWidgets.QSpacerItem(40, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        self.eng_horizontalLayout.addItem(self.eng_spacerItem)
        self.eng_gridLayout.addLayout(self.eng_horizontalLayout, 0, 0, 1, 1)
        
        
        self.eng_verticalLayout = QtWidgets.QVBoxLayout()
        self.eng_verticalLayout.setObjectName("verticalLayout")
        self.eng_spacerItem1 = QtWidgets.QSpacerItem(40, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        self.eng_verticalLayout.addItem(self.eng_spacerItem1)
        self.eng_gridLayout.addLayout(self.eng_verticalLayout, 1, 0, 1, 1)
        
        self.eng_canv = MatplotlibCanvas(self)
        self.eng_toolbar = Navi(self.eng_canv,self.tab_4)
        self.eng_horizontalLayout.addWidget(self.eng_toolbar)
        
        self.eng_pushButton_5.clicked.connect(self.eng_plot_table)
        self.eng_pushButton_5.clicked.connect(self.send_data)
        
        self.tabWidget.addTab(self.tab_4, "")
        #######################################################################
        
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        
        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
        
    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.label.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" color:#000000;\">Target name</span></p></body></html>"))
        self.coords_label.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" color:#000000;\">Coordinates</span></p></body></html>"))
        self.label_filename.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" color:#000000;\">File name</span></p></body></html>"))
        self.label_2.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" color:#000000;\">Slit Configurations <br> (Enter in arcsec) </span></p></body></html>"))
        self.target_list_fname.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" color:#000000;\">Target list</span></p></body></html>"))
        
        self.label_Temp.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" color:black;\">Temperature</span></p><p><br/></p></body></html>"))
        
        self.field_rotation.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" color:#000000;\">Field Rotation</span></p></body></html>"))

        item = self.tableWidget.verticalHeaderItem(0)
        item.setText(_translate("MainWindow", "0"))
        item = self.tableWidget.verticalHeaderItem(1)
        item.setText(_translate("MainWindow", "1"))
        item = self.tableWidget.verticalHeaderItem(2)
        item.setText(_translate("MainWindow", "2"))

        item = self.tableWidget.verticalHeaderItem(3)
        item.setText(_translate("MainWindow", "3"))
        item = self.tableWidget.verticalHeaderItem(4)
        item.setText(_translate("MainWindow", "4"))
        
        item = self.tableWidget.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "Slit Width"))
        item = self.tableWidget.horizontalHeaderItem(1)
        item.setText(_translate("MainWindow", "Offset"))
        self.pushButton_5.setText(_translate("MainWindow", "Load FOV"))

        
        #######################################################################

        self.eng_label_2.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" color:#000000;\">Slit Configurations <br> (Enter in mm) </span></p></body></html>"))
        self.eng_ser_port_fname.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" color:#000000;\">Serial Port Input</span></p></body></html>"))
        self.eng_ser_port_out.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" color:black;\">Serial Port Output</span></p><p><br/></p></body></html>"))
        self.eng_ser_port_out_edit.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" color:red;\">Serial Port Output</span></p><p><br/></p></body></html>"))
        self.eng_ser_portno_fname.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" color:#000000;\">Serial Port ID</span></p></body></html>"))


        item = self.eng_tableWidget.verticalHeaderItem(0)
        item.setText(_translate("MainWindow", "1"))
        item = self.eng_tableWidget.verticalHeaderItem(1)
        item.setText(_translate("MainWindow", "2"))
        item = self.eng_tableWidget.verticalHeaderItem(2)
        item.setText(_translate("MainWindow", "3"))

        item = self.eng_tableWidget.verticalHeaderItem(3)
        item.setText(_translate("MainWindow", "4"))
        item = self.eng_tableWidget.verticalHeaderItem(4)
        item.setText(_translate("MainWindow", "5"))
        
        item = self.eng_tableWidget.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "L Position"))
        item = self.eng_tableWidget.horizontalHeaderItem(1)
        item.setText(_translate("MainWindow", "Slit Width"))
        self.eng_pushButton_5.setText(_translate("MainWindow", "Configure Slit"))
               
        #########################################
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3), _translate("MainWindow", "Astronomers"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_4), _translate("MainWindow", "Engineers"))
        
        



class MainWindow(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.threadpool = QtCore.QThreadPool()
        #self.worker_thread = WorkerThread()
        self.setupUi(self)

        try:   
            pn = os.popen("ls /dev/ttyUSB*").read()[:-1]
            self.eng_ser_port_no.setText(pn)
            port_no = self.eng_ser_port_no.text()
            self.serial_port = serial.Serial(port_no, 9600)
            self.serial_thread = SerialReaderThread(self.serial_port)
            self.serial_thread.received_data.connect(self.receive_data)
            self.serial_thread.start()

        except:
             msg = QMessageBox()
             msg.setIcon(QMessageBox.Critical)
             msg.setText("Serial Port error")
             msg.setInformativeText(f"<html><head/><body><p><span style=\" color:black;\"> There is some issue in serial port communication.<br> Please try giving it permissions.<br> Use `sudo chmod o+rw /dev/ttyUSB*`</span></p><p><br/></p></body></html>")
             msg.setWindowTitle("Error")
             msg.exec_()


    def open_file(self):
            self.filename = QtWidgets.QFileDialog.getOpenFileName(None, 'Open File', '*.fits')
            print(self.filename[0])
            self.textEdit_filename.setText(self.filename[0])
        
    def open_target(self):
            self.filename = QtWidgets.QFileDialog.getOpenFileName(None, 'Open File', '*.csv')
            print(self.filename[0])
            self.target_list.setText(self.filename[0])


    def eng_open_file(self):
            self.filename = QtWidgets.QFileDialog.getOpenFileName(None, 'Open File', '*.fits')
            print(self.filename[0])
            self.eng_textEdit_filename.setText(self.filename[0])
        
    def eng_open_target(self):
            self.filename = QtWidgets.QFileDialog.getOpenFileName(None, 'Open File', '*.csv')
            print(self.filename[0])
            self.eng_target_list.setText(self.filename[0])

    def plot_table(self, vmin=None, vmid=None, vmax=None, pmin=0.25, pmax=99.75, 
                   stretch='linear', exponent=2):
        try:
            source_check = self.textEdit.text()
            coord_check = self.coords.text()
            file_name = self.textEdit_filename.text()
            
            if len(coord_check)!=0:
                coord_string = self.coords.text()
                try:
                    coord = SkyCoord(coord_string, frame="icrs")
                    
                except:
                    coord = SkyCoord(coord_string, frame="icrs", unit=(u.hourangle, u.deg))
                                    
            else:
                source = self.textEdit.text()                           
                coord = SkyCoord.from_name(source, frame = 'icrs')
                print('Coordinates of the requested objects are: ',coord.ra, coord.dec)
                #self.coords.setText(f"{str(coord.ra), str(coord.dec)}")
    
    
            if len(file_name) != 0:
                hdu = fits.open(file_name)
                
                if len(hdu)==0:
                        hdu = fits.open(file_name)[0]
                else:
                        hdu = fits.open(file_name)[1]
                
                wcs = WCS(hdu.header)
                
            elif len(file_name) !=0 and len(source_check)==0:
                # hdu = fits.open(file_name)[0]

                hdu = fits.open(file_name)

                if len(hdu)==0:
                        hdu = fits.open(file_name)[0]
                else:
                        hdu = fits.open(file_name)[1]
                        
                wcs = WCS(hdu.header).celestial
                coord_string = self.coords.text()
                coord = SkyCoord(coord_string, frame="icrs")
                
            if len(file_name) == 0 and (len(source_check) != 0 or len(coord_check) != 0):
                # Need to download from SkyView asynchronously
                if len(coord_check) != 0:
                    position = coord
                else:
                    position = source
                self.pushButton_5.setEnabled(False)
                self.statusbar.showMessage("Downloading image from SkyView...")
                self.skyview_worker = SkyViewWorker(position, 14*u.arcmin, 14*u.arcmin, ['2MASS-K'])
                self.skyview_worker.finished.connect(self.on_skyview_downloaded)
                self.skyview_worker.start()
                return

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

            self.canv = MatplotlibCanvas(self)
            self.toolbar = Navi(self.canv,self.tab_3)
            self.horizontalLayout.addWidget(self.toolbar)
            self.verticalLayout.addWidget(self.canv)
            self.canv.fig
            ax = self.canv.fig.add_subplot(111,)
            
            min_auto = vmin is None
            max_auto = vmax is None

            if min_auto or max_auto:
                interval = AsymmetricPercentileInterval(pmin, pmax,
                                                        n_samples=10000)

                try:
                    vmin_auto, vmax_auto = interval.get_limits(hdu.data)
                except (IndexError, TypeError):  # no valid values
                    vmin_auto = vmax_auto = 0

                vmin = vmin_auto
                vmax = vmax_auto

            # Prepare normalizer object
            if stretch == 'arcsinh':
                stretch = 'asinh'

            if stretch == 'log':
                if vmid is None:
                    if vmin < 0:
                        raise ValueError("When using a log stretch, if vmin < 0, then vmid has to be specified")
                    else:
                        vmid = 0.
                if vmin < vmid:
                    raise ValueError("When using a log stretch, vmin should be larger than vmid")
                log_a = (vmax - vmid) / (vmin - vmid)
                norm_kwargs = {'log_a': log_a}
            elif stretch == 'asinh':
                if vmid is None:
                    vmid = vmin - (vmax - vmin) / 30.
                asinh_a = (vmid - vmin) / (vmax - vmin)
                norm_kwargs = {'asinh_a': asinh_a}
            else:
                norm_kwargs = {}

            normalizer = simple_norm(hdu.data, stretch=stretch, power=exponent,
                                     min_cut=vmin, max_cut=vmax, clip=False,
                                     **norm_kwargs)

            # Adjust vmin/vmax if auto
            if stretch == 'linear':
                vmin = -0.1 * (vmax - vmin) + vmin
            # log.info("Auto-setting vmin to %10.3e" % vmin)
            if stretch == 'linear':
                vmax = 0.1 * (vmax - vmin) + vmax
            # log.info("Auto-setting vmax to %10.3e" % vmax)
            # Update normalizer object
            normalizer.vmin = vmin
            normalizer.vmax = vmax


            targetlist = self.target_list.text()

            if len(targetlist) != 0:
                df = pd.read_csv(targetlist)
                x, y = df['RA'], df['Dec']
                target_coords = SkyCoord(x, y, unit='deg')
                x, y = wcs.world_to_pixel(target_coords)

            else:
                pass


            if len(self.field_rotation_value.text()) != 0 and len(targetlist) != 0:
                image_rotation = float(self.field_rotation_value.text())    #in degrees
                rotated_array, (x1, y1) = ndimage.rotate(hdu.data, np.array([x, y].T), image_rotation, reshape=False)

            elif len(self.field_rotation_value.text()) != 0 and len(targetlist) == 0:
                image_rotation = float(self.field_rotation_value.text())    #in degrees
                rotated_array  = ndimage.rotate(hdu.data, image_rotation, reshape=False)
                
            else:
                image_rotation = 0
                rotated_array = ndimage.rotate(hdu.data, image_rotation, reshape=False)

            fov_center_x = hdu.header['CRPIX1']
            fov_center_y = hdu.header['CRPIX2']

            ax.imshow(rotated_array, origin='lower',
                      cmap="gist_earth", norm=normalizer)

            if len(targetlist) != 0:
                ax.scatter(x1, y1, c='r')

            else:
                pass
            
            # Slit configuration
            fov_width = 3.1*60  # in arcsec
            delta_fov = 9.1  # in arcmin
            coord_fov = SkyCoord(coord.ra+(0/3600)*u.deg,
                                 coord.dec+(0)*u.arcsec)

            reg_detector = RectangleSkyRegion(coord_fov,
                                            width=delta_fov*u.arcmin,
                                            height=delta_fov*u.arcmin)
            pixel_region = reg_detector.to_pixel(wcs)
            artist = pixel_region.as_artist(color='gray', lw=1)
            ax.add_artist(artist)     
                                            
            sky_region = RectangleSkyRegion(coord_fov,
                                            width=fov_width*u.arcsec,
                                            height=delta_fov*u.arcmin)
            pixel_region = sky_region.to_pixel(wcs)
            artist = pixel_region.as_artist(color='white', lw=1)
            ax.add_artist(artist)

            c = ['r', 'tab:orange', 'yellow', 'lime', 'pink']
            slit_fov_check = []
            
            for i in range(5):
                slit0_width = float(self.tableWidget.item(0+i, 0).text())
                delta_x0 = float(self.tableWidget.item(0+i, 1).text())
                coord_slit0 = SkyCoord(coord.ra+(delta_x0/3600)*u.deg,
                                       coord.dec+9.1*(2-i)/5*u.arcmin)
                sky_region = RectangleSkyRegion(coord_slit0,
                                                width=slit0_width*u.arcsec,
                                                height=9.1/5*u.arcmin)
                pixel_region = sky_region.to_pixel(wcs)
                artist = pixel_region.as_artist(color=c[i], lw=1)
                ax.add_artist(artist)
                ax.text(0.05, 0.9-i*0.05, f'Slit {i}', transform=ax.transAxes,
                        c=c[i],  weight="bold")
                if abs(delta_x0) > fov_width/2:
                    slit_fov_check.append(delta_x0)
            
            if len(slit_fov_check) > 0:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText("Slit Out of FOV")
                msg.setInformativeText('Atleast one slit is out of FOV')
                msg.setWindowTitle("Error")
                msg.exec_()

            #ax.plot_coord(coord)
            ax.set_xlabel('Pixel')
            ax.set_ylabel('Pixel')
            ax.set_title("2MASS Image (band Ks)")
            #ax.scatter(coord.ra.value, coord.dec.value,
            #           transform=ax.get_transform('world'), marker='o',
            #           c='None', edgecolors='white', s=100)
            self.canv.draw()
            
        except:
            print("Object not found")
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Error")
            msg.setInformativeText('Object could not be found or slit configuration invalid')
            msg.setWindowTitle("Error")
            msg.exec_()
    
    def on_skyview_downloaded(self, hdu, wcs):
        self.pushButton_5.setEnabled(True)
        self.statusbar.clearMessage()
        if hdu is None or wcs is None:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Download Error")
            msg.setInformativeText('Failed to download image from SkyView')
            msg.setWindowTitle("Error")
            msg.exec_()
            return 
        
    
    def eng_plot_table(self, vmin=None, vmid=None, vmax=None, pmin=0.25, pmax=99.75, 
                   stretch='linear', exponent=2):

        try:    
            try:
                self.eng_horizontalLayout.removeWidget(self.eng_toolbar)
                self.eng_verticalLayout.removeWidget(self.eng_canv)
                sip.delete(self.eng_toolbar)
                sip.delete(self.eng_canv)
                self.eng_toolbar = None         
                self.eng_canv = None
                self.eng_verticalLayout.removeItem(self.eng_spacerItem1)
                
            except Exception as e:
                print(e)
                pass
            
            self.eng_canv = MatplotlibCanvas(self)
            self.eng_toolbar = Navi(self.eng_canv,self.tab_4)
            self.eng_horizontalLayout.addWidget(self.eng_toolbar)
            self.eng_verticalLayout.addWidget(self.eng_canv)
            self.eng_canv.fig

            eng_ax = self.eng_canv.fig.add_subplot(111)
 
            fov_width = 8572 # in mm*100
 
            c = ['r', 'tab:orange', 'skyblue', 'lime', 'pink']
            slit_fov_check = []
            
            for i in range(5):
                l_pos = float(self.eng_tableWidget.item(0+i, 0).text())*100
                slit0_width = l_pos
                delta_x0 = float(self.eng_tableWidget.item(0+i, 1).text())*100

                if delta_x0==0:
                    eng_ax.barh(i, slit0_width, height=0.98, left=0, align='edge', color='gray', alpha = 0.8)  # left slit
                    eng_ax.barh(i, (fov_width-(slit0_width+delta_x0)), height=0.98, left=fov_width-(fov_width-(slit0_width+delta_x0)), align='edge', color='gray',  alpha = 0.8)
                else:
                    eng_ax.barh(i, slit0_width, height=0.98, left=0, align='edge', color=c[i], alpha = 0.8)  # left slit
                    eng_ax.barh(i, (fov_width-(slit0_width+delta_x0)), height=0.98, left=fov_width-(fov_width-(slit0_width+delta_x0)), align='edge', color=c[i],  alpha = 0.8)

                eng_ax.set_xlim(0, fov_width)
                eng_ax.set_ylim(0, 5)
                eng_ax.text(1.02, 0.1+i*0.2, f'Slit {i+1}', transform=eng_ax.transAxes,
                        c=c[i],  weight="bold")
                if slit0_width+(fov_width-(slit0_width+delta_x0)) > fov_width or l_pos>fov_width or l_pos<0 or l_pos+delta_x0>fov_width:
                    slit_fov_check.append(delta_x0)
            
            if len(slit_fov_check) > 0:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText("Slit Configuration Error")
                msg.setInformativeText('Slit Out of FOV or Slit Collision Scenario')
                msg.setWindowTitle("Error")
                msg.exec_()

            eng_ax.set_xlabel('mm$\\times 100$')
            eng_ax.set_ylabel('Slit number')
            eng_ax.set_title("Slit Configuration")
            self.eng_canv.draw()

        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Slit Configuration Error")
            msg.setInformativeText('Invalid inputs')
            msg.setWindowTitle("Error")
            msg.exec_()
            
    def send_data(self):
        data = self.eng_ser_port.text().encode()
        if len(data)!=0:
            self.serial_port.write(data)
        else:
            self.eng_ser_port_out_edit.setText(f"<html><head/><body><p><span style=\" color:red;\">Serial Port connection failure</span></p><p><br/></p></body></html>")

    def receive_data(self, received_data):
        self.eng_ser_port_out_edit.setText(f"<html><head/><body><p><span style=\" color:red;\">{received_data}</span></p><p><br/></p></body></html>")
        print(received_data)

    
class WorkerThread(QtCore.QThread):

    job_done = QtCore.pyqtSignal('QString')

    def __init__(self, parent=None):
        super(WorkerThread, self).__init__(parent)
        self.gui_text = None

    def do_work(self):
        source_check = "SSTc2d J034432.0+321143"
        coord_check = ""
        if len(coord_check) != 0:        
            paths = SkyView.get_images(position=coord_check, height=10*u.arcmin,  width=10*u.arcmin,
                                    survey=['2MASS-K'])
        else:
            paths = SkyView.get_images(position=source_check, height=10*u.arcmin,  width=10*u.arcmin,
                                    survey=['2MASS-K'])
        self.job_done.emit(paths[0])

    def run(self):
        self.do_work()

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(app.exec_())