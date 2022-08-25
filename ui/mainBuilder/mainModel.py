from curses import meta
import time
from datetime import datetime
from datetime import timedelta
import sys 
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import * 
from pyqtgraph import PlotWidget, plot
import pyqtgraph as pg
from functools import partial
import os
import pandas as pd
from functools import partial
import seaborn as sns
import csv

from .mainChroma import *

import corems
from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters

__version__ = '0.1'
__author__ = 'Christian Dewey'

'''
LC ICP/ESIMS data GUI

2022-08-24
Christian Dewey
'''


class MainModel:

	def __init__(self, mainview):
		"""Controller initializer."""
		self._mainview = mainview
		self.intColors = sns.color_palette(n_colors = 6, as_cmap = True)

		self.esi_parser = self._mainview.esi_parser
		

	def importData_generic(self,fdir):
		'''imports LCICPMS .csv file'''
		data = pd.read_csv(fdir,sep=';',skiprows = 0, header = 1)

		return data 

	def plotActiveMetalsMP(self):
		'''plots active metals for selected file'''
		activeMetalsPlot = ICPMS_Data_Class(self._data,self._view.activeMetals)
		activeMetalsPlot.chroma().show()
	
	def plotActiveMetals(self):
		'''plots active metals for selected file'''
		self._mainview.chroma = plotChroma(self._mainview, self._mainview.icpms_elements_for_plotting, self.icpms_data, self._mainview.activeMetals)._plotChroma()
		if self.minline != None:
			self._view.plotSpace.addItem(self.minline)
		if self.maxline != None:
			self._view.plotSpace.addItem(self.maxline)

	def integrate(self, intRange):
		'''integrates over specified x range'''
		self.intRange = intRange
		time_holders = {'start_time': 0, 'stop_time' : 0}
		metalList = ['55Mn','56Fe','59Co','60Ni','63Cu','66Zn','111Cd', '208Pb']
		metal_dict= {key: None for key in metalList}
		corr_dict = {'correction': None}
		tstamp = {'timestamp': None}
		metalConcs = {**tstamp,**time_holders,**corr_dict,**metal_dict}
		peakAreas = {**tstamp,**time_holders,**corr_dict,**metal_dict}

		print(self._view.normAvIndium)
		if self._view.normAvIndium > 0:
			indium_col_ind = self._data.columns.get_loc('115In')
			if len(self._data['Time 115In']) > 2000:
				corr_factor = np.average(self._data.iloc[550:2500,indium_col_ind]) / self._view.normAvIndium  #550:2500 indices correspond to ~ 150 to 350 sec
			else:
				corr_factor = np.average(self._data.iloc[:,indium_col_ind]) / self._view.normAvIndium  #550:2500 indices correspond to ~ 150 to 350 sec
			print('\ncorrection factor: %.4f' % corr_factor)
		else:
			corr_factor = 1

		for metal in self._view.activeMetals:
			if metal != '115In':
				time = self._data['Time ' + metal] / 60
				range_min = self.intRange[0]
				range_max = self.intRange[1]
				min_delta = min(abs(time - range_min))
				max_delta = min(abs(time - range_max))
				i_tmin = int(np.where(abs(time - range_min) == min_delta )[0][0])
				i_tmax = int(np.where(abs(time - range_max) == max_delta )[0][0])
				minval = self._data.iloc[i_tmin]
				minval = minval['Time ' + metal]
				#print( i_tmin, minval/60, range_min)

				maxval = self._data.iloc[i_tmax]
				maxval = maxval['Time ' + metal]
				#print( i_tmax, maxval/60, range_max)

				#print(icpms_dataToSum)
				metalConcs['start_time'] = '%.2f' % range_min
				metalConcs['stop_time'] = '%.2f' % range_max
				peakAreas['start_time'] = '%.2f' % range_min
				peakAreas['stop_time'] = '%.2f' % range_max

				metalConcs['correction'] = '%.3f' % corr_factor 
				peakAreas['correction'] = '%.3f' % corr_factor 

				dateTimeObj = datetime.now()
				timestampStr = dateTimeObj.strftime("%d-%b-%Y (%H:%M:%S)")
				metalConcs['timestamp'] = timestampStr
				peakAreas['timestamp'] = timestampStr

				me_col_ind = self._data.columns.get_loc(metal)
				summed_area = 0
				timeDelta = 0
				for i in range(i_tmin, i_tmax+1):
					icp_1 = self._data.iloc[i,me_col_ind] / corr_factor# cps
					icp_2 = self._data.iloc[i+1,me_col_ind] / corr_factor
					min_height = min([icp_1,icp_2])
					max_height = max([icp_1,icp_2])
					#print('min height: %.2f' % min_height) 
					#print('max height: %.2f' % max_height) 
					
					timeDelta = (self._data.iloc[i+1,me_col_ind - 1] - self._data.iloc[i,me_col_ind - 1])/60 # minutes; time is always to left of metal signal
					#print('time step: %.4f' % timeDelta) 
					#print(i, i+1, timeDelta)
					#print(min_height, max_height)
					rect_area = timeDelta * min_height
					top_area = timeDelta * (max_height - min_height) * 0.5
					An = rect_area + top_area
					#print('An: %.2f' % An )
					#print('rect area: %.2f' % rect_area)
					#print('top area: %.2f' % top_area)
					#print('dArea: %.2f' % An)
					summed_area = summed_area + An  # area =  cps * sec = counts
				
				#print('yes base subtract')
				if self._view.baseSubtract == True:
					#print('yes base subtract')
					baseline_height_1 = self._data.iloc[i_tmin,me_col_ind] / corr_factor
					baseline_height_2 =  self._data.iloc[i_tmax,me_col_ind] / corr_factor
					baseline_timeDelta = (self._data.iloc[i_tmax,me_col_ind - 1] - self._data.iloc[i_tmin,me_col_ind - 1])/60 #minutes
					#print('baseline_height_1: %.2f' % baseline_height_1)
					#print('baseline_height_2: %.2f' % baseline_height_2)
					#print('timeDelta: %.2f' % baseline_timeDelta)

					min_base_height = min([baseline_height_1, baseline_height_2])
					max_base_height = max([baseline_height_1, baseline_height_2])
					#print('min_base_height: %.2f' % min_base_height)
					#print('max_base_height: %.2f' % max_base_height)
					baseline_area_1 = min_base_height * baseline_timeDelta
					baseline_area_2 = (max_base_height - min_base_height) * baseline_timeDelta * 0.5
					#print('baseline_area_1: %.2f' % baseline_area_1)
					#print('baseline_area_2: %.2f' % baseline_area_2)
					
					#print('summed_area: %.2f' % summed_area)
					baseline_area = baseline_area_1 + baseline_area_2
					summed_area = summed_area - baseline_area
					summed_area = max(summed_area,0)
					#print('baseline_area: %.2f' % baseline_area)
					

				cal_curve = self._view.calCurves[metal]	
				slope = cal_curve['m']
				intercept = cal_curve['b']
				conc_ppb = slope * summed_area + intercept
				conc_uM = conc_ppb / self._view.masses[metal]
				
				peakAreas[metal] = '%.1f' % summed_area
				metalConcs[metal] = '%.3f' % conc_uM
				print('\n' + metal + ' uM: %.3f' % conc_uM)
				print(metal  + ' peak area: %.1f' % summed_area)

		
		if self._view.singleOutputFile == False:
			filename =  self._view.homeDir + 'peaks_uM_' + self.fdir.split('/')[-1].split(',')[0]

			if os.path.exists(filename):
				with open(filename, 'a', newline = '') as csvfile:
					fwriter = csv.DictWriter(csvfile, fieldnames=metalConcs.keys())
					fwriter.writerow(metalConcs) 		
			else:
				csv_cols = ['start_time', 'stop_time','correction'] + metalList
				with open(filename, 'w', newline = '') as csvfile:
					fwriter = csv.writer(csvfile, delimiter = ',', quotechar = '|')
					if self._view.normAvIndium > 0:
						fwriter.writerow(['115In correction applied: %.3f' % corr_factor,''])
					fwriter.writerow(['concentrations in uM',''])
					fwriter.writerow(['time in minutes',''])
					fwriter.writerow(csv_cols)
				with open(filename, 'a', newline = '') as csvfile:
					fwriter = csv.DictWriter(csvfile, fieldnames=metalConcs.keys())
					fwriter.writerow(metalConcs) 	
		else:
			filename =  self._view.homeDir + 'concentrations_uM_all.csv' 

			metalConcs = {**{'filename':self.fdir.split('/')[-1].split(',')[0]},**metalConcs}
			if os.path.exists(filename):
				with open(filename, 'a', newline = '') as csvfile:
					fwriter = csv.DictWriter(csvfile, fieldnames=metalConcs.keys())
					fwriter.writerow(metalConcs) 		
			else:
				csv_cols = ['filename','tstamp','start_time', 'stop_time','correction'] + metalList
				with open(filename, 'w', newline = '') as csvfile:
					fwriter = csv.writer(csvfile, delimiter = ',', quotechar = '|')
				#	if self._view.normAvIndium > 0:
			#		fwriter.writerow(['115In correction applied: %.3f' % corr_factor,''])
					fwriter.writerow(['concentrations in uM',''])
					fwriter.writerow(['time in minutes',''])
					fwriter.writerow(csv_cols)
				with open(filename, 'a', newline = '') as csvfile:
					fwriter = csv.DictWriter(csvfile, fieldnames=metalConcs.keys())
					fwriter.writerow(metalConcs) 

			filename =  self._view.homeDir + 'peakareas_counts_all.csv' 

			peakAreas = {**{'filename':self.fdir.split('/')[-1].split(',')[0]},**peakAreas}
			if os.path.exists(filename):
				with open(filename, 'a', newline = '') as csvfile:
					fwriter = csv.DictWriter(csvfile, fieldnames=peakAreas.keys())
					fwriter.writerow(peakAreas) 		
			else:
				csv_cols = ['filename','tstamp','start_time', 'stop_time', 'correction'] + metalList
				with open(filename, 'w', newline = '') as csvfile:
					fwriter = csv.writer(csvfile, delimiter = ',', quotechar = '|')
					if self._view.normAvIndium > 0:
						fwriter.writerow(['115In correction applied: %.3f' % corr_factor,''])
					fwriter.writerow(['time in minutes',''])
					fwriter.writerow(csv_cols)
				with open(filename, 'a', newline = '') as csvfile:
					fwriter = csv.DictWriter(csvfile, fieldnames=peakAreas.keys())
					fwriter.writerow(peakAreas) 
			#print('Intercept: %.4f' % intercept)
			#print('Slope: %.8f' % slope) 

	def plotLowRange(self,xmin,n):
		'''plots integration range'''
		col = self.intColors[0]
		self.minline = pg.InfiniteLine(xmin, pen = col, angle = 90)
		self._view.plotSpace.addItem(self.minline) #InfiniteLine(minInt,angle = 90)
		
	def plotHighRange(self,xmax,n):
		col = self.intColors[0]
		self.maxline = pg.InfiniteLine(xmax, pen=col,angle = 90)
		self._view.plotSpace.addItem(self.maxline)

	def removeIntRange(self):
		self._view.plotSpace.removeItem(self.maxline)
		self._view.plotSpace.removeItem(self.minline)


