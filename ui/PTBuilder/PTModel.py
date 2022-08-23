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

__version__ = '0.1'
__author__ = 'Christian Dewey'

'''
LCICPMS data GUI

2022-04-21
Christian Dewey
'''


class PTModel:
	''' model class for LCICPMS functions'''

	def __init__(self, ptview, mainview):
		"""Controller initializer."""
		self._ptview = ptview
		self._mainview = mainview
		self.ntime = True
		self.intColors = sns.color_palette(n_colors = 6, as_cmap = True)
		
	def importData(self):
		'''imports cal .csv file'''
		fdir = self._calview.calibrationDir + self._calview.listwidget.currentItem().text()
		self._data = pd.read_csv(fdir,sep=';',skiprows = 0, header = 1)

	def plotActiveMetals(self):
		'''plots active metals for selected file'''
		self._calview.chroma = plotChroma(self._calview, self._calview.metalOptions, self._data, self._calview.activeMetals)._plotChroma()

	def integrate(self, intRange):
		'''integrates over specified x range'''
		self.intRange = intRange
		pa_dict = {}
		for metal in self._calview.activeMetals:
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

			icpms_dataToSum = self._data[metal].iloc[i_tmin:i_tmax]
			#print(icpms_dataToSum)

			me_col_ind = self._data.columns.get_loc(metal)
			summed_area = 0
			for i in range(i_tmin, i_tmax):
				icp_1 = self._data.iloc[i,me_col_ind] # cps
				icp_2 = self._data.iloc[i+1,me_col_ind]
				min_height = min([icp_1,icp_2])
				max_height = max([icp_1,icp_2])
				timeDelta = self._data.iloc[i+1,me_col_ind - 1] - self._data.iloc[i,me_col_ind - 1] # seconds; time is always to left of metal signal
				#print(i, i+1, timeDelta)
				#print(min_height, max_height)
				rect_area = timeDelta * min_height
				top_area = timeDelta * (max_height - min_height) * 0.5
				An = rect_area + top_area
				summed_area = summed_area + An  # area =  cps * sec = counts
			print(metal + ': ' + str(summed_area/60))
	
			pa_dict[metal] = summed_area/60
				
			self._calview.n_area = pa_dict

		filename =  self._calview.calibrationDir + 'calibration_areas.txt' 
		with open(filename, 'a', newline = '') as csvfile:
			fwriter = csv.DictWriter(csvfile, fieldnames=pa_dict.keys())
			if self.ntime == True:
				fwriter.writeheader()
				self.ntime = False
			fwriter.writerow(pa_dict) 

	def plotLowRange(self,xmin,n):
		'''plots integration range'''
		col = self.intColors[n]
		minline = pg.InfiniteLine(xmin, pen = col, angle = 90)
		self._calview.plotSpace.addItem(minline) #InfiniteLine(minInt,angle = 90)
		
	def plotHighRange(self,xmax,n):
		col = self.intColors[n]
		maxline = pg.InfiniteLine(xmax, pen=col,angle = 90)
		self._calview.plotSpace.addItem(maxline)

	def calcLinearRegression(self):
		calCurve_dict = {}
		saveDict = {}
		metals = self._calview.activeMetals
		blank_value = 0
		blank_dict = {}

		for m in metals:
			pas = []
			concs = [] 
			for std in self._calview.standards.keys():
				std_list_n = self._calview.standards[std]
				if (std == 'Blank') and (len(std_list_n) > 0):
					blank_dict = std_list_n[0] 
					blank_value = blank_dict[m]
					print('blank PA for ' + m + ' = %.2f' % blank_value)
				if len(std_list_n) > 0:
					std_dict = std_list_n[0]
					pas.append(std_dict[m]-blank_value)
					concs.append(std_list_n[1])
			print(pas, concs)
			X = np.array(pas).reshape(-1, 1)
			y = np.array(concs)
			regr = linear_model.LinearRegression(fit_intercept=False)
			regr.fit(X, y)

			y_pred = regr.predict(X)

			# Print the Intercept:
			print("Metal: " + m)
			print('Intercept:', regr.intercept_)

			# Print the Slope:
			print('Slope:', regr.coef_[0]) 
			# The mean squared error
			mse = mean_squared_error(y, y_pred)
			print("Mean squared error: %.2f" % mse)
			# The coefficient of determination: 1 is perfect prediction
			r2 = r2_score(y, y_pred)
			print("Coefficient of determination: %.2f" % r2)
			
			fig, host = plt.subplots()
			host.scatter(X/1000, y, color="black")
			host.plot(X/1000, y_pred, color="blue", linewidth=3)

			host.set_xlabel(r'Peak Area ($10^3$ ICP-MS counts)')
			host.set_ylabel('Standard Conc. (ppb)')
			host.text(0.8,0.5,'$R^2$ = %.4f' % r2, transform = host.transAxes)
			host.text(0.8,0.4,'MSE = %.2f' % mse, transform = host.transAxes)
		
			host.set_title(m)

			fname = self._calview.calibrationDir + m + '_calibration.png'
			plt.savefig(fname, dpi = 300)
			plt.show()
		
			calCurve_dict[m] = [regr,(mean_squared_error(X, y),r2_score(X, y))]
			saveDict[m] = {'m': regr.coef_[0], 'b': regr.intercept_, 'r2': r2, 'mse': mse}
			
		self._mainview.calCurves = saveDict

		savefile = self._mainview.homeDir + 'calibration_curve.calib'
		with open(savefile, 'w') as file:
			file.write(json.dumps(saveDict))

		savefile = self._calview.calibrationDir + 'calibration_curve.calib'
		with open(savefile, 'w') as file:
			file.write(json.dumps(saveDict))
		
		

