import pandas as pd
#import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, MaxNLocator,PercentFormatter)
import pyqtgraph as pg

class plotChroma:
	def __init__(self,view = None,metalList=None,icpms_data=None,activeMetals=None,plt_title = None):
		self._view = view
		self.metalList= metalList
		self.activeMetals = metalList
		self.icpms_data = icpms_data
		self.max_time = float(self._view.maxx.text())# #max(self.icpms_data['Time ' + self.activeMetals[0]]) / 60
		self.min_time = float(self._view.minx.text())	
		self.max_icp = 0
		self.min_icp = 0
		self.plt_title = plt_title

	def _plotICPMSChroma(self):
		offset = self._view.offset_float
		colors = sns.color_palette(n_colors = len(self.metalList),as_cmap = True)
		c_keys = self.metalList
		color_dict = {c_keys[i]: colors[i] for i in range(len(c_keys))}
		self._view.ICPMS_PlotSpace.clear()
		
		
		for m in self.activeMetals:
			icpms_time = self.icpms_data['Time ' + m] / 60
			#print(icpms_time.head())
			icpms_signal = self.icpms_data[m] / 1000

			if self.max_icp <= max(icpms_signal):
				self.max_icp = max(icpms_signal)

			msize = len(m)
			mass = m[:msize-2]
			element = m[msize-2:]
			self._view.ICPMS_PlotSpace.setBackground('w')
			pen = color_dict[m]
			self._view.ICPMS_PlotSpace.addLegend(offset = [-1,1])
			#chromaPlot = pg.PlotItem(icpms_time, icpms_signal,pen=pen,width = 2,name = m)
			self._view.ICPMS_PlotSpace.setXRange(self.min_time,self.max_time,padding = None)
			self._view.ICPMS_PlotSpace.setYRange(0,self.max_icp*1.1,padding = None)

			chromaPlot = self._view.ICPMS_PlotSpace.plot(icpms_time + offset/60, icpms_signal,pen=pen,width = 2,name = m)
			#viewPlot = self._view.ICPMS_PlotSpace.addItem(chromaPlot)
			#self._view.ICPMS_PlotSpace.plot(icpms_time, icpms_signal,pen=pen,width = 2,name = m)
			#chromaplots.append(chromaPlot)
		return chromaPlot

	def _plotTIC(self,esi_parser):

		tic = esi_parser.get_tic(ms_type = 'MS')[0]

		self._view.ESIMS_PlotSpace.setBackground('w')
		
		pen = 'red'
		tic_np = np.array(tic.tic)
		time_np = np.array(tic.time)

		print('min time: ' + str(self.min_time))
		print('max time: ' + str(self.max_time))

		if (np.max(time_np) >= self.max_time) and (np.min(time_np) <= self.min_time):
			tic_sub = tic_np[np.where(time_np > self.min_time)] # and (time_np > self.max_time))]
			tic_sub = tic_sub[np.where(time_np < self.max_time)]
		
		elif (np.max(time_np) >= self.max_time) and (np.min(time_np) > self.min_time):
			tic_sub = tic_np[np.where(time_np < self.max_time)]

		elif (np.max(time_np) < self.max_time) and (np.min(time_np) <= self.min_time):
			tic_sub = tic_np[np.where(time_np > self.min_time)]

		else: 
			tic_sub = tic_np

		maxTic = np.max(tic_sub) * 1.1
		msPlot = self._view.EIC_PlotSpace.plot(tic.time, tic.tic,pen=pen,width = 2,name ='TIC')

		self._view.EIC_PlotSpace.setXRange(self.min_time,self.max_time,padding = None)

		self._view.EIC_PlotSpace.setYRange(0,maxTic,padding = None)
		#ax.axhline(y=upperlimit, c='r')
		#ax.axhline(y=lowerlimit, c='r')
		return msPlot

	def _plotEIC(self,esi_parser,mzs):

		tic = esi_parser.get_tic(ms_type = 'MS')[0]

		eics = esi_parser.get_eics(mzs, tic_data = tic, ms_type = 'MS')[0]

		self._view.EIC_PlotSpace.clear()

		self._view.EIC_PlotSpace.setBackground('w')

		colors = sns.color_palette(n_colors = len(mzs),as_cmap = True)
		c_keys = mzs
		color_dict = {c_keys[i]: colors[i] for i in range(len(c_keys))}
		
		maxEic = 0 

		for mz in mzs:
			eic = eics[mz]

			eic_np = np.array(eic.eic)
			time_np = np.array(eic.time)

			if (np.max(time_np) >= self.max_time) and (np.min(time_np) <= self.min_time):
				eic_sub = eic_np[np.where(time_np > self.min_time) and (np.where(time_np < self.max_time))] # and (time_np > self.max_time))]
			
			elif (np.max(time_np) >= self.max_time) and (np.min(time_np) > self.min_time):
				eic_sub = eic_np[np.where(time_np < self.max_time)]

			elif (np.max(time_np) < self.max_time) and (np.min(time_np) <= self.min_time):
				eic_sub = eic_np[np.where(time_np > self.min_time)]

			else: 
				eic_sub = eic_np

			if maxEic >= np.max(eic_sub):
				maxEic = maxEic
			elif maxEic < np.max(eic_sub):
				maxEic = 1.1 * np.max(eic_sub)

			pen = color_dict[mz]
			self._view.EIC_PlotSpace.addLegend(offset = [-1,1])

			mzl = f"{mz:.5f}"
			msPlot = self._view.EIC_PlotSpace.plot(np.array(eic.time), eic.eic,pen=pen,width = 2,name ='EIC '+ mzl) # offset is defined for icpms data; however, because icpms data can be plotted before offset is known, eic is adjusted 
			
			
		#ax.axhline(y=upperlimit, c='r')
		self._view.EIC_PlotSpace.setXRange(self.min_time,self.max_time,padding = None)
		self._view.EIC_PlotSpace.setYRange(0,maxEic,padding = None)
		#ax.axhline(y=lowerlimit, c='r')
		return msPlot


			


		
