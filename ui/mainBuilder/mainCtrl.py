import sys 
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtWidgets import * 
import pyqtgraph as pg
from functools import partial
import os
import pandas as pd
from functools import partial
import json

from ui.PTBuilder.PTView import *
from ui.PTBuilder.PTCtrl import *
from ui.PTBuilder.PTModel import *

from ui.mainBuilder.mainChroma import *

from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters

from corems.lc_icpms_ftms.factory.Aligned_ICP_ESI import Aligned_ICP_ESI as alignedMS

__version__ = '0.1'
__author__ = 'Christian Dewey'

'''

LC - ESI/ICP - MS data processing 

2022-08-22
Christian Dewey

'''

# Create a Controller class to connect the GUI and the model
class MainCtrl:
	"""Main Control class """

	window_closed = pyqtSignal()

	def __init__(self, mainmodel, mainview, ptview):
		"""Controller initializer."""
		self._mainmodel = mainmodel
		self._mainview = mainview
		self._ptview = ptview
		
		self._peakRange = []	

		self._range = []	

		self.n_clicks = 0
		self._n = 0
		self._xMin = 0
		self._xMax = 0
		self.button_is_checked = False

		# Connect signals and slots
		self._connectSignals()

	def _buildExpression(self): #, sub_exp):
		"""Build expression."""
		expression = None #sub_exp
		#self._mainview.setDisplayText(expression)

	def _showPeriodicTable(self):
		''' opens periodic table '''
		self._ptview = PTView(mainview=self._mainview)
		ptmodel = PTModel(ptview=self._ptview, mainview = self._mainview)
		PTCtrl(model=ptmodel, mainview=self._mainview, ptview= self._ptview,mainctrl=self)
		self._ptview.show()

	def _selectESIMSFile(self):
		''' opens window to select normalization file for 115In correction; saves average 115In signal from norm file'''
		if self.n_clicks < 2:
			dialog = QFileDialog()
			dialog.setWindowTitle("Select ESI-MS File")
			dialog.setViewMode(QFileDialog.Detail)
			self._mainview.esims_filepath = dialog.getOpenFileName(self._mainview,"Openfile")[0]
			MSParameters.mass_spectrum.threshold_method = 'signal_noise'
			MSParameters.mass_spectrum.s2n_threshold = 6
			MSParameters.ms_peak.peak_min_prominence_percent = 0.1
			self._mainview.esi_parser = rawFileReader.ImportMassSpectraThermoMSFileReader(self._mainview.esims_filepath)
			self._mainview._makeTICPlot()
			print(self._mainview.esims_filepath)
			self.n_clicks = self.n_clicks + 1 
			
			self._mainview.esi_file_selected = True

		self._mainview.replotBtn.setEnabled(True)
		

		#if self._mainview.icp_file_selected is True:
		#	self._mainview.buttons['Assign'].setEnabled(True)
					
	def _selectICPMSFile(self):
		''' opens window to select normalization file for 115In correction; saves average 115In signal from norm file'''
		if self.n_clicks < 2:
			dialog = QFileDialog()
			dialog.setWindowTitle("Select ICP-MS File")
			dialog.setViewMode(QFileDialog.Detail)
			self._mainview.icpms_filepath = dialog.getOpenFileName(self._mainview,"Openfile")[0]
			print(self._mainview.icpms_filepath)

			'''imports LCICPMS .csv file'''
			
			self._mainview.icpms_data = pd.read_csv(self._mainview.icpms_filepath,sep=';|,|\t',skiprows = 0, header = 0)
			self._mainview.icpms_data.dropna(inplace=True)	
			print(self._mainview.icpms_data.head())
			self._mainview.icpms_data.set_index(self._mainview.icpms_data['Number'], inplace=True)
			print(self._mainview.icpms_data.head())
			icpms_header = list(self._mainview.icpms_data.columns.values)
			print(icpms_header)
			icpms_elements = []
			for v in icpms_header:
				if ('Time ' in v) or ('Number' in v) or (' ' in v):
					continue
				else:
					icpms_elements.append(v)
			
			self._mainview.all_icpms_elements = icpms_elements

			self._mainview._createICPElementCheckBoxes()
			self.n_clicks = self.n_clicks + 1 
			
			for cbox in self._mainview.icp_checkBoxes.values():
				cbox.stateChanged.connect(partial(self._mainview.ICPMSClickBox, cbox))
			
			self._mainview.icp_file_selected = True

			#self._mainview._createPlotBox()
			self._mainview.replotBtn.setEnabled(True)
			self._mainview.buttons['Select peak'].setEnabled(True)

			#if self._mainview.esi_file_selected is True:
			#	self._mainview.buttons['Assign'].setEnabled(True)
				
	def _assignFormulas(self):

		esif = self._mainview.esims_filepath
		icpf = self._mainview.icpms_filepath 

		icpesi_obj = alignedMS(icpf, esif,self._mainview.heteroatom)

		icpesi_obj.offset = self._mainview.offset_float
	#	icpesi_obj.threshold = self._mainview.threshold_float
		icpesi_obj.timestart =  self._range[0]  #9.7413  #
		icpesi_obj.timestop = self._range[1]  # #10.57657 #

		#self._mainview.heteroatom = '63Cu'

		print('\n...assiging with ' + self._mainview.heteroatom + ' as heteroatom')

		icpesi_obj.getMSData(self._mainview.icpms_data,self._mainview.esi_parser)

	#	icpesi_obj.timestart = 8.2
#		icpesi_obj.timestop = 9.2

		#elementDict = {'C':(1,50), 'H':(4,100), 'O':(1,20), 'N':(0,4), 'S':(0,0), 'Cl':(0,0), 'Br':(0,0), 'P':(0,0), 'Na':(0,0), 'Cu':(0,1) }

		self._getThreshold()
		self._getOffset()

		elementDict = {}

		for element in self._mainview._currentElements:
			elementDict[element] = (int(self._mainview._elementMins[element].text()), int(self._mainview._elementMaxs[element].text()))

		self.best_results = icpesi_obj.assignFormulas(elementDict,0.4)

		print(self.best_results)

		self._mainview.best_results = self.best_results
		self._mainview.mzs = list(self.best_results['m/z'])
		
		self._mainview.EIC_PlotSpace.clear()
		#self.no_dups_mzs = pd.DataFrame(self._mainview.mzs).drop_duplicates()
		#self._mainview._makeEICPlot(self.no_dups_mzs.values.tolist())
		self._mainview._makeEICPlot(self._mainview.mzs)

		self._mainview.generateResultsReport()

	def _replotAll(self):

		print('Current offset: ' + str(self._mainview.offset_float) + ' s') 

		#try:
		self._mainview._makeICPMSPlot()
		#except:
		#	print('No ICPMS data to plot')
		
		if self._mainview.eic_plotted is True:
		
			#self._mainview._makeEICPlot(self.no_dups_mzs.values.tolist())
			self._mainview._makeEICPlot(self._mainview.mzs)

		else:
			try: 
				self._mainview._makeTICPlot()
			except: 
				print('No ESIMS data to plot')

		
	def _getOffset(self):
		
		offset = self._mainview.offset.text()
		
		self._mainview.offset_float = float(offset)

	def _getThreshold(self):
		
		threshold = self._mainview.threshold.text()
		
		self._mainview.threshold_float = float(threshold)

	
	def _onClick(self, event):

		#self._act_pos = self._mainview.ICPMS_PlotSpace.mapFromScene(event[0].scenePos())

		self._act_pos = self._mainview.ICPMS_PlotSpace.getPlotItem().vb.mapSceneToView(event[0].scenePos())
		cc = len(self._range)
		cc = cc + 1 

		if cc == 1:
			xt = self._act_pos.x()
			self._range.append(xt)
			self._plotRange(xt)
			print('selection: ' + str(cc) + ' ' +  str(xt))
		
		elif ( cc == 2 ):
			xt = self._act_pos.x()
			self._range.append(xt)
			self._plotRange(xt)
			print('selection: ' + str(cc) + ' ' +  str(xt))

			#self._mainview.buttons['Select peak'].setStyleSheet('background-color : gray')

		

	def _plotRange(self, x):
		col = 'gray'
		vline = pg.InfiniteLine(x, pen = col, angle = 90)
		self._mainview.ICPMS_PlotSpace.addItem(vline)


	def _getRange(self):

		#self._mainview.buttons['Select peak'].setStyleSheet('background-color : yellow')

		self._mainview.proxy = pg.SignalProxy(self._mainview.ICPMS_PlotSpace.scene().sigMouseClicked, rateLimit=60, slot=self._onClick)

		self._range.sort()

		self._mainview.buttons['Assign'].setEnabled(True)

		


	def _connectSignals(self):
		"""Connect signals and slots."""
		self._mainview.buttons['Select Elements'].clicked.connect(partial(self._buildExpression, ''))


		self._mainview.buttons['Select Elements'].setEnabled(True)
		self._mainview.buttons['Load ESI-MS Data'].setEnabled(True)
		self._mainview.buttons['Load ICP-MS Data'].setEnabled(True)
		self._mainview.buttons['Assign'].setEnabled(False)
		self._mainview.buttons['Select peak'].setEnabled(False)

		
		self._mainview.buttons['Select Elements'].clicked.connect(self._showPeriodicTable)
		
		self._mainview.buttons['Load ICP-MS Data'].clicked.connect(self._selectICPMSFile)
		self._mainview.buttons['Load ESI-MS Data'].clicked.connect(self._selectESIMSFile)

		self._mainview.buttons['Assign'].clicked.connect(self._assignFormulas)

		self._mainview.buttons['Select peak'].clicked.connect(self._getRange)
		
		self._mainview.replotBtn.setEnabled(False)

		self._mainview.replotBtn.clicked.connect(self._replotAll)
		
		self._mainview.offset.editingFinished.connect(self._getOffset)

		self._mainview.threshold.editingFinished.connect(self._getThreshold)

		self._mainview.setMouseTracking(True)

		
		#print(self._mainview.icp_checkBoxes.values())
		#for cbox in self._mainview.icp_checkBoxes.values():	
			

		#self._mainview.intbox.stateChanged.connect(self._selectIntRange)
		#self._mainview.oneFileBox.stateChanged.connect(self._selectOneFile)
		#self._mainview.baseSubtractBox.stateChanged.connect(self._baselineSubtraction)



