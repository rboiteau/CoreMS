import sys 
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtWidgets import * 
from pyqtgraph import PlotWidget, plot
import pyqtgraph as pg
from functools import partial
import os
import pandas as pd
from functools import partial
import json

from ui.PTBuilder.PTView import *
from ui.PTBuilder.PTCtrl import *
from ui.PTBuilder.PTModel import *

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
LC - ESI/ICP - MS data processing 

2022-08-22
Christian Dewey
'''


# Create a Controller class to connect the GUI and the model
class MainCtrl:
	"""Main Control class """
	def __init__(self, mainmodel, mainview, ptview):
		"""Controller initializer."""
		self._mainmodel = mainmodel
		self._mainview = mainview
		self._ptview = ptview
		
		self._peakRange = []		

		self.n_clicks = 0
		self._n = 0
		self._xMin = 0
		self._xMax = 0
		self.button_is_checked = False

		# Connect signals and slots
		self._connectSignals()

	def _buildExpression(self, sub_exp):
		"""Build expression."""
		expression = sub_exp
		self._mainview.setDisplayText(expression)

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
			self._connectSignals()

	def _selectICPMSFile(self):
		''' opens window to select normalization file for 115In correction; saves average 115In signal from norm file'''
		if self.n_clicks < 2:
			dialog = QFileDialog()
			dialog.setWindowTitle("Select ICP-MS File")
			dialog.setViewMode(QFileDialog.Detail)
			self._mainview.icpms_filepath = dialog.getOpenFileName(self._mainview,"Openfile")[0]
			print(self._mainview.icpms_filepath)

			'''imports LCICPMS .csv file'''
			self._mainview.icpms_data = pd.read_csv(self._mainview.icpms_filepath,sep=';|,',skiprows = 0, header = 0)
			print(self._mainview.icpms_data)
			icpms_header = list(self._mainview.icpms_data.columns.values)
			icpms_elements = []
			for v in icpms_header:
				if ('Time ' in v) or ('Number' in v) or (' ' in v):
					print(v)
					continue
				else:
					icpms_elements.append(v)
			
			self._mainview.all_icpms_elements = icpms_elements

			self._mainview._createICPElementCheckBoxes()
			self.n_clicks = self.n_clicks + 1 
			self._connectSignals()

	def _connectSignals(self):
		"""Connect signals and slots."""
		self._mainview.buttons['Select Elements'].clicked.connect(partial(self._buildExpression, ''))


		self._mainview.buttons['Select Elements'].setEnabled(True)
		self._mainview.buttons['Load ESI-MS Data'].setEnabled(True)
		self._mainview.buttons['Load ICP-MS Data'].setEnabled(True)
		self._mainview.buttons['Reset'].setEnabled(True)

		
		self._mainview.buttons['Select Elements'].clicked.connect(self._showPeriodicTable)

		print(self.n_clicks)
		
		self._mainview.buttons['Load ICP-MS Data'].clicked.connect(self._selectICPMSFile)
		self._mainview.buttons['Load ESI-MS Data'].clicked.connect(self._selectESIMSFile)
		
		for cbox in self._mainview.icp_checkBoxes.values():
			cbox.stateChanged.connect(partial(self._mainview.ICPMSClickBox, cbox))
		
		#print(self._mainview.icp_checkBoxes.values())
		#for cbox in self._mainview.icp_checkBoxes.values():	
			

		#self._mainview.intbox.stateChanged.connect(self._selectIntRange)
		#self._mainview.oneFileBox.stateChanged.connect(self._selectOneFile)
		#self._mainview.baseSubtractBox.stateChanged.connect(self._baselineSubtraction)



