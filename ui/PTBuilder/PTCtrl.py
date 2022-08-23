from posixpath import split
import sys 
from PyQt5.QtCore import Qt, QTimer, QCoreApplication
from PyQt5.QtWidgets import * 
from pyqtgraph import PlotWidget, plot
import pyqtgraph as pg
from functools import partial
import os
import pandas as pd
from functools import partial
import json

from PTBuilder.PTView import *
from PTBuilder.PTCtrl import *
from PTBuilder.PTModel import *

__version__ = '0.1'
__author__ = 'Christian Dewey'

'''
LCICPMS data GUI

2022-04-21
Christian Dewey
'''
from PTBuilder.PTView import *

# Create a Controller class to connect the GUI and the model
class PTCtrl:
	"""Main Control class """
	def __init__(self, model, mainview, ptview, mainctrl):
		"""Controller initializer."""
		self._model = model
		self._view = ptview
		self._mainview = mainview
		self._mainctrl = mainctrl
		self._activeElementList = self._mainview._activeElementList
		
		# Connect signals and slots
		self._connectSignals()

	def _saveElementList(self):
		self._mainview._activeElementList = self._activeElementList
		#QCoreApplication.instance().quit()
		self._view.close()

	def _resetPeriodicTable(self):
		self._mainview._activeElementList = []
		for element, btn in self._view.periodicTable.items():
			col = self._mainview.periodicTableDict[element][2]
			self._view.periodicTable[element].setStyleSheet('background-color : '+ col)
			self._mainview.periodicTableDict[element][3] = 0

	def _alterElementList(self, element):
		
		split_el = element.split('\n')[1]

		isActive = self._mainview.periodicTableDict[element][3]

		if isActive == 0:
			
			self._mainview.periodicTableDict[element][3] = 1
			#self._view._elementLabels.append(element)
			self._view.periodicTable[element].setStyleSheet('background-color : lightgray')
			self._mainview._activeElementList.append(split_el)
			print(self._mainview._activeElementList)

		elif isActive == 1:

			self._mainview.periodicTableDict[element][3] = 0
			#self._view._elementLabels.remove(element)
			col = self._mainview.periodicTableDict[element][2]
			self._view.periodicTable[element].setStyleSheet('background-color : ' + col)
			self._mainview._activeElementList.remove(split_el)
			print(self._mainview._activeElementList)

	def _connectSignals(self):
		"""Connect signals and slots."""
		for btnText, btn in self._view.periodicTable.items():
			if btnText  in {'H'}:
				if self._view.listwidget.currentItem() is None:
					text = ''
				else:
					text = self._view.listwidget.currentItem().text()
		
				btn.clicked.connect(partial(self._buildExpression, text))

		for element, btn in self._view.periodicTable.items():
			self._view.periodicTable[element].clicked.connect(partial(self._alterElementList, element))

		self._view.buttons['Reset'].clicked.connect(self._resetPeriodicTable)
		self._view.buttons['Save'].clicked.connect(self._saveElementList)

