import sys 
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import * 
from pyqtgraph import PlotWidget, plot
import pyqtgraph as pg
from functools import partial
import os
import pandas as pd
from functools import partial

__version__ = '0.1'
__author__ = 'Christian Dewey'

'''
LCICPMS data GUI

2022-04-21
Christian Dewey
'''

# Create a subclass of QMainWindow to setup the periodic table gui
class PTView(QWidget):
	
	def __init__(self,mainview):
		"""View initializer."""
		super().__init__()
		# Set some main window's properties
		self._mainview = mainview
		self._elementLabels = []
    	# self.setLayout(mainLayout)
		self.setWindowTitle('Periodic Table')
		self.setGeometry(100, 60, 800, 600)
		# Set the central widget
		self.generalLayout = QVBoxLayout()
		self.topLayout = QVBoxLayout()
		self.bottomLayout = QHBoxLayout()
		self._centralWidget = QWidget(self)

		self.setLayout(self.generalLayout)

		self.generalLayout.addLayout(self.topLayout)
		self.generalLayout.addLayout(self.bottomLayout)

		self._createPeriodicTable()
		self._createButtons()

	def _createButtons(self):
		"""Create the buttons."""
		self.buttons = {}
		buttonsLayout = QGridLayout()
		# Button text | position on the QGridLayout
		buttons = {'Save': (0, 0),
			'Reset': (0, 1)}
		# Create the buttons and add them to the grid layout
		for btnText, pos in buttons.items():
			self.buttons[btnText] = QPushButton(btnText)
			self.buttons[btnText].setFixedSize(100, 40)
			buttonsLayout.addWidget(self.buttons[btnText], pos[0], pos[1])
		# Add buttonsLayout to the general layout
		self.bottomLayout.addLayout(buttonsLayout)

	def _createResizeHandle(self):
		handle = QSizeGrip(self)
		#self.generalLayout.addWidget(handle)
		self.generalLayout.addWidget(handle, 0, Qt.AlignBottom | Qt.AlignRight)
	   # self.__corner = Qt.BottomRightCorner

		self.resize(self.sizeHint())

	def _createPeriodicTable(self):
		"""Create the buttons."""
		self.periodicTable = {}
		ptLayout = QGridLayout()
		# Button text | position on the QGridLayout

		# Create the buttons and add them to the grid layout
		for element, attr in self._mainview.periodicTableDict.items():
			self.periodicTable[element] = QPushButton(element)
			self.periodicTable[element].setFixedSize(50, 50)
			if attr[3] == 0:
				self.periodicTable[element].setStyleSheet('background-color : ' + attr[2])
			elif attr[3] == 1:
				self.periodicTable[element].setStyleSheet('background-color : lightgray')
			ptLayout.addWidget(self.periodicTable[element], attr[0], attr[1])
		# Add buttonsLayout to the general layout
		self.topLayout.addLayout(ptLayout)
	
	def clicked(self):
		item = self.listwidget.currentItem()
		print('file: ' + item.text())
		return self.listwidget.currentItem()
	
	def setDisplayText(self, text):
		"""Set display's text."""
		self.display.setText(text)
		self.display.setFocus()

	def displayText(self):
		"""Get display's text."""
		return self.display.text()

	def clearChecks(self):
		"""Clear the display."""
		for cbox in self.checkBoxes:
			cbox.setCheckState(Qt.Unchecked)

	def clickBox(self, cbox, state):
		if state == Qt.Checked:
			print('checked: ' + cbox.text())
			self.activeMetals.append(cbox.text())
		   # print(self.activeMetals)
			return self.activeMetals
		elif state == Qt.Unchecked:
			print('Unchecked: ' + cbox.text())
			self.activeMetals.remove(cbox.text())
			#print(self.activeMetals)
			return self.activeMetals
		else:
			print('Unchecked')
			return self.activeMetals


