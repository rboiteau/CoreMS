import sys 
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import * 
from pyqtgraph import PlotWidget, plot
import pyqtgraph as pg
from functools import partial
import os
import re
import pandas as pd
from functools import partial

from .mainChroma import *

__version__ = '0.1'
__author__ = 'Christian Dewey'

'''
LCICPMS data GUI

2022-04-21
Christian Dewey
'''

# Create a subclass of QMainWindow to setup the calculator's GUI
class MainView(QMainWindow):
	
	def __init__(self):
		"""View initializer."""
		super().__init__()
		# Set some main window's properties
		self.setWindowTitle('LC-ICP/ESI-MS Data Analysis')
		self.setGeometry(100, 60, 1000, 600)
		# Set the central widget
		self.generalLayout = QVBoxLayout()
		self.topLayout = QFormLayout()
		self.ICPMS_layout = QVBoxLayout()
		self._centralWidget = QWidget(self)
		self.setCentralWidget(self._centralWidget)
		self._centralWidget.setLayout(self.generalLayout)

		self._activeElementList = []
		self.icpms_filepath = None
		self.icpms_data = None
		self.icpms_elements_for_plotting = []
		self.all_icpms_elements = []
		self.esims_filepath = None
		
		self.filepath = ''
		self.normAvIndium = -999.99
		self.homeDir = '' #/Users/christiandewey/'# '/Users/christiandewey/presentations/DOE-PI-22/day6/day6/'
		self.activeMetals = []
		
		self.icp_checkBoxes = {}

		self.singleOutputFile = False		
		self.baseSubtract = False 

		self._createPTDict()

		self._createButtons()
		
		self._createDisplay()
		self._createICPMSPlot()
		self._createIntegrateCheckBoxes()
		self._createResizeHandle()

	def _createResizeHandle(self):
		handle = QSizeGrip(self)
		#self.generalLayout.addWidget(handle)
		self.generalLayout.addWidget(handle, 0, Qt.AlignBottom | Qt.AlignRight)
	   # self.__corner = Qt.BottomRightCorner

		self.resize(self.sizeHint())

	def _createPTDict(self):
		self.periodicTableDict = {'1\nH': [0, 0, 'salmon',0],
				'2\nHe': [0, 18, 'salmon',0],
				'3\nLi': [1,0, 'salmon',0],
				'4\nBe': [1, 1, 'salmon',0],
				'11\nNa': [2, 0, 'salmon',0],
				'12\nMg': [2, 1, 'salmon',0],
				'19\nK': [3, 0, 'salmon',0],
				'20\nCa': [3, 1, 'salmon',0],
				'37\nRb': [4, 0, 'salmon',0],
				'38\nSr': [4, 1, 'salmon',0],
				'55\nCs': [5, 0, 'salmon',0],
				'56\nBa': [5, 1, 'salmon',0],
				'87\nFr': [6, 0, 'salmon',0],
				'88\nRa': [6, 1, 'salmon',0],
				'21\nSc': [3, 3, 'lightblue',0],
				'22\nTi': [3, 4, 'lightblue',0],
				'23\nV': [3, 5, 'lightblue',0],
				'24\nCr': [3, 6, 'lightblue',0],
				'25\nMn': [3, 7, 'lightblue',0],
				'26\nFe': [3, 8, 'lightblue',0],
				'27\nCo': [3, 9, 'lightblue',0],
				'28\nNi': [3, 10, 'lightblue',0],
				'29\nCu': [3, 11, 'lightblue',0],
				'30\nZn': [3, 12, 'lightblue',0],
				'39\nY': [4, 3, 'lightblue',0],  ##
				'40\nZr': [4, 4, 'lightblue',0],
				'41\nNb': [4, 5, 'lightblue',0],
				'42\nMo': [4, 6, 'lightblue',0],
				'43\nTc': [4, 7, 'lightblue',0],
				'44\nRu': [4, 8, 'lightblue',0],
				'45\nRh': [4, 9, 'lightblue',0],
				'46\nPd': [4, 10, 'lightblue',0],
				'47\nAg': [4, 11, 'lightblue',0],
				'48\nCd': [4, 12, 'lightblue',0],
				'71\nLu': [5, 3, 'lightblue',0],  ##
				'72\nHf': [5, 4, 'lightblue',0],
				'73\nTa': [5, 5, 'lightblue',0],
				'74\nW': [5, 6, 'lightblue',0],
				'75\nRe': [5, 7, 'lightblue',0],
				'76\nOs': [5, 8, 'lightblue',0],
				'77\nIr': [5, 9, 'lightblue',0],
				'78\nPt': [5, 10, 'lightblue',0],
				'79\nAu': [5, 11, 'lightblue',0],
				'80\nHg': [5, 12, 'lightblue',0],
				'103\nLr': [6, 3, 'lightblue',0],  ##
				'104\nRf': [6, 4, 'lightblue',0],
				'105\nDb': [6, 5, 'lightblue',0],
				'106\nSg': [6,6, 'lightblue',0],
				'107\nBh': [6, 7, 'lightblue',0],
				'108\nHs': [6, 8, 'lightblue',0],
				'109\nMt': [6, 9, 'lightblue',0],
				'110\nDs': [6, 10, 'lightblue',0],
				'111\nRg': [6, 11, 'lightblue',0],
				'112\nCn': [6, 12, 'lightblue',0],
				'5\nB': [1, 13, 'green',0],  ##
				'6\nC': [1, 14, 'green',0],
				'7\nN': [1, 15, 'green',0],
				'8\nO': [1,16, 'green',0],
				'9\nF': [1, 17, 'green',0],
				'10\nNe': [1, 18, 'yellow',0],
				'13\nAl': [2, 13, 'green',0],  ##
				'14\nSi': [2, 14, 'green',0],
				'15\nP': [2, 15, 'green',0],
				'16\nS': [2, 16, 'green',0],
				'17\nCl': [2, 17, 'green',0],
				'18\nAr': [2, 18, 'yellow',0],
				'31\nGa': [3, 13, 'green',0],  ##
				'32\nGe': [3, 14, 'green',0],
				'33\nAs': [3, 15, 'green',0],
				'34\nSe': [3, 16, 'green',0],
				'35\nBr': [3, 17, 'green',0],
				'36\nKr': [3, 18, 'yellow',0],
				'49\nIn': [4, 13, 'green',0],  ##
				'50\nSn': [4, 14, 'green',0],
				'51\nSb': [4, 15, 'green',0],
				'52\nTe': [4, 16, 'green',0],
				'53\nI': [4, 17, 'green',0],
				'54\nXe': [4, 18, 'yellow',0],
				'81\nTl': [5, 13, 'green',0],  ##
				'82\nPb': [5, 14, 'green',0],
				'83\nBi': [5, 15, 'green',0],
				'84\nPo': [5, 16, 'green',0],
				'85\nAt': [5, 17, 'green',0],
				'86\nRn': [5, 18, 'yellow',0],
				'113\nNh': [6, 13, 'green',0],  ##
				'114\nFl': [6, 14, 'green',0],
				'115\nMc': [6, 15, 'green',0],
				'116\nLv': [6, 16, 'green',0],
				'117\nTs': [6, 17, 'green',0],
				'118\nOg': [6, 18, 'yellow',0],
				'57\nLa': [7, 2, 'orange',0], ##
				'58\nCe': [7, 3, 'orange',0],
				'59\nPr': [7, 4, 'orange',0],  ##
				'60\nNd': [7, 5, 'orange',0],
				'61\nPm': [7, 6, 'orange',0],
				'62\nSm': [7, 7, 'orange',0],
				'63\nEu': [7, 8, 'orange',0],
				'64\nGd': [7, 9, 'orange',0],
				'65\nTb': [7, 10, 'orange',0],  ##
				'66\nDy': [7, 11, 'orange',0],
				'67\nHo': [7, 12, 'orange',0],
				'68\nEr': [7, 13, 'orange',0],
				'69\nTm': [7, 14, 'orange',0],
				'70\nYb': [7, 15, 'orange',0],
				'89\nAc': [8, 2, 'orange',0], ##
				'90\nTh': [8, 3, 'orange',0],
				'91\nPa': [8, 4, 'orange',0],  ##
				'92\nU': [8, 5, 'orange',0],
				'93\nNp': [8, 6, 'orange',0],
				'94\nPu': [8, 7, 'orange',0],
				'95\nAm': [8, 8, 'orange',0],
				'96\nCm': [8, 9, 'orange',0],
				'97\nBk': [8, 10, 'orange',0],  ##
				'98\nCf': [8, 11, 'orange',0],
				'99\nEs': [8, 12, 'orange',0],
				'100\nFm': [8, 13, 'orange',0],
				'101\nMd': [8, 14, 'orange',0],
				'102\nNo': [8, 15, 'orange',0]
			}

	def _createICPMSPlot(self):
		self.plotSpace = pg.PlotWidget()
		self.plotSpace.setBackground('w')
		styles = { 'font-size':'15px'}
		self.plotSpace.setLabel('left', 'ICP-MS signal (1000s cps)', **styles)
		self.plotSpace.setLabel('bottom', "Retention time (min)", **styles)
		self.icpms_chroma = self.plotSpace
		self.ICPMS_layout.addWidget(self.plotSpace)
		self.generalLayout.addLayout(self.ICPMS_layout)

	def _createESIMSPlot(self):
		self.plotSpace = pg.PlotWidget()
		self.plotSpace.setBackground('w')
		styles = { 'font-size':'15px'}
		self.plotSpace.setLabel('left', 'ICP-MS signal (1000s cps)', **styles)
		self.plotSpace.setLabel('bottom', "Retention time (min)", **styles)
		self.icpms_chroma = self.plotSpace
		self.ICPMS_layout.addWidget(self.plotSpace)
		self.generalLayout.addLayout(self.ICPMS_layout)

	def _createEICPlot(self):
		self.plotSpace = pg.PlotWidget()
		self.plotSpace.setBackground('w')
		styles = { 'font-size':'15px'}
		self.plotSpace.setLabel('left', 'ESI-MS signal (1000s cps)', **styles)
		self.plotSpace.setLabel('bottom', "Retention time (min)", **styles)
		self.icpms_chroma = self.plotSpace
		self.EIC_layout.addWidget(self.plotSpace)
		self.generalLayout.addLayout(self.ESI_layout)

	def _createDirEntry(self):
		self.DirEntry = QLineEdit()
		self.DirEntry.setFixedHeight(35)
		self.DirEntry.setAlignment(Qt.AlignRight)
		self.topLayout.addRow("Enter directory:", self.DirEntry)
		self.topLayout.addWidget(self.DirEntry)

	def _createDisplay(self):
		'''Create the display'''
		# Create the display widget
		self.display = QLineEdit()
		self.display.setFixedHeight(35)
		self.display.setAlignment(Qt.AlignRight)
		self.display.setReadOnly(True)
		self.generalLayout.addWidget(self.display)

	def _createICPElementCheckBoxes(self):
		icp_checkBoxes = {} 
		optionsLayout = QHBoxLayout()

		if self._activeElementList == []:
			icpms_elements = self.all_icpms_elements
			for m in icpms_elements:
				cbox = QCheckBox(m)
				icp_checkBoxes[m] = cbox
				optionsLayout.addWidget(cbox)
		# optionwidget.stateChanged.connect(self.clickBox)

		else:
			icpms_elements = []
			print(self._activeElementList)
			for e in self.all_icpms_elements:
				match = re.findall(r'[A-Za-z]+|\d+', e)
				print(match[1])
				if match[1] in self._activeElementList:
					icpms_elements.append(e)
			
			if '115In' in self.all_icpms_elements:
				icpms_elements.append('115In')

			for m in icpms_elements:
				print(m)
				cbox = QCheckBox(m)
				icp_checkBoxes[m] = cbox
				optionsLayout.addWidget(cbox)
		self.icp_checkBoxes = icp_checkBoxes
		self.icpms_elements_for_plotting = icpms_elements
		self.ICPMS_layout.addLayout(optionsLayout)

	def _createIntegrateCheckBoxes(self):
		# Add some checkboxes to the layout  
		#self.integrateBox= []      
		self.integrateLayout = QHBoxLayout()
		checkboxLayout =QVBoxLayout()
		self.intbox = QCheckBox('Select integration range?')
		self.oneFileBox = QCheckBox('Single output file?')
		self.baseSubtractBox = QCheckBox('Baseline subtraction?')
		checkboxLayout.addWidget(self.intbox)
		checkboxLayout.addWidget(self.oneFileBox)
		checkboxLayout.addWidget(self.baseSubtractBox)
		self.integrateLayout.addLayout(checkboxLayout)




	def _createButtons(self):
		"""Create the buttons."""
		self.buttons = {}
		buttonsLayout = QGridLayout()
		# Button text | position on the QGridLayout
		buttons = {'Select Elements': (0, 0, 200),
				   'Load ESI-MS Data': (0, 1, 200),
				   'Load ICP-MS Data': (0, 2, 200),
				   'Reset': (0,3,70)
				  }
		# Create the buttons and add them to the grid layout
		for btnText, pos in buttons.items():
			self.buttons[btnText] = QPushButton(btnText)
			self.buttons[btnText].setFixedSize(pos[2], 40)
			buttonsLayout.addWidget(self.buttons[btnText], pos[0], pos[1])
		# Add buttonsLayout to the general layout
		self.generalLayout.addLayout(buttonsLayout)
	
	def clicked(self):
		item = self.listwidget.currentItem()
		print('\nfile: ' + item.text())
		return self.listwidget.currentItem()
	
	def setDisplayText(self, text):
		"""Set display's text."""
		self.display.setText(text)
		self.display.setFocus()

	def clearChecks(self):
		"""Clear the display."""
		for cbox in self.checkBoxes.values():
			cbox.setCheckState(Qt.Unchecked)

	def _makeICPMSPlot(self):
		self.icpms_chroma = plotChroma(self, self.icpms_elements_for_plotting, self.icpms_data, self.activeMetals)._plotICPMSChroma()

	def ICPMSClickBox(self, cbox, state):
		if state == Qt.Checked:
			print('checked: ' + cbox.text())
			if cbox.text() not in self.activeMetals:
				self.activeMetals.append(cbox.text())
				self._makeICPMSPlot()
				print(self.activeMetals)
				#return self.activeMetals
		elif state == Qt.Unchecked:
			print('Unchecked: ' + cbox.text())
			print(self.activeMetals)
			self.activeMetals.remove(cbox.text())
			if self.activeMetals == []:
				self.plotSpace.clear()
			else:
				self._makeICPMSPlot()
			
		else:
			print('Unchecked')

