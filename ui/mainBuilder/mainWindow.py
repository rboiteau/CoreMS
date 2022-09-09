import sys
from PyQt5.QtCore import Qt, QAbstractTableModel, QVariant
from PyQt5.QtWidgets import * 
from pyqtgraph import PlotWidget, plot
import pyqtgraph as pg
from functools import partial
import os
import re
import pandas as pd
from functools import partial

from ui.mainBuilder.mainChroma import *

from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters

__version__ = '0.1'
__author__ = 'Christian Dewey'

'''
LCICPMS data GUI

2022-04-21
Christian Dewey
'''


class TableModel(QAbstractTableModel):
	def __init__(self, data, parent=None):
		QAbstractTableModel.__init__(self, parent)
		self._data = data
	def headerData(self, section, orientation, role=Qt.DisplayRole):
		if orientation == Qt.Horizontal and role == Qt.DisplayRole:
			header = list(self._data.columns)
			header[1] = 'R2'
			return  header[section]
			#return 'Column {}'.format(section + 1)
		return super().headerData(section, orientation, role)
	def rowCount(self, parent=None):
		return len(self._data.values)

	def columnCount(self, parent=None):
		return self._data.columns.size

	def data(self, index, role=Qt.DisplayRole):
		if index.isValid():
			if role == Qt.DisplayRole:
				return QVariant(str(
					self._data.iloc[index.row()][index.column()]))
		return QVariant()


# Create a subclass of QMainWindow to setup the calculator's GUI
class MainView(QMainWindow):
	
	def __init__(self):
		"""View initializer."""
		super().__init__()
		# Set some main window's properties
		self.setWindowTitle('LC-ICP/ESI-MS Data Analysis')
		self.setGeometry(100, 60, 1000, 600)
		# Set the central widget
		self.generalLayout = QHBoxLayout()
		self.topLayout = QFormLayout()
		self.plotLayout = QVBoxLayout()
		self.assignLayout = QGridLayout()
		self.ICPMS_layout = QVBoxLayout()
		self.ESIMS_layout = QVBoxLayout()
		self.EIC_layout = QVBoxLayout()



		#self._ElementLayout = QGridLayout()
		self._ReportLayout = QVBoxLayout()
		#self._ParamsLayout = QGridLayout()
		self._centralWidget = QWidget(self)

		self.setCentralWidget(self._centralWidget)
		self._centralWidget.setLayout(self.generalLayout)

		self.generalLayout.addLayout(self.plotLayout)
		self.generalLayout.addLayout(self.assignLayout)

		self._activeElementList = ['C', 'H', 'O']
		self._currentElements = []
		self.icpms_filepath = None # '/Users/christiandewey/CoreMS/tests/tests_data/icpms/161220_soils_hypercarb_3_kansas_qH2O.csv'
		self.icpms_data = None
		self.icpms_elements_for_plotting = []
		self.all_icpms_elements = []
		
		self.esims_filepath = None #'/Users/christiandewey/CoreMS/tests/tests_data/ftms/rmb_161221_kansas_h2o_2.raw'
		self.esi_parser = None

		self.firstSave = True

		self.esi_file_selected = False
		self.icp_file_selected = False

		self._requiredElements = ['C', 'H', 'O']

		self._elementListFlag = False
		self.eic_plotted = False 

		self.best_results = None
		self.mzs = None
		self.offset = None
		'''
		## for testing 
		MSParameters.mass_spectrum.threshold_method = 'signal_noise'
		MSParameters.mass_spectrum.s2n_threshold = 6
		MSParameters.ms_peak.peak_min_prominence_percent = 0.1
		self.esi_parser = rawFileReader.ImportMassSpectraThermoMSFileReader(self.esims_filepath)

		self.icpms_data = pd.read_csv(self.icpms_filepath,sep=';|,',skiprows = 0, header = 0)
		icpms_header = list(self.icpms_data.columns.values)
		icpms_elements = []
		for v in icpms_header:
			if ('Time ' in v) or ('Number' in v) or (' ' in v):
				print(v)
				continue
			else:
				icpms_elements.append(v)
		
		self.all_icpms_elements = icpms_elements

		self._createICPElementCheckBoxes()
		## end of testing section 
		'''


		self.filepath = ''
		self.homeDir = '' #/Users/christiandewey/'# '/Users/christiandewey/presentations/DOE-PI-22/day6/day6/'
		self.activeMetals = []
		
		self.icp_checkBoxes = {}

		self.singleOutputFile = False		
		self.baseSubtract = False 

		self._createPTDict()

		self._createButtons()
		
		#self._createDisplay()
		self._createICPMSPlot()
		self._createEICPlot()
		self._createESIMSPlot()
		self._createResizeHandle()

		self._elementBox = QGroupBox()
		self._elementBox.setMaximumWidth(200)
		self._paramBox = QGroupBox()
		self._paramBox.setMaximumSize(200,100)
		self._plotParamBox = QGroupBox()
		self._plotParamBox.setMaximumSize(200,120)


		self.offset = QLineEdit()
		self.offset_float = 0.0
		self.threshold = QLineEdit()
		self.threshold_float = 0.0

		self.replotBtn = QPushButton('Replot')

		self.maxx = QLineEdit()
		self.maxx.setText('40')
		self.minx = QLineEdit()
		self.minx.setText('0')

		self.assignedResults = QTableView()

	def _createResizeHandle(self):
		handle = QSizeGrip(self)
		#self.generalLayout.addWidget(handle)
		self.generalLayout.addWidget(handle, 0, Qt.AlignBottom | Qt.AlignRight)
	   # self.__corner = Qt.BottomRightCorner

		self.resize(self.sizeHint())

	def clearLayout(self, layout):
		if layout is not None:
			while layout.count():
				item = layout.takeAt(0)
				widget = item.widget()
				if widget is not None:
					widget.deleteLater()
				else:
					self.clearLayout(item.layout())

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
		self.ICPMS_PlotSpace = pg.PlotWidget()
		self.ICPMS_PlotSpace.setBackground('w')
		styles = { 'font-size':'15px'}
		self.ICPMS_PlotSpace.setLabel('left', 'ICP-MS signal (1000s cps)', **styles)
		self.ICPMS_PlotSpace.setLabel('bottom', "Retention time (min)", **styles)
		self.icpms_chroma = self.ICPMS_PlotSpace
		self.ICPMS_layout.addWidget(self.ICPMS_PlotSpace)
		#self.generalLayout.addLayout(self.ICPMS_layout)
		self.plotLayout.addLayout(self.ICPMS_layout)

	def _createESIMSPlot(self):
		self.ESIMS_PlotSpace = pg.PlotWidget()
		self.ESIMS_PlotSpace.setBackground('w')
		styles = { 'font-size':'15px'}
		self.ESIMS_PlotSpace.setLabel('left', 'ESI-MS signal (1000s cps)', **styles)
		self.ESIMS_PlotSpace.setLabel('bottom', "m/z", **styles)
		self.esims = self.ESIMS_PlotSpace
		self.ESIMS_layout.addWidget(self.ESIMS_PlotSpace)
		#self.generalLayout.addLayout(self.ESIMS_layout)
		#self.plotLayout.addLayout(self.ESIMS_layout)

	def _createEICPlot(self):
		self.EIC_PlotSpace = pg.PlotWidget()
		self.EIC_PlotSpace.setBackground('w')
		styles = { 'font-size':'15px'}
		self.EIC_PlotSpace.setLabel('left', 'ESI-MS signal (1000s cps)', **styles)
		self.EIC_PlotSpace.setLabel('bottom', "Retention time (min)", **styles)
		self.eic = self.EIC_PlotSpace
		self.EIC_layout.addWidget(self.EIC_PlotSpace)
		#self.generalLayout.addLayout(self.EIC_layout)
		self.plotLayout.addLayout(self.EIC_layout)

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
	
	def _createAssignBoxes(self):	
		self.assignLayout.addWidget(self._elementBox, 1 , 0)
		self.assignLayout.addWidget(self._paramBox, 0, 0)
		self.assignLayout.addLayout(self._ReportLayout, 2 , 0)
		self.assignLayout.addWidget(self._plotParamBox, 3 , 0)

		
	def _createICPElementCheckBoxes(self):

		icp_checkBoxes = {} 
		optionsLayout = QHBoxLayout()
		self.icp_checkbox_layout = optionsLayout

		if self._activeElementList == []:
			icpms_elements = self.all_icpms_elements
			for m in icpms_elements:
				cbox = QCheckBox(m)
				cbox.setChecked(True)
				icp_checkBoxes[m] = cbox
				optionsLayout.addWidget(cbox)

		else:
			icpms_elements = []
			#print(self._activeElementList)
			for e in self.all_icpms_elements:
				match = re.findall(r'[A-Za-z]+|\d+', e)
			#	print(match[1])
				if match[1] in self._activeElementList:
					icpms_elements.append(e)
			
			if '115In' in self.all_icpms_elements:
				icpms_elements.append('115In')

			for m in icpms_elements:
				#print(m)
				cbox = QCheckBox(m)
				cbox.setChecked(True)
				icp_checkBoxes[m] = cbox
				optionsLayout.addWidget(cbox)
		self.icp_checkBoxes = icp_checkBoxes
		self.icpms_elements_for_plotting = icpms_elements
		self.activeMetals = icpms_elements

		if self.icp_file_selected is True:	
			self._makeICPMSPlot()

		optionsLayout.setAlignment(Qt.AlignLeft)
		self.ICPMS_layout.addLayout(optionsLayout)


	def _updateICPMS_elementList(self):
		
		icpms_elements = []
		
		for e in self.all_icpms_elements:
			match = re.findall(r'[A-Za-z]+|\d+', e)
			#print(match[1])
			if match[1] in self._activeElementList:
				icpms_elements.append(e)
		

		#for m in icpms_elements:
		#	print(m)
		#	cbox = QCheckBox(m)
		#	cbox.setChecked(True)
		#	icp_checkBoxes[m] = cbox
		#	optionsLayout.addWidget(cbox)
		#self.icp_checkBoxes = icp_checkBoxes
		self.icpms_elements_for_plotting = icpms_elements

		for e in self.icp_checkBoxes.keys():
			
			self.icp_checkBoxes[e].setParent(None)

		self.icp_checkBoxes = {}

		self._createICPElementCheckBoxes()

		for cbox in self.icp_checkBoxes.values():
			cbox.stateChanged.connect(partial(self.ICPMSClickBox, cbox))

	def _createElementOrder(self, elements):
	
		sortOrder = ['C', 'H', 'O', 'N', 'P', 'S', 'Cl', 'Na', 'Br']
		addedN = False
		addedP = False
		addedS = False
		addedCl = False
		addedNa = False
		addedBr = False

		#elements = self._activeElementList

		for e in elements:
			if e == 'N':
				addedN = True
			elif e == 'P':
				addedP = True
			elif e == 'S':
				addedS = True
			elif e == 'Cl':
				addedCl = True
			elif e == 'Na':
				addedNa = True
			elif e == 'Br':
				addedBr = True
			
		orderTF = [True, True, True, addedN, addedP, addedS, addedCl, addedNa, addedBr]     

		orderSub = [val for is_good, val in zip(orderTF, sortOrder) if is_good]

		otherElements = [el for el in elements if el not in sortOrder ]

		finalOrder = orderSub + otherElements
		
		return finalOrder

	def _assignHeteroAtom(self, element,it):

		print(it)


		try: 

			max_signal = 0
			
			dominant_isotope = None

			available_isotopes = []

			print(element)
			
			for e in self.all_icpms_elements:
				
				match = re.findall(r'[A-Za-z]+|\d+', e)
				
				if match[1] == element:
					
					imax = max(self.icpms_data[e])

					available_isotopes.append(e)

					if imax > max_signal:		

						dominant_isotope = match[0]

						heteroatom = dominant_isotope + element
			
			strlist = ' '.join(map(str,available_isotopes))

			if heteroatom:
		
				print('\nAvailable isotopes in ICPMS data: ' + strlist)

				print('Most abundant isotope in ICPMS data: ' + heteroatom)

				self.heteroatom = heteroatom

		except:

			print('\nError assigning dominant isotope. Is the ICPMS data loaded?')

	def _createAssignForm_active(self):

	#if self._elementListFlag is False:

		it = 0 

		elements = self._activeElementList
		
		assignOrder = self._createElementOrder(elements)

		self._elementMins = {}
		self._elementMaxs = {}
		self._elementLabs = {}
		self._elementRadios = {}

		#self._elementBox = QGroupBox()
		self._eBoxLayout = QGridLayout()
		layout = self._eBoxLayout
		minLabel = QLabel('Min')
		maxLabel = QLabel('Max')

		layout.addWidget(minLabel,0,1)
		layout.addWidget(maxLabel,0,2)

		commonElements = ['C', 'H', 'O', 'N', 'P', 'S', 'Cl', 'Na', 'Br']
		defaultMins = [1, 4, 1, 0, 0, 0, 0, 0, 0]
		defaultMaxs = [50, 100, 20, 4, 0, 0, 0, 0, 0]
		defaultMinMax = {e: (min1, max1) for e, min1, max1 in zip(commonElements, defaultMins, defaultMaxs)}

		for element, row in zip(assignOrder, range(0,len(assignOrder))):
			self._elementMins[element] = QLineEdit()
			self._elementMins[element].setFixedSize(40, 20)
			self._elementMaxs[element] = QLineEdit()
			self._elementMaxs[element].setFixedSize(40, 20)
			

			if element in defaultMinMax.keys():
				default_min = str(defaultMinMax[element][0])
				default_max = str(defaultMinMax[element][1])
			else:
				default_min = '0'
				default_max = '1'

			self._elementMins[element].setText(default_min)
			self._elementMaxs[element].setText(default_max)
			elMin = self._elementMins[element]
			elMax = self._elementMaxs[element]
			label = QLabel(element)
			self._elementLabs[element] = label
			layout.addWidget(label, row+1, 0)
			layout.addWidget(elMin, row+1, 1)
			layout.addWidget(elMax, row+1, 2)

			chon_elements = ['C', 'H', 'O', 'N']

			print(it)

			if element not in chon_elements:
				self._elementRadios[element] = QRadioButton()
				layout.addWidget(self._elementRadios[element], row + 1 , 3)
				self._elementRadios[element].pressed.connect(partial(self._assignHeteroAtom,element, it))
				print(element)
				

		layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)

		self._elementBox.setLayout(layout)

		
		#self.assignLayout.addWidget(self._elementBox)

		#self._elementListFlag = True
		
		for e in self._activeElementList:
			self._currentElements.append(e)


		self.param_layout = QGridLayout()

		param_layout = self.param_layout

		self.offset.setText('0')
		self.offset.setFixedSize(40, 20)
		self.offsetLabel = QLabel('Offset:')  

		self.threshold.setText('0.8')
		self.threshold.setFixedSize(40, 20)
		self.thresholdLabel = QLabel('Min R2:')
		
		param_layout.addWidget(self.offsetLabel, 0, 0)
		param_layout.addWidget(self.offset, 0, 1)

		param_layout.addWidget(self.thresholdLabel, 1, 0)
		param_layout.addWidget(self.threshold, 1, 1)

		self._paramBox.setLayout(param_layout)

		print('here')


		#self.assignLayout.addWidget(self._paramBox)
	


			
	def _updateElementList(self):

		new_elements = self._activeElementList

		newAssignOrder = self._createElementOrder(new_elements)
		old_elements = self._currentElements
		elements_to_hide = np.setdiff1d(old_elements, new_elements)
		elements_to_add = np.setdiff1d(new_elements, old_elements)
		
		commonElements = ['C', 'H', 'O', 'N', 'P', 'S', 'Cl', 'Na', 'Br']
		defaultMins = [1, 4, 1, 0, 0, 0, 0, 0, 0]
		defaultMaxs = [50, 100, 20, 4, 0, 0, 0, 0, 0]
		defaultMinMax = {e: (min1, max1) for e, min1, max1 in zip(commonElements, defaultMins, defaultMaxs)}

		#self._elementBox.setParent(None)

		layout = self._eBoxLayout

		for element, row in zip(self._currentElements, range(0,len(self._currentElements))):
			
			if (element in elements_to_hide) and (element not in self._requiredElements):
				self._elementMins[element].setParent(None)
				del self._elementMins[element]
				self._elementLabs[element].setParent(None)
				del self._elementLabs[element]
				self._elementMaxs[element].setParent(None)
				del self._elementMaxs[element]
			else:
				pass

		for element, row in zip(newAssignOrder, range(0,len(newAssignOrder))):
			
			if (element in self._currentElements) and (element not in elements_to_add) and (element not in elements_to_hide):

				self._elementMins[element].setParent(None)
				del self._elementMins[element]
				self._elementLabs[element].setParent(None)
				del self._elementLabs[element]
				self._elementMaxs[element].setParent(None)
				del self._elementMaxs[element]

				self._elementMins[element] = QLineEdit()
				self._elementMins[element].setFixedSize(40, 20)
				self._elementMaxs[element] = QLineEdit()
				self._elementMaxs[element].setFixedSize(40, 20)

				if element in defaultMinMax.keys():
					default_min = str(defaultMinMax[element][0])
					default_max = str(defaultMinMax[element][1])
				else:
					default_min = '0'
					default_max = '1'

				self._elementMins[element].setText(default_min)
				self._elementMaxs[element].setText(default_max)
				elMin = self._elementMins[element]
				elMax = self._elementMaxs[element]
				label = QLabel(element)
				self._elementLabs[element] = label
				layout.addWidget(label, row+1, 0)
				layout.addWidget(elMin, row+1, 1)
				layout.addWidget(elMax, row+1, 2)

				#self._elementBox.setLayout(layout)

				#self.assignLayout.addWidget(self._elementBox)

			elif element in elements_to_add:
				self._elementMins[element] = QLineEdit()
				self._elementMins[element].setFixedSize(40, 20)
				self._elementMaxs[element] = QLineEdit()
				self._elementMaxs[element].setFixedSize(40, 20)

				if element in defaultMinMax.keys():
					default_min = str(defaultMinMax[element][0])
					default_max = str(defaultMinMax[element][1])
				else:
					default_min = '0'
					default_max = '1'

				self._elementMins[element].setText(default_min)
				self._elementMaxs[element].setText(default_max)
				elMin = self._elementMins[element]
				elMax = self._elementMaxs[element]
				label = QLabel(element)
				self._elementLabs[element] = label
				layout.addWidget(label, row+1, 0)
				layout.addWidget(elMin, row+1, 1)
				layout.addWidget(elMax, row+1, 2)

				#self._elementBox.setLayout(layout)

				#self.assignLayout.addWidget(self._elementBox)

			else:
				pass 


		self._currentElements.clear()
		for e in self._activeElementList:
				self._currentElements.append(e)


		#self._paramBox.setParent(None)

		param_layout = self.param_layout
		temp_offset = self.offset.text()
		self.offset.setParent(None)
		self.offsetLabel.setParent(None)
		#del self.offsetLabel
		#del self.offset
		
		temp_threshold = self.threshold.text()
		self.threshold.setParent(None)
		self.thresholdLabel.setParent(None)
		#del self.thresholdLabel
		#del self.threshold

		self.offset.setText(temp_offset)
		self.offset.setFixedSize(40, 20)
		self.offsetLabel = QLabel('Offset:')  

		self.threshold.setText(temp_threshold)
		self.threshold.setFixedSize(40, 20)
		self.thresholdLabel = QLabel('Min R2:')
		
		param_layout.addWidget(self.offsetLabel, 0, 0)
		param_layout.addWidget(self.offset, 0, 1)

		param_layout.addWidget(self.thresholdLabel, 1, 0)
		param_layout.addWidget(self.threshold, 1, 1)

		#self._paramBox.setLayout(param_layout)

		#self.assignLayout.addWidget(self._paramBox)


	def _plotParams(self):
		 
		layout2 = QGridLayout()
		'''
		self.offset = QLineEdit()
		self.offset.setText('0')
		self.offset.setFixedSize(40, 20)
		self.offsetLabel = QLabel('Offset:')  

		layout.addWidget(self.offsetLabel, 0, 0)
		layout.addWidget(self.offset, 0, 1)

		self.threshold = QLineEdit()
		self.threshold.setText('0.8')
		self.threshold.setFixedSize(40, 20)
		self.thresholdLabel = QLabel('Min R2:')		

		layout.addWidget(self.thresholdLabel, 1, 0)
		layout.addWidget(self.threshold, 1, 1)
		'''

		self.minx.setText('0')
		self.minx.setFixedSize(40, 20)
		self.minxLabel = QLabel('X min:')		
		
		layout2.addWidget(self.minxLabel, 0, 0)
		layout2.addWidget(self.minx, 0, 1)

		self.maxx.setText('40')
		self.maxx.setFixedSize(40, 20)
		self.maxxLabel = QLabel('X max:')		

		layout2.addWidget(self.maxxLabel, 1, 0)
		layout2.addWidget(self.maxx, 1, 1)

		#self.replotBtn = QPushButton('Replot')
		self.replotBtn.setFixedSize(100,40)
		layout2.addWidget(self.replotBtn,2,0)

		self._plotParamBox.setLayout(layout2)
		#self.assignLayout.addLayout(layout)


	def _createButtons(self):
		"""Create the buttons."""
		self.buttons = {}
		buttonsLayout = QGridLayout()
		# Button text | position on the QGridLayout
		buttons = {'Select Elements': (0, 0, 200),
				   'Load ESI-MS Data': (0, 1, 200),
				   'Load ICP-MS Data': (0, 2, 200),
				   'Select peak': (0, 3, 200),
				   'Assign': (0, 4, 70)
				  }
		# Create the buttons and add them to the grid layout
		for btnText, pos in buttons.items():
			self.buttons[btnText] = QPushButton(btnText)
			self.buttons[btnText].setFixedSize(pos[2], 40)
			buttonsLayout.addWidget(self.buttons[btnText], pos[0], pos[1])
		# Add buttonsLayout to the general layout
		self.plotLayout.addLayout(buttonsLayout)
	
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

	def _makeTICPlot(self):
		self.tic_plot = plotChroma(self, self.icpms_elements_for_plotting, self.icpms_data, self.activeMetals)._plotTIC(self.esi_parser)

	def _makeEICPlot(self,mzs):
		self.eic_plotted = True 
		self.eic_plot = plotChroma(self, self.icpms_elements_for_plotting, self.icpms_data, self.activeMetals)._plotEIC(self.esi_parser,mzs)

	def ICPMSClickBox(self, cbox, state):
		if state == Qt.Checked:
			print('checked: ' + cbox.text())
			if cbox.text() not in self.activeMetals:
				self.activeMetals.append(cbox.text())
				self._makeICPMSPlot()
			#	print(self.activeMetals)
				#return self.activeMetals
		elif state == Qt.Unchecked:
			print('Unchecked: ' + cbox.text())
		#	print(self.activeMetals)
			self.activeMetals.remove(cbox.text())
			if self.activeMetals == []:
				self.ICPMS_PlotSpace.clear()
			else:
				self._makeICPMSPlot()
			
		else:
			print('Unchecked')

	def generateResultsReport(self):
		
		nrows = np.shape(self.best_results)[0]
		print(nrows)
		print(np.shape(self.best_results))
		print(self.best_results.iloc[0,:])

		bestResults_dis = self.best_results.iloc[:,:-1]

		for r in range(nrows): # m/z formatting
			bestResults_dis.iloc[r,0] = "{:.5f}".format(self.best_results.iloc[r,0])

		for r in range(nrows): # corr formatting
			bestResults_dis.iloc[r,1] = "{:.4f}".format(self.best_results.iloc[r,1])
		
		for r in range(nrows): # peak height formatting
			bestResults_dis.iloc[r,2] = "{:.0f}".format(self.best_results.iloc[r,2])
	
		for r in range(nrows): # confidence formatting
			bestResults_dis.iloc[r,3] = "{:.4f}".format(self.best_results.iloc[r,3])

		model = TableModel(bestResults_dis)
		self.assignedResults.setModel(model)
		widths = [100, 70, 100, 100, 150]

		for i,w in zip(range(model.columnCount()),widths):
			self.assignedResults.setColumnWidth(i, w)

		#for c in range(0,ncols):
		#	for r in range(0,nrows):
		#		if r == 0:
		#			self.assignedResults.setItem(0,c, QTableWidgetItem(header[c]))
		##			self.assignedResults.setItem(r,c, QTableWidgetItem(self.best_results.iloc[r-1,c]))

		self._ReportLayout.addWidget(self.assignedResults)






