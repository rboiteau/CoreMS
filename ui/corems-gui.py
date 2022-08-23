import sys 
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import * 
from pyqtgraph import PlotWidget, plot
import pyqtgraph as pg
from functools import partial
import os
import pandas as pd
from functools import partial

from mainBuilder.mainWindow import *
from mainBuilder.mainCtrl import *
from mainBuilder.mainModel import *

from PTBuilder.PTView import *
from PTBuilder.PTCtrl import *
from PTBuilder.PTModel import *

# Client code
def main():
    """Main function."""
    # Create an instance of QApplication
    main_app = QApplication(sys.argv)
    # Show the calculator's GU
    mainView = MainView()
    mainModel = MainModel(mainview = mainView)
    ptView = PTView(mainview = mainView)
    ptModel = PTModel(mainview = mainView, ptview = ptView)
    
    mainView.show()

    # Create instances of the model and the controller
    MainCtrl(mainmodel=mainModel, mainview=mainView, ptview=ptView)
    ptCtrl = PTCtrl(ptview = ptView, model = ptModel, mainview = mainView, mainctrl= MainCtrl)
    # Execute the main loop
    if (sys.flags.interactive != 1) or not hasattr(Qt, 'PYQT_VERSION'):
        main_app.exec_()
    #sys.exit(pycalc.exec_())

if __name__ == '__main__':
    main()