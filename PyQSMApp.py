# -*- coding: utf-8 -*-
# Created on Mar 28, 2020
# Author: Maxwell4444
# All right reserved
from PyQSM.PyQSMGui import QSMClass
from PyQt5 import QtWidgets
import sys

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    qsmapp = QSMClass()
    qsmapp.show()
    sys.exit(app.exec_())