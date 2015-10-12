#!/usr/bin/env python2

import sys
from PyQt4 import QtGui, QtCore
from ui.main_window import Ui_MainWindow
import solvers
from grid.finite_difference import FdGrid2D

class Ui(Ui_MainWindow):
    def __init__(self):
        """
        Initialize the GUI
        :return:
        """
        self.main_window = QtGui.QMainWindow()
        self.setupUi(self.main_window)
        self.setupBindings()

    def setupBindings(self):
        """
        Bind the buttons to the apropriate function calls
        :return:
        """
        self.run_btn.clicked.connect(self.run)
        self.browse_btn.clicked.connect(self.browse_open_file)

    def run(self):
        solver_type = str(self.solver_combo_box.currentText())
        shape = self.mesh_xresolution_spin_box.value(), self.mesh_yresolution_spin_box.value()
        dimensions = self.mesh_xdimensions_spin_box.value(), self.mesh_ydimensions_spin_box.value()

        grid = FdGrid2D(shape, dimensions)

        boundaries = {'East': {'type': str(self.east_boundary_combo_box.currentText()), 'refval': self.east_boundary_spin_box.value()},
                      'West': {'type': str(self.west_boundary_combo_box.currentText()), 'refval': self.west_boundary_spin_box.value()},
                      'North': {'type': str(self.north_boundary_combo_box.currentText()), 'refval': self.north_boundary_spin_box.value()},
                      'South': {'type': str(self.south_boundary_combo_box.currentText()), 'refval': self.south_boundary_spin_box.value()}}

        grid.boundaries = boundaries

        print grid.boundaries

    def browse_open_file(self):
        file_name = QtGui.QFileDialog.getOpenFileName(QtGui.QFileDialog())
        self.mesh_line_edit.setText(file_name)

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    ui = Ui()
    ui.main_window.show()
    sys.exit(app.exec_())
