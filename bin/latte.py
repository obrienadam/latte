#!/usr/bin/env python2

import sys
from PyQt4 import QtGui, QtCore
from ui.main_window import Ui_MainWindow
from grid.finite_volume import FvGrid2D
from solvers.poisson import Poisson

class Latte(QtGui.QApplication, Ui_MainWindow):
    def __init__(self, args):
        """
        Initialize the GUI
        :return: The Latte appication object
        """
        QtGui.QApplication.__init__(self, args)
        self.main_window = QtGui.QMainWindow()
        self.setupUi(self.main_window)
        self.setupBindings()
        self.main_window.show()
        print 'Latte finished with exit code {}'.format(self.exec_())

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

        grid = FvGrid2D(shape, dimensions)

        boundaries = {'East': {'type': str(self.east_boundary_combo_box.currentText()), 'refval': self.east_boundary_spin_box.value()},
                      'West': {'type': str(self.west_boundary_combo_box.currentText()), 'refval': self.west_boundary_spin_box.value()},
                      'North': {'type': str(self.north_boundary_combo_box.currentText()), 'refval': self.north_boundary_spin_box.value()},
                      'South': {'type': str(self.south_boundary_combo_box.currentText()), 'refval': self.south_boundary_spin_box.value()}}

        solver_input = {'solver_type': str(self.solver_combo_box.currentText()),
                        'max_iters': self.max_iters_spin_box.value(),
                        'time_accurate': self.unsteady_on_radio_btn.isChecked()}

        print 'Solver input'
        print solver_input

        poisson = Poisson(grid, **solver_input)
        poisson.solve(self.run_progress_bar)
        self.run_progress_bar.setValue(0)

    def browse_open_file(self):
        file_name = QtGui.QFileDialog.getOpenFileName(QtGui.QFileDialog())
        self.mesh_line_edit.setText(file_name)

if __name__ == '__main__':
    app = Latte(sys.argv)
