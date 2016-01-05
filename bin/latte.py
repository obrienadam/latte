#!/usr/bin/env python2

import sys
from PyQt4 import QtGui, QtCore
from ui.main_window import Ui_MainWindow
from grid.finite_volume import FvRectilinearGrid
from solvers import simple, poisson

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

        grid = FvRectilinearGrid(shape, dimensions)

        bcs = {
            'type': [self.east_boundary_combo_box.currentText(),
                     self.north_boundary_combo_box.currentText(),
                     self.west_boundary_combo_box.currentText(),
                     self.south_boundary_combo_box.currentText()],
            'value': [self.east_boundary_spin_box.value(),
                      self.north_boundary_spin_box.value(),
                      self.west_boundary_spin_box.value(),
                      self.south_boundary_spin_box.value()]
        }

        solver_type = self.solver_combo_box.currentText()
        if solver_type == 'Poisson':
            solver_input = {'gamma': [],
                            'bcs': bcs,
                            'time_accurate': False,
                            }

            solver = poisson.Poisson(grid, **solver_input)

        elif solver_type == 'Simple':
            solver_input = {'max_iters': self.max_iters_spin_box.value(),
                            'time_accurate': self.unsteady_on_radio_btn.isChecked(),
                            'rho': self.density_spin_box.value(),
                            'mu': self.viscosity_spin_box.value()}

            solver = simple.Simple(grid, **solver_input)

        else:
            raise NotImplementedError

        print solver_input

    def browse_open_file(self):
        file_name = QtGui.QFileDialog.getOpenFileName(QtGui.QFileDialog())
        self.mesh_line_edit.setText(file_name)

if __name__ == '__main__':
    app = Latte(sys.argv)
