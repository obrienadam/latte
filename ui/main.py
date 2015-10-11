#!/usr/bin/env python2

import sys
from PyQt4 import QtGui, QtCore
from ui.main_window import Ui_MainWindow

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
        Bind widgets
        :return:
        """
        self.run_btn.clicked.connect(self.run)
        self.browse_btn.clicked.connect(self.browse_open_file)

    def run(self):
        solver_type = str(self.solver_combo_box.currentText())
        density = self.density_spin_box.value()
        mesh_file_name = self.mesh_line_edit.text()
        self.run_progress_bar.setValue(self.run_progress_bar.value() + 25)

        print 'Opening mesh {}'.format(mesh_file_name)

    def browse_open_file(self):
        file_name = QtGui.QFileDialog.getOpenFileName(QtGui.QFileDialog())
        self.mesh_line_edit.setText(file_name)

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    ui = Ui()
    ui.main_window.show()
    sys.exit(app.exec_())
