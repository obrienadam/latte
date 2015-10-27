# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'main_window.ui'
#
# Created: Tue Oct 27 15:00:52 2015
#      by: PyQt4 UI code generator 4.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(800, 620)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.horizontalLayout_2 = QtGui.QHBoxLayout(self.centralwidget)
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.formLayout = QtGui.QFormLayout()
        self.formLayout.setFieldGrowthPolicy(QtGui.QFormLayout.AllNonFixedFieldsGrow)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.solver_label = QtGui.QLabel(self.centralwidget)
        self.solver_label.setObjectName(_fromUtf8("solver_label"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.solver_label)
        self.solver_combo_box = QtGui.QComboBox(self.centralwidget)
        self.solver_combo_box.setObjectName(_fromUtf8("solver_combo_box"))
        self.solver_combo_box.addItem(_fromUtf8(""))
        self.solver_combo_box.addItem(_fromUtf8(""))
        self.formLayout.setWidget(0, QtGui.QFormLayout.FieldRole, self.solver_combo_box)
        self.density_label = QtGui.QLabel(self.centralwidget)
        self.density_label.setObjectName(_fromUtf8("density_label"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.density_label)
        self.density_spin_box = QtGui.QDoubleSpinBox(self.centralwidget)
        self.density_spin_box.setAccelerated(True)
        self.density_spin_box.setDecimals(4)
        self.density_spin_box.setSingleStep(0.1)
        self.density_spin_box.setProperty("value", 1.0)
        self.density_spin_box.setObjectName(_fromUtf8("density_spin_box"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.FieldRole, self.density_spin_box)
        self.viscosity_label = QtGui.QLabel(self.centralwidget)
        self.viscosity_label.setObjectName(_fromUtf8("viscosity_label"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.LabelRole, self.viscosity_label)
        self.viscosity_spin_box = QtGui.QDoubleSpinBox(self.centralwidget)
        self.viscosity_spin_box.setAccelerated(True)
        self.viscosity_spin_box.setDecimals(4)
        self.viscosity_spin_box.setSingleStep(0.1)
        self.viscosity_spin_box.setProperty("value", 0.1)
        self.viscosity_spin_box.setObjectName(_fromUtf8("viscosity_spin_box"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.FieldRole, self.viscosity_spin_box)
        self.max_iters_label = QtGui.QLabel(self.centralwidget)
        self.max_iters_label.setObjectName(_fromUtf8("max_iters_label"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.LabelRole, self.max_iters_label)
        self.max_iters_spin_box = QtGui.QSpinBox(self.centralwidget)
        self.max_iters_spin_box.setAccelerated(True)
        self.max_iters_spin_box.setSuffix(_fromUtf8(""))
        self.max_iters_spin_box.setMinimum(1)
        self.max_iters_spin_box.setMaximum(999999)
        self.max_iters_spin_box.setObjectName(_fromUtf8("max_iters_spin_box"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.FieldRole, self.max_iters_spin_box)
        self.unsteady_label = QtGui.QLabel(self.centralwidget)
        self.unsteady_label.setObjectName(_fromUtf8("unsteady_label"))
        self.formLayout.setWidget(4, QtGui.QFormLayout.LabelRole, self.unsteady_label)
        self.horizontalLayout_10 = QtGui.QHBoxLayout()
        self.horizontalLayout_10.setObjectName(_fromUtf8("horizontalLayout_10"))
        self.unsteady_on_radio_btn = QtGui.QRadioButton(self.centralwidget)
        self.unsteady_on_radio_btn.setChecked(False)
        self.unsteady_on_radio_btn.setObjectName(_fromUtf8("unsteady_on_radio_btn"))
        self.horizontalLayout_10.addWidget(self.unsteady_on_radio_btn)
        self.unsteady_off_radio_btn = QtGui.QRadioButton(self.centralwidget)
        self.unsteady_off_radio_btn.setChecked(True)
        self.unsteady_off_radio_btn.setObjectName(_fromUtf8("unsteady_off_radio_btn"))
        self.horizontalLayout_10.addWidget(self.unsteady_off_radio_btn)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_10.addItem(spacerItem)
        self.formLayout.setLayout(4, QtGui.QFormLayout.FieldRole, self.horizontalLayout_10)
        self.verticalLayout.addLayout(self.formLayout)
        self.line_2 = QtGui.QFrame(self.centralwidget)
        self.line_2.setFrameShape(QtGui.QFrame.HLine)
        self.line_2.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_2.setObjectName(_fromUtf8("line_2"))
        self.verticalLayout.addWidget(self.line_2)
        self.formLayout_2 = QtGui.QFormLayout()
        self.formLayout_2.setObjectName(_fromUtf8("formLayout_2"))
        self.mesh_label = QtGui.QLabel(self.centralwidget)
        self.mesh_label.setObjectName(_fromUtf8("mesh_label"))
        self.formLayout_2.setWidget(2, QtGui.QFormLayout.LabelRole, self.mesh_label)
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.mesh_line_edit = QtGui.QLineEdit(self.centralwidget)
        self.mesh_line_edit.setObjectName(_fromUtf8("mesh_line_edit"))
        self.horizontalLayout_3.addWidget(self.mesh_line_edit)
        self.browse_btn = QtGui.QPushButton(self.centralwidget)
        self.browse_btn.setObjectName(_fromUtf8("browse_btn"))
        self.horizontalLayout_3.addWidget(self.browse_btn)
        self.formLayout_2.setLayout(2, QtGui.QFormLayout.FieldRole, self.horizontalLayout_3)
        self.mesh_resolution_label = QtGui.QLabel(self.centralwidget)
        self.mesh_resolution_label.setObjectName(_fromUtf8("mesh_resolution_label"))
        self.formLayout_2.setWidget(0, QtGui.QFormLayout.LabelRole, self.mesh_resolution_label)
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setObjectName(_fromUtf8("horizontalLayout_4"))
        self.mesh_xresolution_spin_box = QtGui.QSpinBox(self.centralwidget)
        self.mesh_xresolution_spin_box.setAccelerated(True)
        self.mesh_xresolution_spin_box.setMinimum(1)
        self.mesh_xresolution_spin_box.setMaximum(10000)
        self.mesh_xresolution_spin_box.setProperty("value", 11)
        self.mesh_xresolution_spin_box.setObjectName(_fromUtf8("mesh_xresolution_spin_box"))
        self.horizontalLayout_4.addWidget(self.mesh_xresolution_spin_box)
        self.mesh_yresolution_spin_box = QtGui.QSpinBox(self.centralwidget)
        self.mesh_yresolution_spin_box.setAccelerated(True)
        self.mesh_yresolution_spin_box.setMinimum(1)
        self.mesh_yresolution_spin_box.setMaximum(10000)
        self.mesh_yresolution_spin_box.setProperty("value", 11)
        self.mesh_yresolution_spin_box.setObjectName(_fromUtf8("mesh_yresolution_spin_box"))
        self.horizontalLayout_4.addWidget(self.mesh_yresolution_spin_box)
        self.formLayout_2.setLayout(0, QtGui.QFormLayout.FieldRole, self.horizontalLayout_4)
        self.mesh_dimensions_label = QtGui.QLabel(self.centralwidget)
        self.mesh_dimensions_label.setObjectName(_fromUtf8("mesh_dimensions_label"))
        self.formLayout_2.setWidget(1, QtGui.QFormLayout.LabelRole, self.mesh_dimensions_label)
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setObjectName(_fromUtf8("horizontalLayout_5"))
        self.mesh_xdimensions_spin_box = QtGui.QDoubleSpinBox(self.centralwidget)
        self.mesh_xdimensions_spin_box.setAccelerated(True)
        self.mesh_xdimensions_spin_box.setDecimals(4)
        self.mesh_xdimensions_spin_box.setSingleStep(0.1)
        self.mesh_xdimensions_spin_box.setProperty("value", 1.0)
        self.mesh_xdimensions_spin_box.setObjectName(_fromUtf8("mesh_xdimensions_spin_box"))
        self.horizontalLayout_5.addWidget(self.mesh_xdimensions_spin_box)
        self.mesh_ydimensions_spin_box = QtGui.QDoubleSpinBox(self.centralwidget)
        self.mesh_ydimensions_spin_box.setAccelerated(True)
        self.mesh_ydimensions_spin_box.setDecimals(4)
        self.mesh_ydimensions_spin_box.setMinimum(0.0)
        self.mesh_ydimensions_spin_box.setSingleStep(0.1)
        self.mesh_ydimensions_spin_box.setProperty("value", 1.0)
        self.mesh_ydimensions_spin_box.setObjectName(_fromUtf8("mesh_ydimensions_spin_box"))
        self.horizontalLayout_5.addWidget(self.mesh_ydimensions_spin_box)
        self.formLayout_2.setLayout(1, QtGui.QFormLayout.FieldRole, self.horizontalLayout_5)
        self.verticalLayout.addLayout(self.formLayout_2)
        self.line_3 = QtGui.QFrame(self.centralwidget)
        self.line_3.setFrameShape(QtGui.QFrame.HLine)
        self.line_3.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_3.setObjectName(_fromUtf8("line_3"))
        self.verticalLayout.addWidget(self.line_3)
        self.formLayout_3 = QtGui.QFormLayout()
        self.formLayout_3.setObjectName(_fromUtf8("formLayout_3"))
        self.east_boundary_label = QtGui.QLabel(self.centralwidget)
        self.east_boundary_label.setObjectName(_fromUtf8("east_boundary_label"))
        self.formLayout_3.setWidget(0, QtGui.QFormLayout.LabelRole, self.east_boundary_label)
        self.horizontalLayout_6 = QtGui.QHBoxLayout()
        self.horizontalLayout_6.setObjectName(_fromUtf8("horizontalLayout_6"))
        self.east_boundary_combo_box = QtGui.QComboBox(self.centralwidget)
        self.east_boundary_combo_box.setObjectName(_fromUtf8("east_boundary_combo_box"))
        self.east_boundary_combo_box.addItem(_fromUtf8(""))
        self.east_boundary_combo_box.addItem(_fromUtf8(""))
        self.east_boundary_combo_box.addItem(_fromUtf8(""))
        self.horizontalLayout_6.addWidget(self.east_boundary_combo_box)
        self.east_boundary_spin_box = QtGui.QDoubleSpinBox(self.centralwidget)
        self.east_boundary_spin_box.setAccelerated(True)
        self.east_boundary_spin_box.setDecimals(4)
        self.east_boundary_spin_box.setSingleStep(0.1)
        self.east_boundary_spin_box.setObjectName(_fromUtf8("east_boundary_spin_box"))
        self.horizontalLayout_6.addWidget(self.east_boundary_spin_box)
        self.formLayout_3.setLayout(0, QtGui.QFormLayout.FieldRole, self.horizontalLayout_6)
        self.west_boundary_label = QtGui.QLabel(self.centralwidget)
        self.west_boundary_label.setObjectName(_fromUtf8("west_boundary_label"))
        self.formLayout_3.setWidget(1, QtGui.QFormLayout.LabelRole, self.west_boundary_label)
        self.horizontalLayout_7 = QtGui.QHBoxLayout()
        self.horizontalLayout_7.setObjectName(_fromUtf8("horizontalLayout_7"))
        self.west_boundary_combo_box = QtGui.QComboBox(self.centralwidget)
        self.west_boundary_combo_box.setObjectName(_fromUtf8("west_boundary_combo_box"))
        self.west_boundary_combo_box.addItem(_fromUtf8(""))
        self.west_boundary_combo_box.addItem(_fromUtf8(""))
        self.west_boundary_combo_box.addItem(_fromUtf8(""))
        self.horizontalLayout_7.addWidget(self.west_boundary_combo_box)
        self.west_boundary_spin_box = QtGui.QDoubleSpinBox(self.centralwidget)
        self.west_boundary_spin_box.setAccelerated(True)
        self.west_boundary_spin_box.setDecimals(4)
        self.west_boundary_spin_box.setSingleStep(0.1)
        self.west_boundary_spin_box.setProperty("value", 1.0)
        self.west_boundary_spin_box.setObjectName(_fromUtf8("west_boundary_spin_box"))
        self.horizontalLayout_7.addWidget(self.west_boundary_spin_box)
        self.formLayout_3.setLayout(1, QtGui.QFormLayout.FieldRole, self.horizontalLayout_7)
        self.north_boundary_label = QtGui.QLabel(self.centralwidget)
        self.north_boundary_label.setObjectName(_fromUtf8("north_boundary_label"))
        self.formLayout_3.setWidget(2, QtGui.QFormLayout.LabelRole, self.north_boundary_label)
        self.horizontalLayout_8 = QtGui.QHBoxLayout()
        self.horizontalLayout_8.setObjectName(_fromUtf8("horizontalLayout_8"))
        self.north_boundary_combo_box = QtGui.QComboBox(self.centralwidget)
        self.north_boundary_combo_box.setObjectName(_fromUtf8("north_boundary_combo_box"))
        self.north_boundary_combo_box.addItem(_fromUtf8(""))
        self.north_boundary_combo_box.addItem(_fromUtf8(""))
        self.north_boundary_combo_box.addItem(_fromUtf8(""))
        self.horizontalLayout_8.addWidget(self.north_boundary_combo_box)
        self.north_boundary_spin_box = QtGui.QDoubleSpinBox(self.centralwidget)
        self.north_boundary_spin_box.setAccelerated(True)
        self.north_boundary_spin_box.setDecimals(4)
        self.north_boundary_spin_box.setSingleStep(0.1)
        self.north_boundary_spin_box.setObjectName(_fromUtf8("north_boundary_spin_box"))
        self.horizontalLayout_8.addWidget(self.north_boundary_spin_box)
        self.formLayout_3.setLayout(2, QtGui.QFormLayout.FieldRole, self.horizontalLayout_8)
        self.south_boundary_label = QtGui.QLabel(self.centralwidget)
        self.south_boundary_label.setObjectName(_fromUtf8("south_boundary_label"))
        self.formLayout_3.setWidget(3, QtGui.QFormLayout.LabelRole, self.south_boundary_label)
        self.horizontalLayout_9 = QtGui.QHBoxLayout()
        self.horizontalLayout_9.setObjectName(_fromUtf8("horizontalLayout_9"))
        self.south_boundary_combo_box = QtGui.QComboBox(self.centralwidget)
        self.south_boundary_combo_box.setObjectName(_fromUtf8("south_boundary_combo_box"))
        self.south_boundary_combo_box.addItem(_fromUtf8(""))
        self.south_boundary_combo_box.addItem(_fromUtf8(""))
        self.south_boundary_combo_box.addItem(_fromUtf8(""))
        self.horizontalLayout_9.addWidget(self.south_boundary_combo_box)
        self.south_boundary_spin_box = QtGui.QDoubleSpinBox(self.centralwidget)
        self.south_boundary_spin_box.setAccelerated(True)
        self.south_boundary_spin_box.setDecimals(4)
        self.south_boundary_spin_box.setSingleStep(0.1)
        self.south_boundary_spin_box.setObjectName(_fromUtf8("south_boundary_spin_box"))
        self.horizontalLayout_9.addWidget(self.south_boundary_spin_box)
        self.formLayout_3.setLayout(3, QtGui.QFormLayout.FieldRole, self.horizontalLayout_9)
        self.verticalLayout.addLayout(self.formLayout_3)
        self.run_btn = QtGui.QPushButton(self.centralwidget)
        self.run_btn.setObjectName(_fromUtf8("run_btn"))
        self.verticalLayout.addWidget(self.run_btn)
        self.run_progress_bar = QtGui.QProgressBar(self.centralwidget)
        self.run_progress_bar.setProperty("value", 0)
        self.run_progress_bar.setAlignment(QtCore.Qt.AlignBottom|QtCore.Qt.AlignHCenter)
        self.run_progress_bar.setInvertedAppearance(False)
        self.run_progress_bar.setObjectName(_fromUtf8("run_progress_bar"))
        self.verticalLayout.addWidget(self.run_progress_bar)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.line = QtGui.QFrame(self.centralwidget)
        self.line.setFrameShape(QtGui.QFrame.VLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.horizontalLayout.addWidget(self.line)
        self.Tab = QtGui.QTabWidget(self.centralwidget)
        self.Tab.setObjectName(_fromUtf8("Tab"))
        self.mesh_tab = QtGui.QWidget()
        self.mesh_tab.setObjectName(_fromUtf8("mesh_tab"))
        self.Tab.addTab(self.mesh_tab, _fromUtf8(""))
        self.run_tab = QtGui.QWidget()
        self.run_tab.setObjectName(_fromUtf8("run_tab"))
        self.Tab.addTab(self.run_tab, _fromUtf8(""))
        self.solution_tab = QtGui.QWidget()
        self.solution_tab.setObjectName(_fromUtf8("solution_tab"))
        self.Tab.addTab(self.solution_tab, _fromUtf8(""))
        self.horizontalLayout.addWidget(self.Tab)
        self.horizontalLayout_2.addLayout(self.horizontalLayout)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 27))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)
        self.actionOpen = QtGui.QAction(MainWindow)
        self.actionOpen.setObjectName(_fromUtf8("actionOpen"))
        self.menuFile.addAction(self.actionOpen)
        self.menubar.addAction(self.menuFile.menuAction())

        self.retranslateUi(MainWindow)
        self.east_boundary_combo_box.setCurrentIndex(1)
        self.north_boundary_combo_box.setCurrentIndex(2)
        self.south_boundary_combo_box.setCurrentIndex(2)
        self.Tab.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "Lightweight Applications for Transport Equations", None))
        self.solver_label.setText(_translate("MainWindow", "Solver", None))
        self.solver_combo_box.setItemText(0, _translate("MainWindow", "Simple", None))
        self.solver_combo_box.setItemText(1, _translate("MainWindow", "Piso", None))
        self.density_label.setText(_translate("MainWindow", "Density", None))
        self.viscosity_label.setText(_translate("MainWindow", "Viscosity", None))
        self.max_iters_label.setText(_translate("MainWindow", "Max. Iterations", None))
        self.unsteady_label.setText(_translate("MainWindow", "Unsteady", None))
        self.unsteady_on_radio_btn.setText(_translate("MainWindow", "On", None))
        self.unsteady_off_radio_btn.setText(_translate("MainWindow", "Off", None))
        self.mesh_label.setText(_translate("MainWindow", "Mesh File", None))
        self.browse_btn.setText(_translate("MainWindow", "Browse", None))
        self.mesh_resolution_label.setText(_translate("MainWindow", "Mesh Resolution", None))
        self.mesh_dimensions_label.setText(_translate("MainWindow", "Mesh Dimensions", None))
        self.east_boundary_label.setText(_translate("MainWindow", "East Boundary", None))
        self.east_boundary_combo_box.setItemText(0, _translate("MainWindow", "Inlet", None))
        self.east_boundary_combo_box.setItemText(1, _translate("MainWindow", "Outlet", None))
        self.east_boundary_combo_box.setItemText(2, _translate("MainWindow", "Wall", None))
        self.west_boundary_label.setText(_translate("MainWindow", "West Boundary", None))
        self.west_boundary_combo_box.setItemText(0, _translate("MainWindow", "Inlet", None))
        self.west_boundary_combo_box.setItemText(1, _translate("MainWindow", "Outlet", None))
        self.west_boundary_combo_box.setItemText(2, _translate("MainWindow", "Wall", None))
        self.north_boundary_label.setText(_translate("MainWindow", "North Boundary", None))
        self.north_boundary_combo_box.setItemText(0, _translate("MainWindow", "Inlet", None))
        self.north_boundary_combo_box.setItemText(1, _translate("MainWindow", "Outlet", None))
        self.north_boundary_combo_box.setItemText(2, _translate("MainWindow", "Wall", None))
        self.south_boundary_label.setText(_translate("MainWindow", "South Boundary", None))
        self.south_boundary_combo_box.setItemText(0, _translate("MainWindow", "Inlet", None))
        self.south_boundary_combo_box.setItemText(1, _translate("MainWindow", "Outlet", None))
        self.south_boundary_combo_box.setItemText(2, _translate("MainWindow", "Wall", None))
        self.run_btn.setText(_translate("MainWindow", "Run", None))
        self.Tab.setTabText(self.Tab.indexOf(self.mesh_tab), _translate("MainWindow", "Mesh", None))
        self.Tab.setTabText(self.Tab.indexOf(self.run_tab), _translate("MainWindow", "Run", None))
        self.Tab.setTabText(self.Tab.indexOf(self.solution_tab), _translate("MainWindow", "Solution", None))
        self.menuFile.setTitle(_translate("MainWindow", "File", None))
        self.actionOpen.setText(_translate("MainWindow", "Open", None))
        self.actionOpen.setShortcut(_translate("MainWindow", "Ctrl+O", None))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    MainWindow = QtGui.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

