# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'main_window.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
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
        MainWindow.resize(800, 600)
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
        self.solver_combo_box = QtGui.QComboBox(self.centralwidget)
        self.solver_combo_box.setObjectName(_fromUtf8("solver_combo_box"))
        self.solver_combo_box.addItem(_fromUtf8(""))
        self.solver_combo_box.addItem(_fromUtf8(""))
        self.formLayout.setWidget(0, QtGui.QFormLayout.FieldRole, self.solver_combo_box)
        self.solver_label = QtGui.QLabel(self.centralwidget)
        self.solver_label.setObjectName(_fromUtf8("solver_label"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.solver_label)
        self.density_label = QtGui.QLabel(self.centralwidget)
        self.density_label.setObjectName(_fromUtf8("density_label"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.density_label)
        self.viscosity_label = QtGui.QLabel(self.centralwidget)
        self.viscosity_label.setObjectName(_fromUtf8("viscosity_label"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.LabelRole, self.viscosity_label)
        self.density_spin_box = QtGui.QDoubleSpinBox(self.centralwidget)
        self.density_spin_box.setDecimals(4)
        self.density_spin_box.setSingleStep(0.1)
        self.density_spin_box.setProperty("value", 1.0)
        self.density_spin_box.setObjectName(_fromUtf8("density_spin_box"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.FieldRole, self.density_spin_box)
        self.viscosity_spin_box = QtGui.QDoubleSpinBox(self.centralwidget)
        self.viscosity_spin_box.setDecimals(4)
        self.viscosity_spin_box.setSingleStep(0.1)
        self.viscosity_spin_box.setProperty("value", 0.1)
        self.viscosity_spin_box.setObjectName(_fromUtf8("viscosity_spin_box"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.FieldRole, self.viscosity_spin_box)
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setObjectName(_fromUtf8("horizontalLayout_4"))
        self.mesh_line_edit = QtGui.QLineEdit(self.centralwidget)
        self.mesh_line_edit.setObjectName(_fromUtf8("mesh_line_edit"))
        self.horizontalLayout_4.addWidget(self.mesh_line_edit)
        self.browse_btn = QtGui.QPushButton(self.centralwidget)
        self.browse_btn.setObjectName(_fromUtf8("browse_btn"))
        self.horizontalLayout_4.addWidget(self.browse_btn)
        self.formLayout.setLayout(3, QtGui.QFormLayout.FieldRole, self.horizontalLayout_4)
        self.mesh_label = QtGui.QLabel(self.centralwidget)
        self.mesh_label.setObjectName(_fromUtf8("mesh_label"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.LabelRole, self.mesh_label)
        self.verticalLayout.addLayout(self.formLayout)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.verticalLayout.addItem(spacerItem)
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
        self.solution_tab = QtGui.QTabWidget(self.centralwidget)
        self.solution_tab.setObjectName(_fromUtf8("solution_tab"))
        self.mesh_tab = QtGui.QWidget()
        self.mesh_tab.setObjectName(_fromUtf8("mesh_tab"))
        self.verticalLayout_3 = QtGui.QVBoxLayout(self.mesh_tab)
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.solution_tab.addTab(self.mesh_tab, _fromUtf8(""))
        self.solution_tab1 = QtGui.QWidget()
        self.solution_tab1.setObjectName(_fromUtf8("solution_tab1"))
        self.solution_tab.addTab(self.solution_tab1, _fromUtf8(""))
        self.horizontalLayout.addWidget(self.solution_tab)
        self.horizontalLayout_2.addLayout(self.horizontalLayout)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 21))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.solution_tab.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "Lightweight Applications for Transport Equations", None))
        self.solver_combo_box.setItemText(0, _translate("MainWindow", "Simple", None))
        self.solver_combo_box.setItemText(1, _translate("MainWindow", "Piso", None))
        self.solver_label.setText(_translate("MainWindow", "Solver", None))
        self.density_label.setText(_translate("MainWindow", "Density", None))
        self.viscosity_label.setText(_translate("MainWindow", "Viscosity", None))
        self.browse_btn.setText(_translate("MainWindow", "Browse", None))
        self.mesh_label.setText(_translate("MainWindow", "Mesh File", None))
        self.run_btn.setText(_translate("MainWindow", "Run", None))
        self.solution_tab.setTabText(self.solution_tab.indexOf(self.mesh_tab), _translate("MainWindow", "Mesh", None))
        self.solution_tab.setTabText(self.solution_tab.indexOf(self.solution_tab1), _translate("MainWindow", "Solution", None))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    MainWindow = QtGui.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

