import os
import sys
import importlib
cwd = os.getcwd()
plugin_path = os.path.dirname(os.path.abspath(__file__))

# print(f"[DEBUG] cwd: {cwd}")
# print(f"[DEBUG] plugin_path: {plugin_path}")
# print(f"[DEBUG] sys.path (before): {sys.path}")

# Ensure both cwd and plugin_path are in sys.path
if cwd not in sys.path:
    sys.path.insert(0, cwd)
if plugin_path not in sys.path:
    sys.path.insert(0, plugin_path)

# print(f"[DEBUG] sys.path (after): {sys.path}")

DEV_MODE = True

if DEV_MODE:
    import ui_layout
    importlib.reload(ui_layout)
    import logic
    importlib.reload(logic)

from PyQt5.QtWidgets import QDialog
from ui_layout import Ui_MainWindow
from logic import ProteinAnalyzerLogic

class ProteinAnalyzerPlugin(QDialog):
	def __init__(self):
		super().__init__()
		self.setWindowTitle("PDB2Graph")
		self.resize(600, 400)
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)
		self.logic = ProteinAnalyzerLogic(self.ui)


		# Wire up UI elements to centrality logic
		self.ui.centrality_button.clicked.connect(self.on_centrality_clicked)


	def on_centrality_clicked(self):
		metric = self.ui.centrality_combo.currentText()
		n = self.ui.centrality_topn.value()
		self.logic.apply_centrality(metric, n)

try:
    window = ProteinAnalyzerPlugin()
    window.show()
    window.raise_()
    window.activateWindow()
    print("[PDB2Graph] GUI launched via direct run.")
except Exception as e:
    print(f"[PDB2Graph] Failed to launch GUI: {e}")