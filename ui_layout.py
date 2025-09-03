# 
# from PyQt5.QtWidgets import (
#     QVBoxLayout, QHBoxLayout, QGridLayout, QPushButton,
#     QLabel, QComboBox, QWidget, QGroupBox
# )
# 
# class Ui_MainWindow:
#     def setupUi(self, Dialog):
#         main_layout = QVBoxLayout(Dialog)
# 
#         # Section 1: Protein Structure Input and Network Generation
#         box1 = QGroupBox("Generate Protein Structure Network:")
#         layout1 = QVBoxLayout()
#         self.load_button = QPushButton("Load protein structures")
#         self.select_structure = QComboBox()
#         self.create_network = QPushButton("Create network")
#         layout1.addWidget(self.load_button)
#         layout1.addWidget(self.select_structure)
#         layout1.addWidget(self.create_network)
#         box1.setLayout(layout1)
# 
#         # Section 2: Visualization Options
#         box2 = QGroupBox("Select Visualization:")
#         layout2 = QVBoxLayout()
#         self.view_in_protein = QPushButton("View network in protein")
#         self.view_as_network = QPushButton("View as network")
#         layout2.addWidget(self.view_in_protein)
#         layout2.addWidget(self.view_as_network)
#         box2.setLayout(layout2)
# 
#         # Section 3: Node Centrality
#         box3 = QGroupBox("Select a Network Centrality:")
#         layout3 = QVBoxLayout()
#         self.node_betweenness = QPushButton("Node betweenness")
#         self.node_closeness = QPushButton("Node closeness")
#         self.node_degree = QPushButton("Node degree")
#         self.node_eigenvector = QPushButton("Node eigenvector")
#         layout3.addWidget(self.node_betweenness)
#         layout3.addWidget(self.node_closeness)
#         layout3.addWidget(self.node_degree)
#         layout3.addWidget(self.node_eigenvector)
#         box3.setLayout(layout3)
# 
#         # Section 4: Edge Analysis
#         box4 = QGroupBox("Edge Analysis:")
#         layout4 = QVBoxLayout()
#         self.edge_betweenness = QPushButton("Edge betweenness")
#         self.edge_weight = QPushButton("Edge weight")
#         self.edge_weight_ebc = QPushButton("Edge weight * EBC")
#         self.reset_edges = QPushButton("Reset edges")
#         layout4.addWidget(self.edge_betweenness)
#         layout4.addWidget(self.edge_weight)
#         layout4.addWidget(self.edge_weight_ebc)
#         layout4.addWidget(self.reset_edges)
#         box4.setLayout(layout4)
# 
#         # Section 5: Dashboard & Comparison
#         box5 = QGroupBox("Guide for customizing the visualization:")
#         layout5 = QVBoxLayout()
#         self.generate_dashboard = QPushButton("Generate Dashboard")
#         self.protein1 = QComboBox()
#         self.protein2 = QComboBox()
#         self.unique_edges1 = QPushButton("Unique edge protein 1")
#         self.unique_edges2 = QPushButton("Unique edges protein 2")
#         layout5.addWidget(self.generate_dashboard)
#         layout5.addWidget(self.protein1)
#         layout5.addWidget(self.protein2)
#         layout5.addWidget(self.unique_edges1)
#         layout5.addWidget(self.unique_edges2)
#         box5.setLayout(layout5)
# 
#         # Layout grouping
#         top_row = QHBoxLayout()
#         top_row.addWidget(box1)
#         top_row.addWidget(box2)
# 
#         middle_row = QHBoxLayout()
#         middle_row.addWidget(box3)
#         middle_row.addWidget(box4)
# 
#         main_layout.addLayout(top_row)
#         main_layout.addLayout(middle_row)
#         main_layout.addWidget(box5)


from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QGridLayout, QPushButton,
    QLabel, QComboBox, QWidget, QGroupBox, QSpinBox
)

class Ui_MainWindow:
    def setupUi(self, Dialog):
        main_layout = QVBoxLayout(Dialog)

        # Section 1: Protein Structure Input and Network Generation
        box1 = QGroupBox("Generate Protein Structure Network:")
        layout1 = QVBoxLayout()
        self.load_button = QPushButton("Load protein structures")
        self.select_structure = QComboBox()
        self.create_network = QPushButton("Create network")
        layout1.addWidget(self.load_button)
        layout1.addWidget(self.select_structure)
        layout1.addWidget(self.create_network)
        box1.setLayout(layout1)

        # Section 2: Visualization Options
        box2 = QGroupBox("Select Visualization:")
        layout2 = QVBoxLayout()
        self.view_in_protein = QPushButton("View network in protein")
        self.view_as_network = QPushButton("View as network")
        layout2.addWidget(self.view_in_protein)
        layout2.addWidget(self.view_as_network)
        box2.setLayout(layout2)

        # Section 3: Node Centrality
        box3 = QGroupBox("Select a Network Centrality:")
        layout3 = QVBoxLayout()
        self.node_betweenness = QPushButton("Node betweenness")
        self.node_closeness = QPushButton("Node closeness")
        self.node_degree = QPushButton("Node degree")
        self.node_eigenvector = QPushButton("Node eigenvector")
        layout3.addWidget(self.node_betweenness)
        layout3.addWidget(self.node_closeness)
        layout3.addWidget(self.node_degree)
        layout3.addWidget(self.node_eigenvector)

        # New Centrality Highlight Controls
        self.centrality_combo = QComboBox()
        self.centrality_combo.addItems(["Degree", "Closeness", "Betweenness", "Eigenvector"])
        self.centrality_topn = QSpinBox()
        self.centrality_topn.setRange(1, 100)
        self.centrality_topn.setValue(5)
        self.centrality_button = QPushButton("Highlight Top Nodes")
        layout3.addWidget(QLabel("Centrality Metric:"))
        layout3.addWidget(self.centrality_combo)
        layout3.addWidget(QLabel("Top N Nodes:"))
        layout3.addWidget(self.centrality_topn)
        layout3.addWidget(self.centrality_button)

        box3.setLayout(layout3)

        # Section 4: Edge Analysis
        box4 = QGroupBox("Edge Analysis:")
        layout4 = QVBoxLayout()
        self.edge_betweenness = QPushButton("Edge betweenness")
        self.edge_weight = QPushButton("Edge weight")
        self.edge_weight_ebc = QPushButton("Edge weight * EBC")
        self.reset_edges = QPushButton("Reset edges")
        layout4.addWidget(self.edge_betweenness)
        layout4.addWidget(self.edge_weight)
        layout4.addWidget(self.edge_weight_ebc)
        layout4.addWidget(self.reset_edges)
        box4.setLayout(layout4)

        # Section 5: Dashboard & Comparison
        box5 = QGroupBox("Guide for customizing the visualization:")
        layout5 = QVBoxLayout()
        self.generate_dashboard = QPushButton("Generate Dashboard")
        self.protein1 = QComboBox()
        self.protein2 = QComboBox()
        self.unique_edges1 = QPushButton("Unique edge protein 1")
        self.unique_edges2 = QPushButton("Unique edges protein 2")
        layout5.addWidget(self.generate_dashboard)
        layout5.addWidget(self.protein1)
        layout5.addWidget(self.protein2)
        layout5.addWidget(self.unique_edges1)
        layout5.addWidget(self.unique_edges2)
        box5.setLayout(layout5)

        # Layout grouping
        top_row = QHBoxLayout()
        top_row.addWidget(box1)
        top_row.addWidget(box2)

        middle_row = QHBoxLayout()
        middle_row.addWidget(box3)
        middle_row.addWidget(box4)

        main_layout.addLayout(top_row)
        main_layout.addLayout(middle_row)
        main_layout.addWidget(box5)
