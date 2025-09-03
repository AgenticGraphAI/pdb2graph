from Bio.PDB import PDBParser, is_aa
from scipy.spatial import cKDTree
import pandas as pd
import numpy as np
from collections import defaultdict
import networkx as nx
import os
import tempfile

# PyMOL
import pymol
from pymol import cmd


class ProteinAnalyzerLogic:
    def __init__(self, ui):
        self.ui = ui
        # you likely set this elsewhere when a file is chosen; keeping attribute here for safety
        if not hasattr(self, "structure_path"):
            self.structure_path = ""  # set by your loader logic
        self.connect_signals()

    def connect_signals(self):
        """Wire only the buttons that exist in the current UI to avoid AttributeErrors."""
        def connect_safe(attr, slot):
            if hasattr(self.ui, attr):
                getattr(self.ui, attr).clicked.connect(slot)

        connect_safe("create_network", self.build_psn_from_loaded_structure)
        connect_safe("view_in_protein", self.visualize_psn_in_protein)
        connect_safe("view_as_network", self.show_as_network)
        connect_safe("edge_weight", self.visualize_weighted_psn_in_protein)

        # Node centrality buttons
        connect_safe("node_degree", lambda: self.apply_centrality("degree", n=self._ui_topn()))
        connect_safe("node_betweenness", lambda: self.apply_centrality("betweenness", n=self._ui_topn()))
        connect_safe("node_closeness", lambda: self.apply_centrality("closeness", n=self._ui_topn()))
        connect_safe("node_eigenvector", lambda: self.apply_centrality("eigenvector", n=self._ui_topn()))

        # Combo + button (optional)
        if hasattr(self.ui, "centrality_button") and hasattr(self.ui, "centrality_combo"):
            self.ui.centrality_button.clicked.connect(
                lambda: self.apply_centrality(self.ui.centrality_combo.currentText().lower(), n=self._ui_topn())
            )

        # Edge betweenness overlays
        connect_safe("edge_betweenness", lambda: self.visualize_topN_ebc_edges(n=self._ui_topn()))
        connect_safe("edge_weight_ebc", self.visualize_edge_weight_times_ebc)
        connect_safe("reset_edges", lambda: (self._clear_prev_edge_overlays("EBC"), self._clear_prev_edge_overlays("EBCW")))

        # When user clicks "Load protein structures", try to populate the combo from current PyMOL objects
        if hasattr(self.ui, "load_button") and hasattr(self.ui, "select_structure"):
            def _refresh_combo_from_pymol():
                try:
                    # If they already loaded objects into PyMOL, offer a quick temp-export option entry
                    objs = cmd.get_names("objects")
                    if objs:
                        self.ui.select_structure.clear()
                        # NOTE: we can’t know real file paths for PyMOL-loaded objects;
                        # they will be temp-exported at build time by _ensure_structure_path().
                        # Put object names in the combo just for user feedback.
                        for o in objs:
                            self.ui.select_structure.addItem(f"[PyMOL object] {o}")
                except Exception:
                    pass
            self.ui.load_button.clicked.connect(_refresh_combo_from_pymol)

        # If select_structure holds a real filesystem path, update structure_path on change
        if hasattr(self.ui, "select_structure"):
            try:
                self.ui.select_structure.currentTextChanged.connect(
                    lambda txt: setattr(self, "structure_path", txt) if os.path.exists(str(txt)) else None
                )
            except Exception:
                pass


    # ---------- Core PSN build ----------
    def generate_weighted_and_unweighted_psn(self, pdb_path, distance_cutoff=4.5):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", pdb_path)

        atom_coords = []
        res_ids = []

        for model in structure:
            for chain in model:
                if chain.id != "A":
                    continue
                for residue in chain:
                    if not is_aa(residue):
                        continue
                    res_id = (chain.id, residue.id[1], residue.get_resname())
                    for atom in residue:
                        if atom.element != "H":
                            atom_coords.append(atom.coord)
                            res_ids.append(res_id)

        atom_coords = np.array(atom_coords)
        if len(atom_coords) == 0:
            print("[PDB2Graph] No atoms collected; check structure/chain filter.")
            return

        kdtree = cKDTree(atom_coords)
        pairs = kdtree.query_pairs(r=distance_cutoff)

        contacts = defaultdict(set)
        for i, j in pairs:
            res1 = res_ids[i]
            res2 = res_ids[j]
            if res1[1] != res2[1]:  # ignore same-residue atom pairs
                contacts[(res1[0], res1[1])].add((res2[0], res2[1]))
                contacts[(res2[0], res2[1])].add((res1[0], res1[1]))

        weighted_edges = []
        unweighted_edges = []

        for (chain_res1), neighbors in contacts.items():
            for chain_res2 in neighbors:
                res1 = (chain_res1[0], chain_res1[1], "")
                res2 = (chain_res2[0], chain_res2[1], "")
                unweighted_edges.append((res1, res2))

        weight_counts = defaultdict(int)
        for (res1, res2) in unweighted_edges:
            key = tuple(sorted((res1, res2)))
            weight_counts[key] += 1

        for key, w in weight_counts.items():
            weighted_edges.append((key[0], key[1], w))

        self.unweighted_df = pd.DataFrame(unweighted_edges, columns=["Residue1", "Residue2"])
        self.weighted_df = pd.DataFrame(weighted_edges, columns=["Residue1", "Residue2", "Weight"])
        print("[PDB2Graph] PSN generated (weighted & unweighted).")

    def build_psn_from_loaded_structure(self):
        try:
            pdb_path = self._ensure_structure_path()
        except RuntimeError as e:
            print(f"[PDB2Graph] Error: {e}")
            return

        try:
            self.generate_weighted_and_unweighted_psn(pdb_path)
            print("[PDB2Graph] PSN successfully built from loaded structure.")
        except Exception as e:
            print(f"[PDB2Graph] Error building PSN: {e}")


        # ---------- Helpers to resolve structure path ----------
    def _ensure_structure_path(self):
        """
        Return a usable PDB path, trying (in order):
        1) self.structure_path (if set and exists)
        2) current text in UI select_structure (if looks like an existing file)
        3) any loaded PyMOL object exported to a temp .pdb
        Raises RuntimeError if nothing is available.
        """
        # 1) explicit path set
        path = getattr(self, "structure_path", "") or ""
        if isinstance(path, str) and path and os.path.exists(path):
            return path

        # 2) UI combobox path
        if hasattr(self.ui, "select_structure"):
            try:
                candidate = self.ui.select_structure.currentText().strip()
                if candidate and os.path.exists(candidate):
                    self.structure_path = candidate
                    return candidate
            except Exception:
                pass

        # 3) Fall back to any loaded object in PyMOL: save a temp PDB
        try:
            objs = cmd.get_names("objects")
            if objs:
                obj = objs[0]  # first object
                fd, tmp_path = tempfile.mkstemp(prefix="pdb2graph_", suffix=".pdb")
                os.close(fd)
                cmd.save(tmp_path, obj)  # export current object to PDB
                self.structure_path = tmp_path
                print(f"[PDB2Graph] Exported '{obj}' to temp PDB: {tmp_path}")
                return tmp_path
        except Exception as e:
            print(f"[PDB2Graph] Could not export loaded object to PDB: {e}")

        raise RuntimeError("No structure source found. Load a PDB in PyMOL or pick a path in the UI.")
    
    # ---------- Visualizations ----------
    def visualize_psn_in_protein(self):
        """Show protein cartoon + CA spheres + PSN edges overlaid as dashed distances."""
        if not hasattr(self, "unweighted_df"):
            print("[PDB2Graph] Error: PSN not yet built.")
            return

        try:
            # Base protein view
            cmd.hide("everything", "all")
            cmd.show("cartoon", "all")
            cmd.color("lightorange", "chain A")
            cmd.show("spheres", "chain A and name CA")
            cmd.set("sphere_scale", 0.4, "chain A and name CA")
            cmd.color("yellow", "chain A and name CA")

            # Clear old overlays and redraw this PSN
            self._clear_overlay_objects(prefix_list=("EBC", "EBCW", "NET"))
            idx = 0
            for _, row in self.unweighted_df.iterrows():
                res1 = row["Residue1"]
                res2 = row["Residue2"]
                sel1 = self._edge_sel_from_tuple(res1)
                sel2 = self._edge_sel_from_tuple(res2)
                idx += 1
                self._draw_edge_distance("NET", idx, sel1, sel2, dash_color="red", dash_width=1.2)

            cmd.bg_color("white")
            print("[PDB2Graph] PSN overlay on protein complete.")
        except Exception as e:
            print(f"[PDB2Graph] Error during embedded PSN visualization: {e}")

    def show_as_network(self):
        """
        Network-only view inside PyMOL:
        - Hide everything (cartoon, sticks, ribbons, labels)
        - Show ONLY PSN nodes (Cα spheres)
        - Draw ALL PSN edges into ONE measurement object (PSN_EDGES)
        """
        if not hasattr(self, "unweighted_df"):
            print("[PDB2Graph] Error: PSN not yet built.")
            return
        try:
            cmd.set("defer_builds_mode", 3)
            # Clean old spam first (fixes the “1917 objects” problem)
            self._nuke_psn_objects()

            # Hide everything, then show PSN nodes
            cmd.hide("everything", "all")
            sel_nodes = self._psn_nodes_selection()
            if sel_nodes != "none":
                cmd.show("spheres", sel_nodes)
                cmd.set("sphere_scale", 0.5, sel_nodes)
                cmd.color("yellow", sel_nodes)

            # Draw all edges into ONE measurement object name
            # Re-using the same name appends segments to that single object.
            for res1, res2 in self.unweighted_df.itertuples(index=False):
                sel1 = f"chain {res1[0]} and resi {res1[1]} and name CA"
                sel2 = f"chain {res2[0]} and resi {res2[1]} and name CA"
                cmd.distance("PSN_EDGES", sel1, sel2)

            # Style the single measurement object like sticks
            cmd.hide("labels", "PSN_EDGES")
            cmd.set("dash_gap", 0.0, "PSN_EDGES")
            cmd.set("dash_color", "gray", "PSN_EDGES")
            cmd.set("dash_width", 1.4, "PSN_EDGES")

            cmd.bg_color("white")
            print("[PDB2Graph] Network-only: Cα spheres + single-object edges drawn.")
        except Exception as e:
            print(f"[PDB2Graph] Error (network-only): {e}")
        finally:
            try:
                cmd.set("defer_builds_mode", 0)
            except Exception:
                pass

    def _nuke_psn_objects(self):
        """Delete prior overlay spam so the object list doesn’t explode."""
        prefixes = ("NET", "EBC", "EBCW")
        kill = {"PSN_EDGES", "PSN_NET"}  # our single-edge object names
        try:
            for obj in cmd.get_names("objects"):
                if obj in kill or any(obj.startswith(p + "_") for p in prefixes):
                    cmd.delete(obj)
        except Exception as e:
            print(f"[PDB2Graph] Cleanup warning: {e}")

    def _psn_nodes_selection(self):
        """Make one selection for PSN nodes only (Cα)."""
        nodes_by_chain = {}
        for res1, res2 in self.unweighted_df.itertuples(index=False):
            nodes_by_chain.setdefault(res1[0], set()).add(int(res1[1]))
            nodes_by_chain.setdefault(res2[0], set()).add(int(res2[1]))
        if not nodes_by_chain:
            return "none"
        parts = []
        for ch, resis in nodes_by_chain.items():
            resi_str = "+".join(str(r) for r in sorted(resis))
            parts.append(f"(chain {ch} and name CA and resi {resi_str})")
        sel = " or ".join(parts)
        cmd.select("PSN_NODES", sel)
        return "PSN_NODES"

    def visualize_weighted_psn_in_protein(self):
        if not hasattr(self, "weighted_df"):
            print("[PDB2Graph] Error: Weighted PSN not yet built.")
            return

        max_weight = float(self.weighted_df["Weight"].max() or 1.0)
        cmd.hide("everything", "all")
        cmd.show("cartoon", "all")
        cmd.color("lightorange", "chain A")

        for _, row in self.weighted_df.iterrows():
            res1 = row["Residue1"]
            res2 = row["Residue2"]
            weight = float(row["Weight"])
            selection1 = self._edge_sel_from_tuple(res1)
            selection2 = self._edge_sel_from_tuple(res2)
            cmd.bond(selection1, selection2)
            scale = 0.1 + 0.4 * (weight / max_weight)
            cmd.set_bond("stick_radius", scale, selection1, selection2)
            cmd.set_bond("stick_color", "blue", selection1, selection2)
            cmd.show("sticks", f"{selection1}, {selection2}")

        cmd.bg_color("white")
        print("[PDB2Graph] Weighted PSN visualization complete.")

    # ---------- Node centrality (spheres) ----------
    def apply_centrality(self, metric="degree", n=10):
        if not hasattr(self, "unweighted_df"):
            print("[PDB2Graph] Error: PSN not built yet.")
            return

        G = nx.Graph()
        for _, row in self.unweighted_df.iterrows():
            G.add_edge(row["Residue1"], row["Residue2"])

        if metric == "degree":
            centrality = nx.degree_centrality(G)
        elif metric == "betweenness":
            centrality = nx.betweenness_centrality(G, normalized=True)
        elif metric == "closeness":
            centrality = nx.closeness_centrality(G)
        elif metric == "eigenvector":
            # numpy method is stable for undirected graphs
            centrality = nx.eigenvector_centrality_numpy(G)
        else:
            print(f"[PDB2Graph] Unknown centrality metric: {metric}")
            return

        top_n = sorted(centrality.items(), key=lambda x: x[1], reverse=True)[:int(n)]
        if not top_n:
            print("[PDB2Graph] No central nodes found.")
            return

        # show spheres for top nodes
        cmd.show("spheres")
        best = top_n[0][1] or 1e-6
        for (chain, resi, _), score in top_n:
            sel = f"chain {chain} and resi {resi} and name CA"
            scale = 0.3 + 0.9 * (score / best)
            cmd.set("sphere_scale", float(scale), sel)
            cmd.color("yellow", sel)

        print(f"[PDB2Graph] Centrality visualization complete: {metric}, top {len(top_n)}")

    # ---------- Edge betweenness helpers ----------
    def _compute_edge_betweenness(self, G: nx.Graph):
        """Compute edge betweenness centrality; return dict with normalized keys (u,v) sorted."""
        if G.number_of_edges() == 0:
            return {}
        ebc = nx.edge_betweenness_centrality(
            G,
            weight="Weight" if any("Weight" in d for *_, d in G.edges(data=True)) else None,
            normalized=True,
        )
        return {tuple(sorted((u, v))): val for (u, v), val in ebc.items()}

    def _normalize_scores(self, d):
        if not d:
            return {}
        m = max(d.values())
        if m == 0:
            return {k: 0.0 for k in d}
        return {k: v / m for k, v in d.items()}

    def _clear_prev_edge_overlays(self, prefix="EBC"):
        try:
            for obj in cmd.get_names("objects"):
                if obj.startswith(prefix + "_"):
                    cmd.delete(obj)
        except Exception:
            pass

    def _edge_sel_from_tuple(self, res_tuple):
        """res_tuple is expected like (chain, resi, resname)"""
        chain, resi = res_tuple[0], res_tuple[1]
        return f"chain {chain} and resi {resi} and name CA"

    # ---------- Edge betweenness overlays ----------
    def visualize_topN_ebc_edges(self, n=20):
        """Draw top‑N edges by edge betweenness as distance objects with width scaled by rank."""
        if not hasattr(self, "unweighted_df"):
            print("[PDB2Graph] Error: PSN not yet built.")
            return

        G = nx.Graph()
        for _, row in self.unweighted_df.iterrows():
            u, v = row["Residue1"], row["Residue2"]
            G.add_edge(u, v)

        ebc = self._compute_edge_betweenness(G)
        ranked = sorted(ebc.items(), key=lambda kv: kv[1], reverse=True)[:max(0, int(n))]
        if not ranked:
            print("[PDB2Graph] No edges selected.")
            return

        self._clear_prev_edge_overlays("EBC")
        cmd.hide("labels", "all")
        cmd.set("dash_gap", 0.0)
        cmd.set("dash_color", "red")

        min_w, max_w = 0.6, 2.6
        total = len(ranked)
        for i, ((u, v), _score) in enumerate(ranked, start=1):
            sel1 = self._edge_sel_from_tuple(u)
            sel2 = self._edge_sel_from_tuple(v)
            obj = f"EBC_{i}"
            try:
                cmd.distance(obj, sel1, sel2)
                # rank-based width
                width = min_w + (max_w - min_w) * (1 - (i - 1) / max(1, total - 1))
                cmd.set("dash_width", float(width), obj)
                cmd.hide("labels", obj)
            except Exception as e:
                print(f"[PDB2Graph] distance creation failed for {u}-{v}: {e}")

        print(f"[PDB2Graph] Drew top {total} edges by EBC.")

    def visualize_edge_weight_times_ebc(self):
        """Overlay weighted edges with width ∝ normalized(weight) × normalized(EBC)."""
        if not hasattr(self, "weighted_df") or not hasattr(self, "unweighted_df"):
            print("[PDB2Graph] Error: Need weighted and unweighted PSN built.")
            return

        G = nx.Graph()
        for _, row in self.unweighted_df.iterrows():
            G.add_edge(row["Residue1"], row["Residue2"])
        ebc = self._normalize_scores(self._compute_edge_betweenness(G))
        if not ebc:
            print("[PDB2Graph] EBC empty.")
            return

        self._clear_prev_edge_overlays("EBCW")
        cmd.hide("labels", "all")
        cmd.set("dash_gap", 0.0)
        cmd.set("dash_color", "blue")

        max_wt = float(self.weighted_df["Weight"].max() or 1.0)
        min_w, max_w = 0.2, 2.8
        i = 0

        for _, row in self.weighted_df.iterrows():
            u, v = row["Residue1"], row["Residue2"]
            key = tuple(sorted((u, v)))
            weight = float(row.get("Weight", 1.0))
            combo = (weight / max_wt) * ebc.get(key, 0.0)
            if combo <= 0:
                continue

            sel1 = self._edge_sel_from_tuple(u)
            sel2 = self._edge_sel_from_tuple(v)
            i += 1
            obj = f"EBCW_{i}"
            try:
                cmd.distance(obj, sel1, sel2)
                width = min_w + (max_w - min_w) * combo
                cmd.set("dash_width", float(width), obj)
                cmd.hide("labels", obj)
            except Exception as e:
                print(f"[PDB2Graph] distance creation failed for {u}-{v}: {e}")

        print(f"[PDB2Graph] Drew {i} edges with weight*EBC overlay.")


    # ---------- PSN render helpers (in-PyMOL) ----------
    def _edge_sel_from_tuple(self, res_tuple):
        """res_tuple is (chain, resi, resname)"""
        chain, resi = res_tuple[0], res_tuple[1]
        return f"chain {chain} and resi {resi} and name CA"

    def _clear_overlay_objects(self, prefix_list=("EBC", "EBCW", "NET")):
        """Remove previously drawn distance objects for PSN overlays."""
        try:
            names = cmd.get_names("objects")
            for name in names:
                if any(name.startswith(p + "_") for p in prefix_list):
                    cmd.delete(name)
        except Exception:
            pass

    def _draw_edge_distance(self, obj_prefix, idx, sel1, sel2, dash_color="gray", dash_width=1.2):
        name = f"{obj_prefix}_{idx}"
        cmd.distance(name, sel1, sel2)
        cmd.hide("labels", name)
        cmd.set("dash_gap", 0.0, name)
        cmd.set("dash_color", dash_color, name)
        cmd.set("dash_width", float(dash_width), name)
        return name


    # ---------- Small UI helper ----------
    def _ui_topn(self, default_val=20):
        if hasattr(self.ui, "centrality_topn"):
            try:
                return int(self.ui.centrality_topn.value())
            except Exception:
                return default_val
        return default_val
