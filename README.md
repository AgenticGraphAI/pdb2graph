# PDB2Graph: PyMOL Plugin for Protein Structure Networks

This repository contains a PyMOL plugin for constructing, analyzing, and visualizing protein structure networks (PSNs) directly within a structural biology workflow.

## Overview

The plugin enables:

- Construction of protein structure networks from PDB structures
- Visualization of residue-level interaction graphs in PyMOL
- Calculation of network centrality metrics (degree, betweenness, closeness, eigenvector)
- Exploration of structure–network relationships in ligand-bound proteins

## Background

Earlier versions of the plugin relied on the RING server for PSN generation. Changes to the RING API resulted in unreliable or empty outputs, motivating the development of an internal PSN generation pipeline.

## Current Direction

The project is transitioning toward a fully self-contained PSN engine, including:

- Distance-based and contact-based graph construction (e.g., Cα–Cα thresholds)
- Edge weighting based on geometric or interaction features
- Integration with downstream analysis (centrality, graph traversal, ML pipelines)

## Status

The plugin is under active redevelopment following the transition away from external dependencies. Previous versions were demonstrated in research settings (e.g., Blind Lab), and current work focuses on rebuilding the PSN generation layer for robustness and reproducibility.

## Long-Term Goal

To provide an interactive structural analysis tool that bridges:

Structure → Graph → Network Analysis → Machine Learning

within a single PyMOL-based interface.
