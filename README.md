# AFAR-WQS (Assimilation Factor Analysis of Rivers for Water Quality Simulation)
AFAR-WQS (Assimilation Factor Analysis of Rivers for Water Quality Simulation), a MATLAB-based tool designed to simulate key water quality determinants using the concept of steady-state assimilation factors. The tool combines graph theory and a recursive algorithm to model assimilation factors, concentrations, and loads for each water quality determinant, accounting for both diffuse and point-source pollution. AFAR-WQS supports large-scale simulations with minimal computational demand and low data requirements, offering researchers, managers, and policymakers an accessible tool for optimizing intervention strategies and resource allocation in river basin management.

## How does AFAR-WQS operate?
AFAR-WQS operates based on a graph representation of the drainage network within a basin. Each river reach is modeled as a pair of nodes connected by a reach, as illustrated in Figure 1-a. The connectivity of the topological network is resolved using a recursive al-gorithm, which begins at the basin’s outlet reach. From this starting point, the algorithm traverses the drainage network upstream until it encounters a headwater reach, defined as one with no incoming connections. Upon reaching such a headwater, the model starts accumulating the loads for various water quality determinants, estimating both assimilation factors and concentrations for that specific reach. The algo-rithm then proceeds downstream to the immediately adjacent reach. If this downstream reach has incoming connections, the model recursively moves upstream, repeating the process described, until it returns to the original confluence point. 
This modeling scheme assumes that the physicochemical conditions are uniform throughout the entire reach, as it represents the smallest modeling unit in AFAR-WQS. Therefore, the assimilation factors are calculated based on the average values of the phys-icochemical characteristics of the reach. Moreover, the recursive approach used by the tool to solve the topo-logical network allows for the analysis of networks with a large number of nodes while using minimal storage resources, something that would not be possible with adjacency matrix-based schemes, which require significantly higher computational resources. 

# Requirements
AFAR-WQS toolbox requires MATLAB™ R2023b or newer version, with the Mapping Toolbox™.

# Documentation
The AFAR-WQS toolbox comes with an user guide for the implementation of the toolbox, fully commented code and interactive functions for first time users. Refer to the user guide for more informations.

# References
- Nogales, J. and Rogéliz-Prada. (2024) AFAR-WQS: A quick and simple toolbox for water quality simulation. Water
- Nogales, J.; Rogéliz-Prada, C.; Cañon, M.A.; Vargas-Luna, A. (2023). An Integrated Methodological Framework for the Durable Conservation of Freshwater Ecosystems: A Case Study in Colombia’s Caquetá River Basin. Front Environ Sci 2023, 11, 1264392, doi:10.3389/fenvs.2023.1264392.
