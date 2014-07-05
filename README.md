RDMMP
=====

Drug Markov Mean Properties in R

It is a tool for the calculation of Markov Mean Properties (MMPs) molecular descriptors for drugs using SMILES formulas and atom physical - chemical properties as inputs. It is using the MInD-Prot tool formulas for drugs.
The atom properties are Mulliken electronegativity (EM), Kang-Jhon polarizability (PKJ), van der Waals area (vdWA) and atom contribution to partition coefficient (AC2P). The file with these properties can be completed with any others or replaced with different ones.

Script steps: 
- Read the inputs: SMILES formulas and atom properties
- Get connectivity matrix (CM), nodes = atoms, edges = chemical bonds
- Get weights vector (w) for each atom property
- Calculate weighted matrix (W) using CM and w
- Calculate transition probability (P) based on W
- Calculate absolute initial probability vectpr (p0j) based on W
- Calculate k powers of P -> Pk matrices
- Calculate matrix products p0j * Pk * w -> MMPk (descriptors for each k and atomic property)
- Calculate MMP for each atom property by averaging MMPk values

The next steps:
- create the same MMPs for each type of atom (Csat, Cinst, Halog, Hetero)
- correct errors
- create a console version
