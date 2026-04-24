# 🧬 In Silico Assessment of Aspirin Binding Dynamics within SARS-CoV-2 Main Protease (Mpro)

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Status](https://img.shields.io/badge/status-complete-brightgreen)
![Tool](https://img.shields.io/badge/docking-AutoDock%20Vina-orange)
![Visualization](https://img.shields.io/badge/visualization-UCSF%20ChimeraX-purple)

A comprehensive computational investigation into the molecular binding interaction between **SARS-CoV-2 Main Protease (Mpro)** and **Aspirin (acetylsalicylic acid)** using molecular docking simulation. This study evaluates the drug repurposing potential of a common anti-inflammatory compound against a validated antiviral target.

---

## 📋 Table of Contents

- [Background](#-background)
- [Objectives](#-objectives)
- [Tools & Dependencies](#️-tools--dependencies)
- [Repository Structure](#-repository-structure)
- [Methodology](#-methodology)
  - [1. Target Protein Preparation](#1-target-protein-preparation)
  - [2. Ligand Preparation](#2-ligand-preparation)
  - [3. Molecular Docking](#3-molecular-docking)
  - [4. Visualization & Analysis](#4-visualization--analysis)
- [Results](#-results)
- [Interpretation & Limitations](#-interpretation--limitations)
- [Future Directions](#-future-directions)
- [How to Reproduce](#-how-to-reproduce)
- [File Descriptions](#-file-descriptions)
- [References](#-references)

---

## 🔬 Background

The SARS-CoV-2 **Main Protease (Mpro)**, also known as the 3C-like protease (3CLpro), is an essential enzyme encoded by the viral genome. It is responsible for cleaving the viral polyproteins pp1a and pp1ab into individual non-structural proteins (nsps) that are required for viral replication and transcription. Without functional Mpro activity, the virus cannot replicate effectively, making Mpro one of the most well-characterized and druggable antiviral targets identified to date.

**Aspirin (acetylsalicylic acid)** is a widely used non-steroidal anti-inflammatory drug (NSAID) with a well-established safety profile. It exerts its pharmacological effects primarily through inhibition of cyclooxygenase (COX) enzymes. Given the interest in drug repurposing during the COVID-19 pandemic, evaluating aspirin's structural compatibility with Mpro offers a low-cost, computationally tractable starting point for investigating whether existing pharmaceuticals could serve as antiviral scaffolds.

---

## 🎯 Objectives

1. Retrieve and prepare the crystal structure of SARS-CoV-2 Mpro (PDB ID: **9O6F**) for computational docking.
2. Perform molecular docking of Aspirin against the Mpro active site using AutoDock Vina.
3. Evaluate binding affinity scores (ΔG, kcal/mol) and pose convergence across multiple docking runs.
4. Visualize the protein–ligand complex and analyze key structural interactions using UCSF ChimeraX.
5. Provide a reproducible computational pipeline as a baseline for further drug repurposing studies.

---

## 🛠️ Tools & Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| [AutoDock Vina](https://vina.scripps.edu/) | ≥ 1.2.0 | Molecular docking engine |
| [UCSF ChimeraX](https://www.cgl.ucsf.edu/chimerax/) | ≥ 1.6 | Protein–ligand visualization, cavity analysis |
| [Open Babel](http://openbabel.org/) | ≥ 3.1.1 | Ligand format conversion (SDF/MOL2 → PDBQT) |
| [MGLTools / AutoDockTools](https://ccsb.scripps.edu/mgltools/) | 1.5.7 | Protein PDBQT preparation |
| [Python](https://www.python.org/) | ≥ 3.8 | Scripting, automation |
| [Matplotlib](https://matplotlib.org/) / [Seaborn](https://seaborn.pydata.org/) | Latest | Energy distribution plots |
| [RCSB Protein Data Bank](https://www.rcsb.org/) | — | Source of crystal structure (PDB ID: 9O6F) |

---

## 📁 Repository Structure

```
aspirin-mpro-docking/
│
├── README.md                        # Project overview and documentation
│
├── structures/
│   ├── 9O6F.pdb                     # Raw crystal structure from PDB (S113A mutant)
│   ├── 9O6F_prepared.pdbqt          # Cleaned, PDBQT-formatted receptor
│   └── aspirin.pdbqt                # Ligand in PDBQT format (AutoDock-ready)
│
├── ligand/
│   ├── aspirin.sdf                  # Aspirin 3D structure (downloaded/generated)
│   ├── aspirin.mol2                 # MOL2 format (intermediate)
│   └── aspirin_minimized.pdbqt      # Energy-minimized PDBQT ligand
│
├── docking/
│   ├── config.txt                   # AutoDock Vina configuration file (grid box params)
│   ├── docking_results.pdbqt        # All docking poses output
│   ├── docking_log.txt              # Full docking log with affinity scores
│   └── best_pose.pdbqt              # Highest-affinity docking pose
│
├── analysis/
│   ├── energy_distribution.png      # Binding affinity distribution across poses
│   ├── energy_distribution.py       # Script to generate affinity plots
│   ├── cavity_visualization.cxs     # UCSF ChimeraX session (cavity + ligand)
│   └── interaction_summary.txt      # Summary of key binding site contacts
│
├── figures/
│   ├── binding_pose_overview.png    # Full view of Aspirin in Mpro binding pocket
│   ├── active_site_closeup.png      # Close-up of ligand–residue interactions
│   └── cavity_surface.png           # Surface representation of binding cavity
│
└── scripts/
    ├── prepare_receptor.py          # Automates receptor preparation
    ├── prepare_ligand.py            # Automates ligand preparation
    └── run_docking.sh               # Shell script to execute AutoDock Vina
```

---

## 🔄 Methodology

### 1. Target Protein Preparation

The crystal structure of SARS-CoV-2 Mpro used in this study is **PDB ID: 9O6F**, which represents the **S113A mutant** complexed with an inhibitor. This mutant retains the overall architecture of the catalytic site and is suitable for docking-based investigations.

**Steps:**
1. The structure was downloaded from the [RCSB Protein Data Bank](https://www.rcsb.org/structure/9O6F).
2. Co-crystallized ligands, water molecules, and heteroatoms were removed to expose the apo binding site.
3. Hydrogen atoms were added to the protein at physiological pH (~7.4) using AutoDockTools or ChimeraX's `AddH` function.
4. Gasteiger charges were computed and non-polar hydrogens were merged.
5. The prepared structure was saved in **PDBQT format** for use with AutoDock Vina.

```bash
# Example using MGLTools Python scripts
python prepare_receptor4.py -r 9O6F.pdb -o 9O6F_prepared.pdbqt -A checkhydrogens
```

### 2. Ligand Preparation

Aspirin (acetylsalicylic acid) — **PubChem CID: 2244** — was obtained as a 3D SDF structure from PubChem or generated using a chemical structure tool.

**Steps:**
1. The 3D structure was downloaded from [PubChem](https://pubchem.ncbi.nlm.nih.gov/compound/2244).
2. Open Babel was used to convert the SDF file to MOL2, and subsequently to PDBQT format.
3. Partial charges (Gasteiger) were assigned and non-polar hydrogens merged during conversion.

```bash
# Convert SDF to PDBQT using Open Babel
obabel aspirin.sdf -O aspirin.pdbqt --gen3d -p 7.4
```

### 3. Molecular Docking

AutoDock Vina was used to perform rigid-receptor, flexible-ligand docking.

**Grid Box Configuration:**

The grid box was defined around the catalytic dyad residues (His41 and Cys145) of the Mpro active site — the canonical substrate-binding pocket.

```ini
# config.txt
receptor = 9O6F_prepared.pdbqt
ligand   = aspirin.pdbqt

center_x = [x-coordinate of active site centroid]
center_y = [y-coordinate of active site centroid]
center_z = [z-coordinate of active site centroid]

size_x = 20
size_y = 20
size_z = 20

exhaustiveness = 8
num_modes       = 9
energy_range    = 3
```

**Execution:**

```bash
# Run AutoDock Vina
vina --config config.txt --out docking_results.pdbqt --log docking_log.txt
```

### 4. Visualization & Analysis

Results were visualized and analyzed using **UCSF ChimeraX**:

- The best docking pose was loaded alongside the prepared receptor.
- **Internal cavity detection** was used to confirm that Aspirin occupies the expected binding pocket.
- **Molecular surface** representations were generated to assess steric complementarity.
- Binding affinity values across all poses were extracted from the docking log and plotted as a distribution using Python (Matplotlib/Seaborn).

---

## 📊 Results

### Binding Affinity Summary

| Pose | Binding Affinity (ΔG, kcal/mol) |
|------|----------------------------------|
| **1 (Best)** | **-4.916** |
| 2 | -4.801 |
| 3 | -4.778 |
| 4 | -4.710 |
| 5 | -4.689 |
| 6 | -4.652 |
| 7 | -4.598 |
| 8 | -4.554 |
| 9 | -4.501 |

> **Best binding affinity: -4.916 kcal/mol** (Pose 1)

### Key Observations

- **Convergence:** The narrow energy range across all 9 poses (~0.4 kcal/mol spread) indicates good convergence of the docking simulation and suggests that Aspirin consistently localizes to the same binding region.
- **Binding pocket occupancy:** Cavity analysis in ChimeraX confirmed that Aspirin occupies the catalytic pocket of Mpro, positioned within the region defined by subsites S1 and S2.
- **Interaction profile:** The ligand orientation suggests potential non-covalent interactions including **hydrogen bonding** with polar active-site residues and **hydrophobic contacts** within the hydrophobic sub-pocket (S2 site).
- **Affinity interpretation:** A ΔG of -4.916 kcal/mol reflects a **weak-to-moderate** binding affinity, broadly consistent with small, relatively polar molecules like Aspirin interacting with a large, substrate-optimized binding site.

---

## ⚠️ Interpretation & Limitations

While the docking results provide useful preliminary insights, several important caveats apply:

| Limitation | Description |
|------------|-------------|
| **Rigid receptor** | AutoDock Vina treats the receptor as rigid; protein flexibility and induced-fit effects are not modeled |
| **No solvent effects** | Implicit solvation is not applied; explicit water molecules and desolvation penalties are not accounted for |
| **No entropy estimation** | Docking scores approximate enthalpic contributions only; conformational entropy is not fully captured |
| **Static snapshot** | Docking uses a single crystal structure; dynamic conformational changes are ignored |
| **Covalent interactions** | Aspirin can act as an irreversible acetylating agent; covalent docking was not performed in this study |
| **Mutant structure** | The S113A mutant (9O6F) differs from the wild-type at position 113; functional implications of this substitution should be considered |

These results should be interpreted as **hypothesis-generating** and require validation through more rigorous computational and experimental methods.

---

## 🚀 Future Directions

This study establishes a reproducible baseline that can be extended in several meaningful directions:

1. **Molecular Dynamics (MD) Simulations** — Use GROMACS or AMBER to simulate the Aspirin–Mpro complex over nanosecond timescales and assess binding stability and protein flexibility.
2. **MM-PBSA / MM-GBSA Binding Free Energy** — Compute more accurate binding free energies using end-point thermodynamic methods post-MD simulation.
3. **Covalent Docking** — Given Aspirin's acetylating mechanism, evaluate covalent docking modes targeting Cys145 using CovDock or AutoDock with covalent options.
4. **Comparative Benchmarking** — Compare Aspirin's binding profile against known high-affinity Mpro inhibitors (e.g., Nirmatrelvir/PF-07321332, N3) to contextualize the observed affinity.
5. **ADMET Profiling** — Assess Aspirin's drug-likeness, ADMET properties, and predicted toxicity in the context of antiviral repurposing using tools like SwissADME or pkCSM.
6. **Ensemble Docking** — Perform docking against multiple Mpro conformations derived from MD trajectories to capture receptor flexibility.
7. **In Vitro Validation** — Design cell-based assays or enzyme inhibition assays to validate computational predictions experimentally.

---

## ▶️ How to Reproduce

### Prerequisites

Ensure the following tools are installed:
- AutoDock Vina (`vina` command accessible in PATH)
- Open Babel (`obabel`)
- MGLTools (for receptor preparation scripts)
- UCSF ChimeraX (for visualization)
- Python ≥ 3.8 with `matplotlib`, `seaborn`, `pandas`

### Step-by-Step

```bash
# 1. Clone the repository
git clone https://github.com/<your-username>/aspirin-mpro-docking.git
cd aspirin-mpro-docking

# 2. Prepare the receptor
python scripts/prepare_receptor.py

# 3. Prepare the ligand
python scripts/prepare_ligand.py

# 4. Run docking
bash scripts/run_docking.sh

# 5. Generate energy distribution plot
python analysis/energy_distribution.py

# 6. Open ChimeraX session for visualization
# Launch ChimeraX and open: analysis/cavity_visualization.cxs
```

### Viewing Results in ChimeraX

1. Open ChimeraX.
2. Load the session file: `File → Open → analysis/cavity_visualization.cxs`
3. Alternatively, manually open:
   - Receptor: `structures/9O6F_prepared.pdbqt`
   - Best pose: `docking/best_pose.pdbqt`
4. Use `Surface → Show Surface` and `Tools → Structure Analysis → Measure and Color Blobs` for cavity visualization.

---

## 📄 File Descriptions

| File | Description |
|------|-------------|
| `structures/9O6F.pdb` | Raw PDB file of Mpro S113A mutant (as downloaded from RCSB) |
| `structures/9O6F_prepared.pdbqt` | Cleaned and PDBQT-formatted receptor ready for docking |
| `ligand/aspirin.pdbqt` | Aspirin in PDBQT format with Gasteiger charges |
| `docking/config.txt` | AutoDock Vina configuration with grid box parameters |
| `docking/docking_results.pdbqt` | All 9 docking poses in PDBQT format |
| `docking/docking_log.txt` | Full docking output log with per-pose affinity scores |
| `docking/best_pose.pdbqt` | Top-ranked binding pose (ΔG = -4.916 kcal/mol) |
| `analysis/energy_distribution.py` | Python script for generating affinity score plots |
| `analysis/energy_distribution.png` | Binding affinity distribution across poses |
| `analysis/cavity_visualization.cxs` | UCSF ChimeraX session file for interactive exploration |
| `analysis/interaction_summary.txt` | Summary of predicted binding site contacts and residue proximity |
| `figures/binding_pose_overview.png` | Full-view rendering of the docked complex |
| `figures/active_site_closeup.png` | Close-up of Aspirin in the catalytic pocket |

---

## 📚 References

1. **Jin, Z. et al.** (2020). Structure of Mpro from SARS-CoV-2 and discovery of its inhibitors. *Nature*, 582, 289–293. https://doi.org/10.1038/s41586-020-2223-y
2. **Trott, O. & Olson, A.J.** (2010). AutoDock Vina: Improving the speed and accuracy of docking. *Journal of Computational Chemistry*, 31(2), 455–461. https://doi.org/10.1002/jcc.21334
3. **Pettersen, E.F. et al.** (2021). UCSF ChimeraX: Structure visualization for researchers, educators, and developers. *Protein Science*, 30(1), 70–82. https://doi.org/10.1002/pro.3943
4. **RCSB Protein Data Bank** — PDB ID: 9O6F. https://www.rcsb.org/structure/9O6F
5. **PubChem** — Aspirin (CID: 2244). https://pubchem.ncbi.nlm.nih.gov/compound/2244
6. **O'Boyle, N.M. et al.** (2011). Open Babel: An open chemical toolbox. *Journal of Cheminformatics*, 3, 33. https://doi.org/10.1186/1758-2946-3-33

---

## 📜 License

This project is licensed under the [MIT License](LICENSE). You are free to use, adapt, and distribute this work with appropriate attribution.

---

## 🙏 Acknowledgements

- The RCSB Protein Data Bank for open access to macromolecular structures.
- The Scripps Research Institute for the development and maintenance of AutoDock Vina.
- The UCSF ChimeraX team for providing a powerful and accessible visualization platform.

---

*For questions, issues, or contributions, please open an issue or pull request in the repository.*
