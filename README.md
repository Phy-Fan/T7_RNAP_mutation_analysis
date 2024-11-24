# Interaction Analysis Scripts
This repository contains three independent Python scripts for analyzing protein interactions: hydrogen bonds, salt bridges, and hydrophobic interactions. Each script is self-contained and can be executed separately.


## Setup Instructions
### Step 1: Create a Python Environment
Before running the scripts, create a Conda environment with Python 3.7:

```conda create -n interaction_env python=3.7 -y```

```conda activate interaction_env```


### Step 2: Install Required Packages
Once the virtual environment is active, install the necessary Python packages using pip by running the following command:

```pip install pymol-open-source mdtraj pandas numpy```

### Step 3: Clone the Repository
Download this repository to your local system. To clone the repository, run:

```git clone https://github.com/Phy-Fan/T7_RNAP_mutation_analysis```

## Script Descriptions
**Hydrogen Bond Analysis**

**File:** ```hydrogen_bond.py```

**Description:** Hydrogen bonds are calculated using MDTraj, defined by a donor-acceptor distance of less than 3.5 Å and a donor-hydrogen-acceptor angle greater than 120°.

**Salt Bridge Analysis**

**File:** ```salt_bridge.py```

**Description:** Identifies salt bridges by calculating distances between charged residues. A salt bridge is defined as a distance of less than 4.5 Å between the opposite charges.

**Hydrophobic Interaction Analysis**

**File:** ```hydrophobic_interaction.py```

**Description:** Calculates hydrophobic interactions between hydrophobic residues. Interactions are defined as distances less than 6.5 Å between hydrophobic groups.

## Usage
Each script is standalone and can be executed independently. Ensure that:

Input files (e.g., PDB structures) are in the correct format and location as specified in the scripts.
File paths are adjusted according to your directory structure.
To run a script, use the following command as an example:

```python hydrogen_bond.py```

