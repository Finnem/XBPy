
from collections import defaultdict

# lengths are taken from pymol: https://github.com/schrodinger/pymol-open-source/blob/9d3061ca58d8b69d7dad74a68fc13fe81af0ff8e/layer2/AtomInfo.cpp

simplified_bond_lengths = defaultdict(
    lambda: defaultdict(lambda : 1.54, {}), {
    "C" : defaultdict(lambda: 1.54, {
                "I": 2.14,
                "Cl": 1.77,
                "Br": 1.94,
                "F": 1.35,
                "S": 1.82,
                "N": 1.47,
                "O": 1.43,
                "P": 1.84,
                "C": 1.54,
    }),
    "N" : defaultdict(lambda: 1.45, {
                "N": 1.25,
                "O": 1.40,
                "S": 1.75,
    }),
    "O" : defaultdict(lambda: 1.35, {
                "O": 1.35,
                "S": 1.44,
    }),
    "S" : defaultdict(lambda: 1.82, {
                "S": 2.05,
    }),
})

# order is Single, Aromatic, Double, Triple
# if no aromatic is possible it should be the same as double
bond_lengths = defaultdict(
    lambda: defaultdict(lambda : [1.54, 1.34, 1.34, 1.2], {}), {
    "H" : defaultdict(lambda: [1.1, 0, 0, 0], {
                "H" :  [0.8, 0, 0, 0],
                "C" :  [1.2, 0, 0, 0],
                "N" :  [1.02, 0, 0, 0],
                "O" :  [0.96, 0, 0, 0],
                "S" :  [1.34, 0, 0, 0],
    }),
    "C" : defaultdict(lambda: [1.54, 1.41, 1.34, 1.2], {
                "I": [2.14, 0, 0, 0],
                "Cl": [1.77, 0, 0, 0],
                "Br": [1.94, 0, 0, 0],
                "F": [1.35, 0, 0, 0],
                "S": [1.82, 1.73, 1.71, 0],
                "N": [1.47, 1.36, 1.29, 1.16],
                "O": [1.43, 1.36, 1.25, 0],
                "P": [1.84, 0, 0, 0],
                "C": [1.54, 1.4, 1.34, 1.2],
    }),
    "N" : defaultdict(lambda: [1.45, 1.25, 1.25, 0], {
                "N": [1.45, 1.35, 1.25, 0],
                "O": [1.40, 1.3, 1.21, 0],
                "S": [1.75, 1.6, 1.53, 0],
    }),
    "O" : defaultdict(lambda: [1.45, 1.35, 1.35, 0], {
                "O" : [1.45, 1.4, 1.35, 0],
                "S" : [1.75, 1.6, 1.44, 0],
    }),
    "S" : defaultdict(lambda: [1.82, 0, 0, 0], {
                "S": [2.05, 0, 0, 0],
    }),
})

# ideal length list from https://en.wikipedia.org/wiki/Bond_length
ideal_bond_lengths = defaultdict(
    lambda: defaultdict(lambda: [1.54, 1.34, 1.34, 1.20], {}),
    {
        "H": defaultdict(
            lambda: [1.1, 0, 0, 0],  # Default H bonds (primarily single bonds)
            {
                "H": [0.74, 0, 0, 0],    # H-H single bond
                "C": [1.09, 0, 0, 0],    # H-C single bond
                "N": [1.01, 0, 0, 0],    # H-N single bond
                "O": [0.96, 0, 0, 0],    # H-O single bond
                "S": [1.34, 0, 0, 0],    # H-S single bond
            }
        ),
        "C": defaultdict(
            lambda: [1.54, 1.34, 1.34, 1.20],  # Default C bonds
            {
                "I": [2.14, 0, 0, 0],         # C-I single bond
                "Cl": [1.77, 0, 0, 0],        # C-Cl single bond
                "Br": [1.94, 0, 0, 0],        # C-Br single bond
                "F": [1.35, 0, 0, 0],         # C-F single bond
                "S": [1.82, 1.73, 1.71, 1.60],# C-S bonds: single, aromatic, double, triple
                "N": [1.47, 1.36, 1.29, 1.16],# C-N bonds: single, aromatic, double, triple
                "O": [1.43, 1.36, 1.25, 0],    # C-O bonds: single, aromatic, double
                "P": [1.84, 0, 0, 0],         # C-P single bond
                "C": [1.54, 1.34, 1.34, 1.20],# C-C bonds: single, aromatic, double, triple
            }
        ),
        "N": defaultdict(
            lambda: [1.45, 1.25, 1.25, 0],  # Default N bonds
            {
                "N": [1.45, 1.35, 1.25, 0],    # N-N bonds: single, aromatic, double
                "O": [1.40, 1.30, 1.21, 0],    # N-O bonds: single, aromatic, double
                "S": [1.75, 1.60, 1.53, 0],    # N-S bonds: single, aromatic, double
            }
        ),
        "O": defaultdict(
            lambda: [1.45, 1.35, 1.35, 0],  # Default O bonds
            {
                "O": [1.45, 1.40, 1.35, 0],    # O-O bonds: single, aromatic, double
                "S": [1.75, 1.60, 1.44, 0],    # O-S bonds: single, aromatic, double
            }
        ),
        "S": defaultdict(
            lambda: [1.82, 0, 0, 0],        # Default S bonds
            {
                "S": [2.05, 0, 0, 0],          # S-S single bond
            }
        ),
    }
)

# make dict kommutative
next_bond_length = {}
for atom in bond_lengths:
    for other_atom in bond_lengths[atom]:
        if other_atom not in next_bond_length:
            next_bond_length[other_atom] = defaultdict(lambda: [1.54, 1.34, 1.34, 1.2], {})
        next_bond_length[other_atom][atom] = bond_lengths[atom][other_atom]
for atom in next_bond_length:
    if atom not in bond_lengths:
        bond_lengths[atom] = next_bond_length[atom]
    else:
        bond_lengths[atom].update(next_bond_length[atom])
# make dict kommutative
next_bond_length = {}
for atom in ideal_bond_lengths:
    for other_atom in ideal_bond_lengths[atom]:
        if other_atom not in next_bond_length:
            next_bond_length[other_atom] = defaultdict(lambda: [1.54, 1.34, 1.34, 1.2], {})
        next_bond_length[other_atom][atom] = ideal_bond_lengths[atom][other_atom]
for atom in next_bond_length:
    if atom not in ideal_bond_lengths:
        ideal_bond_lengths[atom] = next_bond_length[atom]
    else:
        ideal_bond_lengths[atom].update(next_bond_length[atom])




def jump_to_nth_last_line(file_path, n):
    with open(file_path, 'rb') as file:
        file.seek(0, 2)  # Go to the end of the file
        file_size = file.tell()
        lines_found = 0
        buffer_size = 1024
        buffer = b""
        
        for position in range(file_size, 0, -buffer_size):
            seek_position = max(0, position - buffer_size)
            file.seek(seek_position)
            buffer = file.read(min(buffer_size, position))
            
            lines_found += buffer.count(b'\n')
            
            if lines_found >= n:
                break
        
        # Now that we have the buffer, find the exact position of the 4th last line
        lines = buffer.splitlines()
        line_to_jump_to = lines[-n] if len(lines) >= n else lines[0]  # Safeguard against short files
        
        return line_to_jump_to.decode('utf-8')  # Convert byte string to regular string

import numpy as np
from itertools import combinations
from collections import defaultdict


def get_covalent_radius(element_symbol):
    """
    Retrieve the covalent radius for the given element symbol.
    Uses the 'periodictable' library.
    """
    from rdkit.Chem import PeriodicTable
    periodic_table = PeriodicTable.GetPeriodicTable()
    try:
        return periodic_table.GetRCovalent(element_symbol)
    except:
        raise ValueError(f"Unknown element symbol: {element_symbol}")


possible_geometries = {
    "H": {
        1: {
            "angles": [],
            "bond_orders": [1],
        }
    },
    "C": {
        1: {
            "angles": [],
            "bond_orders": [3],
        },
        2: {
            "angles": [180],
            "bond_orders": [1, 2, 3],
        },
        3: {
            "angles": [120],
            "bond_orders": [1, 2],
        },
        4: {
            "angles": [109.5],
            "bond_orders": [1],
        }
    },
    "N": {
        1: {
            "angles": [],
            "bond_orders": [3],
        },
        2: {
            "angles": [180],
            "bond_orders": [2, 3],
        },
        3: {
            "angles": [107],
            "bond_orders": [1],
        },
        4: {
            "angles": [109.5],
            "bond_orders": [1],
        }
    },
    "O": {
        1: {
            "angles": [],
            "bond_orders": [2],
        },
        2: {
            "angles": [104.5],
            "bond_orders": [1],
        },
        3: {
            "angles": [107],
            "bond_orders": [1],
        }
    },
    "F": {
        1: {
            "angles": [],
            "bond_orders": [1],
        }
    },
    "Cl": {
        1: {
            "angles": [],
            "bond_orders": [1],
        }
    },
    "Br": {
        1: {
            "angles": [],
            "bond_orders": [1],
        }
    },
    "I": {
        1: {
            "angles": [],
            "bond_orders": [1],
        }
    },
    "B": {
        3: {
            "angles": [120],
            "bond_orders": [1],
        },
        4: {
            "angles": [109.5],
            "bond_orders": [1],
        }
    },
    "Si": {
        4: {
            "angles": [109.5],
            "bond_orders": [1],
        }
    },
    "P": {
        3: {
            "angles": [93],  # Due to lone pairs
            "bond_orders": [1],
        },
        4: {
            "angles": [109.5],
            "bond_orders": [1],
        },
        5: {
            "angles": [90, 120],
            "bond_orders": [1],
        }
    },
    "S": {
        2: {
            "angles": [119],
            "bond_orders": [2],
        },
        4: {
            "angles": [109.5],
            "bond_orders": [1],
        },
        6: {
            "angles": [90],
            "bond_orders": [1],
        }
    }
}
