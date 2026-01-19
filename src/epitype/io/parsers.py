"""Structure file parsers."""
from __future__ import annotations

from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser as BioPDBParser
from Bio.PDB.Structure import Structure as BioStructure

from epitype.core.structure import Atom, Chain, Residue, Structure


def parse_pdb(filepath: str | Path) -> Structure:
    """
    Parse PDB file into Structure.

    Args:
        filepath: Path to PDB file

    Returns:
        Parsed Structure object
    """
    filepath = Path(filepath)
    parser = BioPDBParser(QUIET=True)
    bio_structure = parser.get_structure(filepath.stem, str(filepath))

    return _convert_biopython_structure(bio_structure, filepath.stem)


def parse_mmcif(filepath: str | Path) -> Structure:
    """
    Parse mmCIF file into Structure.

    Args:
        filepath: Path to mmCIF file

    Returns:
        Parsed Structure object
    """
    filepath = Path(filepath)
    parser = MMCIFParser(QUIET=True)
    bio_structure = parser.get_structure(filepath.stem, str(filepath))

    return _convert_biopython_structure(bio_structure, filepath.stem)


def parse_structure(filepath: str | Path) -> Structure:
    """
    Parse structure file (auto-detect format).

    Args:
        filepath: Path to PDB or mmCIF file

    Returns:
        Parsed Structure object
    """
    filepath = Path(filepath)
    suffix = filepath.suffix.lower()

    if suffix in {".pdb", ".ent"}:
        return parse_pdb(filepath)
    elif suffix in {".cif", ".mmcif"}:
        return parse_mmcif(filepath)
    else:
        # Try PDB first, then mmCIF
        try:
            return parse_pdb(filepath)
        except Exception:
            return parse_mmcif(filepath)


def _convert_biopython_structure(bio_struct: BioStructure, name: str) -> Structure:
    """Convert BioPython structure to our Structure class."""
    structure = Structure(name=name)

    atom_index = 0
    residue_index = 0

    # Take first model only
    model = bio_struct[0]

    for bio_chain in model:
        chain_id = bio_chain.id
        chain = Chain(chain_id=chain_id)

        for bio_residue in bio_chain:
            # Skip water and heteroatoms (unless modified amino acid)
            if bio_residue.id[0] != " ":
                # Keep modified residues (e.g., MSE)
                if bio_residue.id[0] not in {"W", "H_HOH"}:
                    pass  # Include modified residues
                else:
                    continue

            residue = Residue(
                index=residue_index,
                name=bio_residue.resname,
                chain_id=chain_id,
                seq_num=bio_residue.id[1],
                insertion_code=bio_residue.id[2].strip(),
            )

            for bio_atom in bio_residue:
                element = bio_atom.element.strip() if bio_atom.element else bio_atom.name[0]

                atom = Atom(
                    index=atom_index,
                    name=bio_atom.name,
                    element=element,
                    coords=np.array(bio_atom.coord, dtype=np.float64),
                    residue_index=residue_index,
                    chain_id=chain_id,
                    b_factor=bio_atom.bfactor,
                    occupancy=bio_atom.occupancy,
                )
                residue.atoms.append(atom)
                atom_index += 1

            chain.residues.append(residue)
            residue_index += 1

        if chain.residues:
            structure.chains[chain_id] = chain

    return structure
