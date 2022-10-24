#!/usr/bin/env python
# coding: utf-8

# ## Minimize structures
# 
# Add missing hydrogen using pdbfixer. 5' base will be automatically detected so there is no need to change the first residue name.
# 
# **Notes**
# Alternatively, missing hydrogens could be added by pdb4amber. Note that the first residue needs to be modify (e.g. A --> A5).
# >import pdb4amber  
# >pdb4amber.run(arg_pdbin='input.pdb', arg_pdbout='pdb4amber.pdb', arg_add_missing_atoms=True)
# 
# **References**
# - [OpenMM-Tricks-and-Recipes](https://github-wiki-see.page/m/ParmEd/ParmEd/wiki/OpenMM-Tricks-and-Recipes)
# - [OPENMM_TUTORIAL](https://gpantel.github.io/assets/PDF/OpenMM_Tutorial.pdf)  

import os, sys, shutil
import pathlib
import glob as glob
import numpy as np
import re
import warnings
import mdtraj as md
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from openmm.app import PDBFile
from pdbfixer import PDBFixer


def minimize(output_path, file):
    
    """
    load file
    """
    fixer = PDBFixer(filename=file)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(True)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    #fixer.addSolvent(fixer.topology.getUnitCellDimensions())
    #PDBFile.writeFile(fixer.topology, fixer.positions, open('pdbfixer.pdb', 'w'))
    
    
    """
    setup system
    """
    forcefield = ForceField('amber14/RNA.OL3.xml', 'implicit/gbn2.xml')
    system = forcefield.createSystem(fixer.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
    
    # heavy atom restraint
    force = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter("k", 30.0*kilocalories_per_mole/angstroms**2)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    
    for i, atom in enumerate(fixer.topology.atoms()):
        if atom.element.symbol != "H":
            atom_crd = fixer.positions[i]
            force.addParticle(i, atom_crd.value_in_unit(nanometers))
    system.addForce(force)
    
    
    """
    minimize structure
    """
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(fixer.topology, system, integrator)
    simulation.context.setPositions(fixer.positions)
    # http://getyank.org/latest/api/multistate_api/index.html?highlight=minimi#yank.multistate.multistatesampler.MultiStateSampler.minimize
    #simulation.minimizeEnergy(maxIterations=100, tolerance=1.0*kilojoules_per_mole/nanometers)
    simulation.minimizeEnergy(maxIterations=100)
    minpositions = simulation.context.getState(getPositions=True).getPositions()
    
    """
    save pdb
    """    
    try:
        basename = os.path.basename(file)
        PDBFile.writeFile(fixer.topology, minpositions, open(os.path.join(output_path, basename), 'w'))   
        
        # check if simulation can be run properly
        #simulation.reporters.append(PDBReporter('/Users/takabak/Desktop/dump.pdb', 1))
        #simulation.reporters.append(StateDataReporter(stdout, 1, step=True, potentialEnergy=True, temperature=True))
        #simulation.reporters.append(StateDataReporter(stdout, 1, step=False, potentialEnergy=False, temperature=False))
        #simulation.step(10)
    except:
        print("{}: Check structure!!".format(basename))

    del simulation, system


def test(f):
    fixer = PDBFixer(filename=f)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.findMissingAtoms()
    
    if fixer.missingAtoms:
        print("{}: missing atoms".format(f))
        #shutil.move(f, f + ".warning")
        status = "failed"
    else:
        status = "success"

    return status




if __name__ == "__main__":
    base_path = os.path.dirname(os.path.abspath("__file__")).strip('scripts')
    output_path = os.path.join(base_path, "minimized")

    # motif
    _path = os.path.join(base_path, "pdb", "motif", "cluster", "triplebase")
    files = glob.glob(_path + "/*/centroid/rep*.pdb")    
    print(">{} files found".format(len(files)))
    
    for file in files:
        print(os.path.basename(file))
        status = test(file)
        if status == "success":
            minimize(output_path, file)
        else:
            print(status)
