#!/usr/bin/env python
# coding: utf-8

##
## Check all minimized structures regarding the [Isomorphic AssertionError](https://github.com/choderalab/rna_bgsu/issues/1) issue
##


import os, sys, shutil
import glob as glob
import numpy as np
import logging
import warnings
import pickle
import MDAnalysis as mda
from pprint import pprint
from itertools import product

#from openeye import oechem
from openff.qcsubmit.common_structures import QCSpec, PCMSettings, DriverEnum, SCFProperties
from openff.qcsubmit.factories import OptimizationDatasetFactory, BasicDatasetFactory
from openff.qcsubmit.workflow_components import ComponentResult
from openff.toolkit.topology import Molecule
from qcelemental.models.results import WavefunctionProtocolEnum
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit import rdBase
IPythonConsole.drawOptions.addAtomIndices = True
IPythonConsole.drawOptions.addStereoAnnotation = True

from openff.toolkit._version import get_versions
# Warnings that tell us we have undefined stereo and charged molecules
logging.getLogger("openff.toolkit").setLevel(logging.ERROR)
warnings.simplefilter("ignore")

from openff.toolkit.utils import GLOBAL_TOOLKIT_REGISTRY
from openff.toolkit.utils import OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper, BuiltInToolkitWrapper

print(rdBase.rdkitVersion)
print(get_versions())
print(GLOBAL_TOOLKIT_REGISTRY.registered_toolkits)



def get_basesets(string, n):
    """
    Generate base sequence
    """
    base_set = string    
    arr = []
    for a in list(product(base_set, repeat=n)):
        seq = ''.join(a)
        arr.append(seq)
            
    return arr



def pdb_rdkit(files, suffix, stereo_reference):
    rdmols = []
    #for file in tqdm(files):
    for f in files:
        flag = 0
        basename = os.path.basename(f).split('.pdb')[0]
        try:
            u = mda.Universe(f)
            u.add_TopologyAttr("elements", u.atoms.types)
            rdmol = u.atoms.convert_to("RDKIT")  # preserves atom order
            flag = 1
        except:
            print("#PDBConvertError: Could not convert {} to rdkit molecule".format(basename))
            shutil.move(f, f + ".PDBConvertError")

        # check stereocenters
        if flag == 1:
            status = check_stereocenters(f, rdmol, stereo_reference)    
            if status == "success":
                rdmol.SetProp('_Name', basename)
                rdmols.append(rdmol)

    #assert len(files) == len(rdmols)
    print("#{}/{} files converted to rdkit molecules".format(len(rdmols), len(files)))
    
    # save rdkit molecule object
    with open('rdmols_{}.pkl'.format(suffix), 'wb') as pkl:
        pickle.dump(rdmols, pkl, protocol=4)
    
    return rdmols



def check_stereocenters(f, rdmol, stereo_reference):
    """
    exclude unexpected stereocenters
    """
    c = Chem.FindMolChiralCenters(rdmol, force=True, includeUnassigned=True, useLegacyImplementation=False)
    #print(x)

    # exclude backbone phosphorouses
    s = [ x[1] for x in c if rdmol.GetAtomWithIdx(x[0]).GetAtomicNum() != 15 ] 
    s = ''.join(s)

    if s != stereo_reference:
        basename = os.path.basename(f).split('.pdb')[0]
        print("#StereoWarning: {} does not match expected stereo pattern (expected {} but got {})".format(basename, stereo_reference, s))
        shutil.move(f, f + ".StereoWarning")
        return "false"
    else:
        return "success"



def check_bondstereo(f, rdmol):
    """
    only used for structures from triplebase database
    """
    for rdb in rdmol.GetBonds():
        tag = rdb.GetStereo()
        if tag != Chem.BondStereo.STEREONONE:
            basename = os.path.basename(f).split('.pdb')[0]
            print("#BondStereoError: {} does not match expected bond stereo (expected {} but got {})".format(basename, "STEREONONE", tag))
            shutil.move(f, f + ".BondStereoError")
            return "false"
            
    return "success"



# Note that we need to some workaround to solve sterecenters of backbone phosphorouses that are defined randomly (because one bond to O is technically single and the other is double)  
# Related issue: [Isomorphic Assertion Error](https://github.com/choderalab/rna_bgsu/issues/1)
def create_mols(rdmols, stereo_reference):
    #mols = [ Molecule.from_rdkit(rdmol, allow_undefined_stereo=False, hydrogens_are_explicit=True) for rdmol in tqdm(rdmols) ]    
    mols = [ Molecule.from_rdkit(rdmol, allow_undefined_stereo=False, hydrogens_are_explicit=True) for rdmol in rdmols ]    
    for mol in mols:
        for atom in mol.atoms:
            if (atom.atomic_number == 15) and (atom.stereochemistry is not None):
                atom.stereochemistry = None
                #print('stripping P stereo')

    # check if stereocenters are same as the reference
    for i, mol in enumerate(mols):
        atom = mol.atoms
        s = [ a.stereochemistry for a in atom if a.stereochemistry is not None ]
        s = ''.join(s)        
        assert s == stereo_reference, "{} does not match expected sequence (expected {} but got {})".format(rdmols[i].GetProp('_Name'), stereo_reference, s)

    return mols



def create_dataset(mols):
    #factory = BasicDatasetFactory(driver=DriverEnum.gradient, 
    #                          qc_specifications={ 'default': QCSpec(method='wb97m-d3bj', 
    #                                                                basis='def2-tzvppd', 
    #                                                                program='psi4', 
    #                                                                spec_name='rna_default', 
    #                                                                spec_description='RNA quantum chemistry specification', 
    #                                                                store_wavefunction=WavefunctionProtocolEnum.orbitals_and_eigenvalues, 
    #                                                                implicit_solvent=None, 
    #                                                                maxiter=200, 
    #                                                                scf_properties=[SCFProperties.Dipole, SCFProperties.Quadrupole, SCFProperties.WibergLowdinIndices, 
    #                                                                                SCFProperties.MayerIndices, SCFProperties.MBISCharges], 
    #                                                                keywords={'wcombine': False}) },
    #                         store_wavefunction=WavefunctionProtocolEnum.orbitals_and_eigenvalues)
    
    factory = BasicDatasetFactory(
        driver=DriverEnum.gradient, 
        qc_specifications = {
            'default': QCSpec(method='b3lyp-d3bj', 
                            basis='dzvp', 
                            program='psi4', 
                            spec_name='default', 
                            spec_description='Standard OpenFF optimization quantum chemistry specification', 
                            store_wavefunction="none", 
                            implicit_solvent=None, 
                            maxiter=200, 
                            scf_properties=[SCFProperties.Dipole, SCFProperties.Quadrupole], 
                            ), 
        }
    )
    
    dataset = factory.create_dataset(dataset_name="RNA Trinucleotide Single Point Dataset v1.0", 
                                     molecules=mols, 
                                     tagline="QM dataset for ML", 
                                     description="This is a single point energy calculations of RNA trinucleotides extracted from loop motifs of experimental structures. Detailed description on how the data is generated can be found at XXXX")



def workflow_baseset(basesets, suffix):   
    for baseset in basesets:
        # minimized
        files = glob.glob("../minimized/rep*{}*.pdb".format(baseset))
        print(">{}: {} files found".format(baseset, len(files)))

        # convert pdb to rdkit object
        print(">convert pdb to rdkit object")
        rdmols = pdb_rdkit(files, suffix)

        # check missing stereocenters
        print(">check missing stereocenters")
        missing_stereocenters(rdmols)

        # convert rdkit molecule to Openff-toolkit molecule
        print(">create openff-toolkit molecule")
        mols = create_mols(rdmols)

        # check inchikey
        print(">check unique inchikey")
        check_inchikey(rdmols, mols)

        # compare stereocenters
        print(">compare stereocenters")
        seq = "RSRRRSRRRSRR"
        compare_stereo(rdmols, mols, seq, files)

        # create dataset
        print(">create dataset")
        create_dataset(mols)



def workflow(filenames, suffix, stereo_reference):   
    # minimized
    files = glob.glob(filenames)
    print(">{}: {} files found".format(suffix, len(files)))

    # convert pdb to rdkit object
    print(">convert pdb to rdkit object")
    rdmols = pdb_rdkit(files, suffix, stereo_reference)

    # convert rdkit molecule to Openff-toolkit molecule
    print(">create openff-toolkit molecule")
    mols = create_mols(rdmols, stereo_reference)

    # create dataset
    print(">create dataset")
    create_dataset(mols)



if __name__ == "__main__":
    # check triple base databse
    print("### LOOP MOTIFS")
    workflow(filenames="../minimized/rep*.pdb", suffix="loopMOTIFS", stereo_reference="RSRRRSRRRSRR")
