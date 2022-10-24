# Create trinucleotide RNA dataset for QCArchive submission


Download json/csv files from [RNA BGSU](https://www.bgsu.edu/research/rna.html)
-----
- copied from [rna-bgsu-dinucleotide](https://github.com/kntkb/rna-bgsu-dinucleotide)


Extract structure data using [RNA BGSU APIs](https://www.bgsu.edu/research/rna/APIs.html)
-----
- copied from [rna-bgsu-dinucleotide](https://github.com/kntkb/rna-bgsu-dinucleotide)


Check downloaded cif/pdb files -> copied from rna-bgsu-dinucleotide
-----
- copied from [rna-bgsu-dinucleotide](https://github.com/kntkb/rna-bgsu-dinucleotide)



Split loops into three consecutive bases and perform clustering for each base set
------

- `split_motif.ipynb`  
    - Hairpins, internal, and junction loops were split into dinucleotides (three connected bases).
    - Atoms `P`, `OP1`, and `OP2` were deleted from 5' base.

- `cluster_motifs.ipynb`
    - Each triple base set was clustered with bottom-up (agglomerative) hierarchical clustering with average linking.
    - Euclidean distances between set of pre-defined atom pairs were used as input features.



Add hydrogen atoms and minimize structure
------

- `script/openmm_implicit_minimizer.py`
    - Added hydrogen prior to minimization
    - Implicit solvent applied (GBN2)
    - Heavy atom restraint applied (30 kcal/mol/A^2)
    - Add suffix to pdb if warnings and/or problems were found and move problematic strucutres to a different directory
        - Minimized structures are stored in `minimized/pdb`
        - PDB structures with problems were moved to `minimized/error`


Test creating dataset prior to QCArchive submission
------

- `qca-dataset-submission_TEST/check_structures.py`
    - convert minimized pdb structures into rdkit molecule object and check stererocenters and stereobonds
    - **pdb structures that passed the check and conveted to rdkit mol object are saved as pickle**
    - `rdmols_loopMOTIFS.pkl` are used for QCArchive submission