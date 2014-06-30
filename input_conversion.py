import user_input_parsing as uip
import atom_attributes as aa
import default_selection as ds
from horton import *

import numpy as np

# this module processes user input in order to initialize certain parts of horton

# sets the verbosity

if ds.user_input_defaults['verbosity'] == 'silent':
    log.set_level(log.silent)
elif ds.user_input_defaults['verbosity'] == 'warning':
    log.set_level(log.warning)
elif ds.user_input_defaults['verbosity'] == 'low':
    log.set_level(log.low)
elif ds.user_input_defaults['verbosity'] == 'medium':
    log.set_level(log.medium)
elif ds.user_input_defaults['verbosity'] == 'high':
    log.set_level(log.high)
elif ds.user_input_defaults['verbosity'] == 'debug':
    log.set_level(log.debug)
else:
    raise ValueError, "The level of verbosity specified is not a valid option."
           
# should the user specify the units that the coordinate matrix is in? This assumes that the units are in angstroms and are converted to au


ds.user_input_defaults['coordinate_array']= angstrom*np.array(ds.user_input_defaults['coordinate_array'])

# takes the atoms and creates a vector with atomic numbers. the counter variable is there in order to access array elements.

i=0
for s in np.nditer(ds.user_input_defaults["atom_array"], ["refs_ok"]):
    ds.user_input_defaults["atom_array"][i] = aa.atom_electrons[str(s)]
    i += 1
ds.user_input_defaults["atom_array"]= ds.user_input_defaults["atom_array"].astype(int)

# Suite of statements for hartree-fock

if ds.user_input_defaults['calculation'] == 'hf':
    mol= Molecule(coordinates = np.squeeze(ds.user_input_defaults['coordinate_array']), energy= float(ds.user_input_defaults['energy']),  numbers= np.squeeze(ds.user_input_defaults['atom_array']),  basis= ds.user_input_defaults['basis'], restriction= ds.user_input_defaults['restriction'], iterations= int(ds.user_input_defaults['iterations']), intgrid= ds.user_input_defaults['grid'])


# Statements for dft

if ds.user_input_defaults["calculation"]== 'dft':

    mol= Molecule(coordinates = np.squeeze(ds.user_input_defaults['coordinate_array']), energy= float(ds.user_input_defaults['energy']),  numbers= np.squeeze(ds.user_input_defaults['atom_array']),  basis= ds.user_input_defaults['basis'], restriction= ds.user_input_defaults['restriction'], iterations= int(ds.user_input_defaults['iterations']), intgrid= ds.user_input_defaults['grid'])
    

    # Initializes the right libxc functional. Does this for hybrid functionals, or seperate exchange/ correlation functionals.
    lib=[]
    hybrid = None
    functional_dict = {}
    functionals_seperated = ds.user_input_defaults['functional'].split(' ')
    ex_corr = False
    if len(functionals_seperated) == 2:
        ex_corr == True
    if len(functionals_seperated) > 2:
        raise ValueError, "You cannot specify more than two functionals" 
    for entry in functionals_seperated:
        functional_temp = entry.split('_')
        functional_dict[''.join(functional_temp[0].lower())] = '_'.join(functional_temp[1:])
    for key in functional_dict.keys():
        if key == 'lda':
            lib.append(RLibXCLDA(functional_dict[key]))
        elif key == 'gga':
            lib.append(RLibXCGGA(functional_dict[key]))
        elif key == 'hyb':
            hybrid = functional_dict[key]
            clean_functional = string.split(functional_dict[key])
            lib.append(RLibXCHybridGGA('_'.join(clean_functional[2:])))
        else:
            raise ValueError, "The type of functional specified is not a valid option currently."

# Defines grid variable (in horton). There are 4 grids in horton: becke, atomic, radial and line?

    if ds.user_input_defaults['integration'] == "becke":
        grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate= False, agspec= mol.intgrid )
    elif ds.user_input_defaults['integration'] == "atomic":
        grid = AtomicGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate= False, agspec= mol.intgrid )
    elif ds.user_input_defaults['integration'] == "line":
        grid = LineGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers)
    else:
        raise ValueError, "The type of integration grid specified is not a valid option currently."

# Calls the appropriate SCF solver for both HF and DFT. Note that some SCF solvers use an expansion, where others use a density matrix

if ds.user_input_defaults['scf'] == 'oda':
    scf_solver = ODASCFSolver(mol.energy, mol.iterations)
    expansion = False
elif ds.user_input_defaults['scf'] == 'scf':
    scf_solver = PlainSCFSolver(mol.energy, mol.iterations)
    expansion = True
elif ds.user_input_defaults['scf'] == 'ediis2':
    scf_solver = EDIIS2SCFSolver(mol.energy, mol.iterations)
    expansion = False
elif ds.user_input_defaults['scf'] == 'ediis':
    scf_solver = EDIISSCFSolver(mol.energy, mol.iterations)
    expansion = False
elif ds.user_input_defaults['scf'] == 'cdiis':
    scf_solver = CDIISSCFSolver(mol.energy, mol.iterations)
    expansion = False
else:
    raise ValueError, "The type of SCF solver specified is not a valid option currently."


# takes the spin multiplicity and the charge, and works out the number of alpha and beta electrons <---- Take another look at this part of the code, there is a built in function to do this in horton

charge = int(ds.user_input_defaults["charge"])
multiplicity = ds.user_input_defaults["spin"]

multiplicity_conversion ={'singlet' : 1, 'doublet' : 2, 'triplet' : 3, 'quartet' : 4, 'quintet' : 5} # add more entries as needed

int_multiplicity = multiplicity_conversion[multiplicity] 

paired_electrons= aa.total_electrons-(int_multiplicity-1)
nalpha= (int_multiplicity -1) + (paired_electrons/2)
nbeta = (paired_electrons/2)

ds.user_input_defaults["nalpha"] = nalpha
ds.user_input_defaults["nbeta"] = nbeta

if uip.test:
    print "nalpha %s" %(nalpha)
    print 'nbeta %s' %(nbeta)
    if 'lib' in locals():
        print "lib_initialized"
    
    


