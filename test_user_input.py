import os, sys
import numpy as np
from user_input_revised import get_filename, read_file, check_value, get_attributes, check_relativistic, select_defaults, set_verbosity, initialize_functionals, get_nalpha_nbeta, initialize_molecule, initialize_grid, initialize_scf_solver, exchange_term, initialize_horton, initialize_restricted_hf, initialize_unrestricted_hf, initialize_restricted_dft, initialize_unrestricted_dft

# TODO Need to come up with a test for this
def test_get_filename():
    pass


def test_read_file():

    assert read_file(filename = 'test_read.txt').keys() == ['restriction', 'basis', 'calculation', 'multiplicity', 'energy', 'scf', 'integration', 'coordinate_array', 'charge', 'grid', 'iterations', 'atom_array', 'verbosity']


def test_get_attributes():

    user_input = {'atom_array': np.array(['h','h'])}
    assert get_attributes(user_input) == {'diatomic' : True, 'even_electrons' : True, 'total_electrons': 2}


# Prints out error message, how to suppress this?
def test_check_relativistic():
    
    devnull = open(os.devnull, 'w')
    sys.stdout, sys.stderr = devnull, devnull

    user_input = {'atom_array': np.array(['u','u'])}
    assert check_relativistic(user_input) == ['u']


# Prints out error message, how to supress?
def test_select_defaults():

    devnull = open(os.devnull, 'w')
    sys.stdout, sys.stderr = devnull, devnull

    user_input = {'calculation': 'dft', 'atom_array': np.array(['h','h']), 'coordinate_array': np.array(['1.00000','1.00000'], dtype = float)}
    diatomic = True
    even_electrons = True
    assert select_defaults(user_input, diatomic, even_electrons).keys() == ['restriction', 'calculation', 'energy', 'scf', 'functional', 'coordinate_array', 'grid', 'iterations', 'spin', 'basis', 'multiplicity', 'verbosity', 'integration', 'charge', 'atom_array']



def test_set_verbosity():

    user_input_defaults = {'verbosity': 'high'}
    assert set_verbosity(user_input_defaults) == 'high'


def test_initialize_functionals():
    
    user_input_defaults = {'functional' : 'hyb_gga_xc_o3lyp lda_x', 'restriction': 'restricted'}
    assert len(initialize_functionals(user_input_defaults)) == 2
    assert initialize_functionals.hybrid
    

def test_get_nalpha_nbeta():
    
    user_input_defaults = {'atom_array' : np.array([1,1]), 'charge' : 0, 'spin' : 'singlet'}
    assert get_nalpha_nbeta(user_input_defaults) == {'nalpha' : 1, 'nbeta' : 1}


def test_initialize_molecule():
    
    user_input_defaults = {"calculation" : 'dft', "basis": "cc-pvtz", "scf": "ediis2", "iterations" : 300, "energy": 10e-7, "functional": 'gga_x_b88 gga_c_p86', "grid": 'medium', "verbosity" : "medium", "charge" : 0, "spin" : "singlet", "restriction": None, "coordinate_array": None, "atom_array": np.array([1,1]), "integration" : "becke", 'coordinate_array': np.array([[0,0,0],[0,0,0]])}
    assert initialize_molecule(user_input_defaults)

def test_initialize_grid():
    
    user_input_defaults = {"calculation" : 'dft', "basis": "cc-pvtz", "scf": "ediis2", "iterations" : 300, "energy": 10e-7, "functional": 'gga_x_b88 gga_c_p86', "grid": 'medium', "verbosity" : "silent", "charge" : 0, "spin" : "singlet", "restriction": None, "coordinate_array": None, "atom_array": np.array([1,1]), "integration" : "becke", 'coordinate_array': np.array([[0,0,0],[0,0,0]])}
    set_verbosity(user_input_defaults) # Stops Horton from outputting anything
    mol = initialize_molecule(user_input_defaults)
    
    assert initialize_grid(user_input_defaults, mol)



def test_initialize_scf_solver():
    
    user_input_defaults = {"calculation" : 'dft', "basis": "cc-pvtz", "scf": "ediis2", "iterations" : 300, "energy": 10e-7, "functional": 'gga_x_b88 gga_c_p86', "grid": 'medium', "verbosity" : "silent", "charge" : 0, "spin" : "singlet", "restriction": None, "coordinate_array": None, "atom_array": np.array([1,1]), "integration" : "becke", 'coordinate_array': np.array([[0,0,0],[0,0,0]])}
    mol = initialize_molecule(user_input_defaults)
    
    assert initialize_scf_solver(user_input_defaults, mol)


def test_exchange_term():

    from horton import RLibXCHybridGGA
    
    hybrid = RLibXCHybridGGA('xc_o3lyp')
    restriction = 'restricted'
    assert exchange_term


# TODO Should the energy be compared to a standard calculation here?
def test_initialize_horton():
    pass

    

test_read_file()
test_get_attributes()
test_check_relativistic()
test_select_defaults()
test_set_verbosity()
test_initialize_functionals()
test_get_nalpha_nbeta()
test_initialize_molecule()
test_initialize_grid()
test_initialize_scf_solver()
test_exchange_term()
