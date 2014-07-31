import sys
import numpy as np

'''
NOTE: (bug) changing maxiter in the plain scf solver doesn't do anything. Not sure if it's a problem with this code or with horton. Changing threshold energy for the plain scf solver works fine though.
'''

# Gets filename from command line
def get_filename():
    try:
        filename = sys.argv[1]
        return filename
    except:
        raise ValueError, "The input filename has not been specified or is not recognized. Please check the name of your input file and try again."


# reads in input from the file, including coordinates
def read_file(filename):
    keywords= ["calculation", "scf", "correlation", "exchange", "functional", "grid", "basis", "verbosity", "charge", "coordinates", "restriction", "iterations", "energy", "gradient", "multiplicity", "integration"]
    atoms = ['h', 'he', 'li', 'be', 'b', 'c', 'n', 'o', 'f', 'ne', 'na', 'mg', 'al', 'si', 'p', 's', 'cl', 'ar', 'k', 'ca', 'sc', 'ti', 'v', 'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr', 'rb', 'sr', 'y', 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in', 'sn', 'sb', 'te', 'i', 'xe', 'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu', 'hf', 'ta', 'w', 're', 'os', 'ir', 'pt', 'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn', 'fr', 'ra', 'ac', 'th', 'pa', 'u', 'np', 'pu', 'am', 'cm', 'bk', 'cf', 'es', 'fm', 'md', 'no', 'lr', 'rf', 'db', 'sg', 'bh', 'hs', 'mt', 'ds', 'rg', 'cn', 'uut', 'fl', 'uup', 'lv', 'uus', 'uuo']
    
    user_input= {}
    atom_list=[]
    coordinate_list=[]
    line_count = 0
    with open(filename) as handle:
        for line in handle:
            line_count += 1
            if line == '\n':
                continue
            tempinput = line.split()
            key= str(tempinput[0]).lower()
            value = ' '.join(tempinput[1:]).lower()
            if "coordinates" in key:
                break
            check_value(user_input, value, key, keywords, line_count)
            user_input[str(key)] = value
	    # wipe leftover values in the variables before the next loop starts   
        key = None
        value = None
    
# Imports coordinates (includes check to make sure only the coordinates are last)
        for line in handle:
            line_count += 1
            if line == '\n':
	            continue
            tempinput = line.split()
            # Without the try/except, adding a space to a blank line doesn't register as \n and gives an error
            try:
                key= str(tempinput[0]).lower()
            except:
                raise SyntaxError, "The input file should not have blank lines after the coordinate matrix. Please check your input file and try again."
            value= tempinput [1:]
            if len(value) != 3:
	            raise KeyError, "Please check line %s of your input file. An error was found in the coordinate matrix." % (line_count)   
            if len(value) == 0:
	            raise KeyError, "The keyword on line %s does not have a corresponding value" % (line_count)
	        # Ensure the coordinates are last, the coordinate keyword isn't repeated, and the user specifed real atoms 	
            if key in atoms:
	            atom_list.append(key)
	            coordinate_list.append(value)
            elif key == 'coordinates':
	            raise KeyError, 'The coordinate keyword has been repeated on line %s.' %(line_count)
            else:
                raise KeyError, "Please check line %s of your input file. An error was found in the coordinate matrix." % (line_count)
                #TODO change this errror message to something more appropriate
    user_input["atom_array"] = np.array(atom_list)
    user_input["coordinate_array"] = np.array(coordinate_list, dtype = float)
    return user_input
    
    
def check_value(user_input, value, key, keywords, line_count):
    if key in user_input.keys(): # checks for repeated keywords (excluding coordinates)
	    raise KeyError, 'One or more of the keywords has been repeated on line %s' % (line_count) 
    if key not in keywords:
	    raise KeyError, "The keyword on line %s is not recognized." % (line_count)
    if len(value) == 0:
	    raise KeyError, "The keyword on line %s does not have a corresponding value" % (line_count)


# Finds attributes of the molecule (whether it is diatomic, has even electrons and the total number of electrons) 
def get_attributes(user_input):
    total_electrons= 0
    atom_electrons= {'h' : 1, 'he' : 2, 'li' : 3, 'be' : 4 , 'b' : 5 , 'c' : 6, 'n' : 7 , 'o' : 8 , 'f' : 9, 'ne': 10, 'na':11, 'mg':12, 'al':13, 'si':14, 'p':15, 's':16, 'cl':17, 'ar':18, 'k':19, 'ca':20, 'sc':21, 'ti':22, 'v':23, 'cr':24, 'mn':25, 'fe':26, 'co':27, 'ni':28, 'cu': 29, 'zn':30, 'ga':31, 'ge':32, 'as':33, 'se':34, 'br':35, 'kr':36, 'rb':37, 'sr':38, 'y':39, 'zr':40, 'nb':41, 'mo':42, 'tc':43, 'ru':44, 'rh':45, 'pd':46, 'ag':47, 'cd':48, 'in':49, 'sn':50, 'sb':51, 'te':52, 'i':53, 'xe':54, 'cs':55, 'ba':56, 'la':57, 'ce':58, 'pr':59, 'nd':60, 'pm':61, 'sm':62, 'eu':63, 'gd':64, 'tb':65, 'dy':66, 'ho':67, 'er':68, 'tm':69, 'yb':70, 'lu':71, 'hf':72, 'ta':73, 'w':74, 're':75, 'os':76, 'ir':77, 'pt':78, 'au':79, 'hg':80, 'tl':81, 'pb':82, 'bi':83, 'po':84, 'at':85, 'rn':86, 'fr':87, 'ra':88, 'ac':89, 'th':90, 'pa':91, 'u':92, 'np':93, 'pu':94, 'am':95, 'cm':96, 'bk':97, 'cf':98, 'es':99, 'fm':100, 'md':101, 'no':102, 'lr':103, 'rf':104, 'db':105, 'sg':106, 'bh':107, 'hs':108, 'mt':109, 'ds':110, 'rg':111, 'cn':112, 'uut':113, 'fl':114, 'uup':115, 'lv':116, 'uus':117, 'uuo':118}
     
    diatomic = len(user_input['atom_array']) == 2 and user_input['atom_array'][0] == user_input['atom_array'][1]
    for atom in np.nditer(user_input['atom_array']):
        electrons_temp = atom_electrons[str(atom)]
        total_electrons += electrons_temp
    even_electrons = total_electrons%2==0
    return {'diatomic': diatomic, 'even_electrons': even_electrons, 'total_electrons': total_electrons}
    
    
# Checks for relativistic atoms and prints warning messages accordingly    
def check_relativistic(user_input):

    relativistic_atoms = ['co', 'ni', 'cu', 'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr', 'rb', 'sr', 'y', 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in', 'sn', 'sb', 'te', 'i', 'xe', 'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu', 'hf', 'ta', 'w', 're', 'os', 'ir', 'pt', 'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn', 'fr', 'ra', 'ac', 'th', 'pa', 'u', 'np', 'pu', 'am', 'cm', 'bk', 'cf', 'es', 'fm', 'md', 'no', 'lr', 'rf', 'db', 'sg', 'bh', 'hs', 'mt', 'ds', 'rg', 'cn', 'uut', 'fl', 'uup', 'lv', 'uus', 'uuo']
    relativistic_print = []
    relativistic_effect= False
    for x in user_input['atom_array']:
        if x in relativistic_atoms and x not in relativistic_print:
            relativistic_print.append(''.join(x))
            relativistic_output = ', '.join(relativistic_print)
            relativistic_effect = True
    if relativistic_effect:
        print "Warning, %s must be treated relativistically to get accurate output" %(relativistic_output)
        
    return relativistic_print


# Selects appropriate default values if the user didn't select them. Also converts atomic symbols to numbers and certain values to their appropriate type (might need to put that in a seperate function)
def select_defaults(user_input, diatomic, even_electrons):
    from horton import angstrom
    
    atom_electrons= {'h' : 1, 'he' : 2, 'li' : 3, 'be' : 4 , 'b' : 5 , 'c' : 6, 'n' : 7 , 'o' : 8 , 'f' : 9, 'ne': 10, 'na':11, 'mg':12, 'al':13, 'si':14, 'p':15, 's':16, 'cl':17, 'ar':18, 'k':19, 'ca':20, 'sc':21, 'ti':22, 'v':23, 'cr':24, 'mn':25, 'fe':26, 'co':27, 'ni':28, 'cu': 29, 'zn':30, 'ga':31, 'ge':32, 'as':33, 'se':34, 'br':35, 'kr':36, 'rb':37, 'sr':38, 'y':39, 'zr':40, 'nb':41, 'mo':42, 'tc':43, 'ru':44, 'rh':45, 'pd':46, 'ag':47, 'cd':48, 'in':49, 'sn':50, 'sb':51, 'te':52, 'i':53, 'xe':54, 'cs':55, 'ba':56, 'la':57, 'ce':58, 'pr':59, 'nd':60, 'pm':61, 'sm':62, 'eu':63, 'gd':64, 'tb':65, 'dy':66, 'ho':67, 'er':68, 'tm':69, 'yb':70, 'lu':71, 'hf':72, 'ta':73, 'w':74, 're':75, 'os':76, 'ir':77, 'pt':78, 'au':79, 'hg':80, 'tl':81, 'pb':82, 'bi':83, 'po':84, 'at':85, 'rn':86, 'fr':87, 'ra':88, 'ac':89, 'th':90, 'pa':91, 'u':92, 'np':93, 'pu':94, 'am':95, 'cm':96, 'bk':97, 'cf':98, 'es':99, 'fm':100, 'md':101, 'no':102, 'lr':103, 'rf':104, 'db':105, 'sg':106, 'bh':107, 'hs':108, 'mt':109, 'ds':110, 'rg':111, 'cn':112, 'uut':113, 'fl':114, 'uup':115, 'lv':116, 'uus':117, 'uuo':118}
    
    dft_profile= {"calculation" : user_input["calculation"], "basis": "cc-pvtz", "scf": "ediis2", "iterations" : 300, "energy": 10e-7, "functional": 'gga_x_b88 gga_c_p86', "grid": 'medium', "verbosity" : "medium", "charge" : 0, "spin" : "singlet", "restriction": None, "coordinate_array": None, "atom_array": None, "integration" : "becke"}  #unify functional
    hf_profile= {"calculation" : user_input["calculation"], "basis": "cc-pvtz", "scf": "ediis2", "iterations" : 300, "energy": 10e-8, "functional": None , "grid": 'medium' , "verbosity" : "medium", "charge" : 0, "spin" : "singlet", "restriction": None, "coordinate_array": None, "atom_array": None }
    
    unspecified_keywords= []
    
    if user_input["calculation"] == "dft":
        user_input_defaults = dft_profile
    elif user_input["calculation"] == "hf":
        user_input_defaults = hf_profile
    else:
        raise ValueError, 'The type of calculation specified is not currently valid.'
        
    for key in user_input_defaults:
        if key in user_input:
            user_input_defaults[key]= user_input[key]
        if key not in user_input:
            unspecified_keywords.append(key)
            
    if diatomic and 'basis' in unspecified_keywords:
        user_input_defaults["basis"] = 'aug-' + user_input_defaults["basis"]
    if  even_electrons and "restriction" in unspecified_keywords:
        user_input_defaults["restriction"] = "restricted"
        print "Warning, type of calculation changed to %s" %(user_input_defaults["restriction"])
    if not even_electrons and "restriction"  in unspecified_keywords:
        user_input_defaults["restriction"] = "unrestricted"
        print "Warning, type of calculation changed to %s" %(user_input_defaults["restriction"])
    if not "multiplicity" in user_input:
        user_input_defaults["multiplicity"]=  "doublet"
        print "Warning, the default spin multiplicity of a(n) %s will be used" %(user_input_defaults["spin"])
        
    user_input_defaults['energy'] = float(user_input_defaults['energy'])
    user_input_defaults['iterations'] = int(user_input_defaults['iterations'])
    user_input_defaults['coordinate_array']= angstrom*user_input_defaults['coordinate_array']
    user_input_defaults['charge'] = int(user_input_defaults["charge"])
    
   
    i = 0
    for s in np.nditer(user_input_defaults["atom_array"], ["refs_ok"]): #TODO: use enumerate and convert in aa
        user_input_defaults["atom_array"][i] = atom_electrons[str(s)]
        i += 1
    user_input_defaults["atom_array"]= user_input_defaults["atom_array"].astype(int)
    
    return user_input_defaults

# Sets the verbosity of the output
def set_verbosity(user_input_defaults):
    from horton import log
    
    if user_input_defaults['verbosity'] == 'silent': #TODO: 1 statement
        log.set_level(log.silent) # TODO ds.user... change to dummy var
        return 'silent'
    elif user_input_defaults['verbosity'] == 'warning':
        log.set_level(log.warning)
        return 'warning'
    elif user_input_defaults['verbosity'] == 'low':
        log.set_level(log.low)
        return 'low'
    elif user_input_defaults['verbosity'] == 'medium':
        log.set_level(log.medium)
        return 'medium'
    elif user_input_defaults['verbosity'] == 'high':
        log.set_level(log.high)
        return 'high'
    elif user_input_defaults['verbosity'] == 'debug':
        log.set_level(log.debug)
        return 'debug'
    else:
        raise ValueError, "The level of verbosity specified is not a valid option."

# Initializes the appropriate libxc functionals based off what the user specified
def initialize_functionals(user_input_defaults):

    from horton import RLibXCLDA, RLibXCGGA, RLibXCHybridGGA, ULibXCLDA, ULibXCGGA, ULibXCHybridGGA
    
    lib=[]
    initialize_functionals.hybrid = None
    functional_dict = {}
    functionals_seperated = user_input_defaults['functional'].split()
    if len(functionals_seperated) > 2:
        raise ValueError, "You cannot specify more than two functionals" 
    for entry in functionals_seperated:
        functional_temp = entry.split('_')
        functional_header = functional_temp[0]
        if functional_temp[0] == 'hyb':
            functional_remainder = '_'.join(functional_temp[2:])
        elif functional_temp[0] == 'gga' or functional_temp[0] == 'lda':
            functional_remainder = '_'.join(functional_temp[1:])
        else:
            raise ValueError, 'Please ensure your specified functionals are supported'# TODO find better error message
            
        if user_input_defaults['restriction'] ==  'restricted':
            if functional_header == 'lda':
                lib.append(RLibXCLDA(functional_remainder))
            elif functional_header == 'gga':
                lib.append(RLibXCGGA(functional_remainder))
            elif functional_header == 'hyb':
                lib.append(RLibXCHybridGGA(functional_remainder))
                initialize_functionals.hybrid = RLibXCHybridGGA(functional_remainder)
            else:
               raise ValueError, "The type of functional specified is not a valid option currently."
        
        if user_input_defaults['restriction'] ==  'unrestricted':
            if functional_header == 'lda':
                lib.append(ULibXCLDA(functional_remainder))
            elif functional_header == 'gga':
                lib.append(ULibXCGGA(functional_remainder))
            elif functional_header == 'hyb':
                lib.append(ULibXCHybridGGA(functional_remainder))
                initialize_functionals.hybrid = ULibXCHybridGGA(functional_remainder)
            else:
                raise ValueError, "The type of functional specified is not a valid option currently."
        
    return lib


# Determines the number of nalpha and nbeta electrons
def get_nalpha_nbeta(user_input_defaults):
    charge = user_input_defaults["charge"]
    total_electrons = np.sum(user_input_defaults['atom_array'])
    multiplicity = user_input_defaults["spin"]

    multiplicity_conversion ={'singlet' : 1, 'doublet' : 2, 'triplet' : 3, 'quartet' : 4, 'quintet' : 5}

    int_multiplicity = multiplicity_conversion[multiplicity]

    paired_electrons= total_electrons-(int_multiplicity-1)
    nalpha= (int_multiplicity -1) + (paired_electrons/2)
    nbeta = (paired_electrons/2)
    
    return {'nalpha' : nalpha, 'nbeta': nbeta}


# Initializes a molecule object with the appropriate information
def initialize_molecule(user_input_defaults):

    from horton import Molecule
    nalpha_nbeta_dict = get_nalpha_nbeta(user_input_defaults)

    mol= Molecule(coordinates = user_input_defaults['coordinate_array'], energy= user_input_defaults['energy'],  numbers= user_input_defaults['atom_array'],  basis= user_input_defaults['basis'], restriction= user_input_defaults['restriction'], iterations= user_input_defaults['iterations'], intgrid= user_input_defaults['grid'], nalpha= nalpha_nbeta_dict['nalpha'], nbeta = nalpha_nbeta_dict['nbeta'])

    return mol


# Initializes the appropriate grid
def initialize_grid(user_input_defaults, mol):
    
    from horton import BeckeMolGrid

    if user_input_defaults['integration'] == "becke":
        grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate= False, agspec= mol.intgrid )
    else:
        raise ValueError, "The type of integration grid specified is not a valid option currently."
        
    return grid
    # The other grid types are commented out since they don't work as of yet (TODO ask about their uses)
    '''
    elif ds.user_input_defaults['integration'] == "atomic":
        grid = AtomicGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate= False, agspec= mol.intgrid )
    elif ds.user_input_defaults['integration'] == "line":
        grid = LineGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers)
    '''


# Initializes the appropriate SCF solver
def initialize_scf_solver(user_input_defaults, mol):
    
    from horton import ODASCFSolver, PlainSCFSolver, EDIIS2SCFSolver, EDIISSCFSolver, CDIISSCFSolver
    
    if user_input_defaults['scf'] == 'oda':
        scf_solver = ODASCFSolver(mol.energy, mol.iterations)
    elif user_input_defaults['scf'] == 'scf':
        scf_solver = PlainSCFSolver(mol.energy, mol.iterations)
    elif user_input_defaults['scf'] == 'ediis2':
        scf_solver = EDIIS2SCFSolver(mol.energy, mol.iterations)
    elif user_input_defaults['scf'] == 'ediis':
        scf_solver = EDIISSCFSolver(mol.energy, mol.iterations)
    elif user_input_defaults['scf'] == 'cdiis':
        scf_solver = CDIISSCFSolver(mol.energy, mol.iterations)
    else:
        raise ValueError, "The type of SCF solver specified is not a valid option currently."
    
    return scf_solver


# For hybrid fuctionals only. Returns the appropriate fraction of hartree-fock exchange energy
def exchange_term(user_input_defaults, mol, er):
    
    from horton import RExchangeTerm, UExchangeTerm
    
    hybrid = initialize_functionals.hybrid
    restriction = mol.restriction
    if restriction == 'restricted'and hybrid:
        return RExchangeTerm(er, 'x_hf', hybrid.get_exx_fraction())
    if restriction == 'unrestricted'and hybrid:
        return UExchangeTerm(er, 'x_hf', hybrid.get_exx_fraction())
        

#Main function to start up Horton        
def initialize_horton(user_input_defaults, mol):
    
    from horton import get_gobasis, DenseLinalgFactory

    obasis = get_gobasis(mol.coordinates, mol.numbers, mol.basis)
    lf= DenseLinalgFactory(obasis.nbasis)
    
    olp = obasis.compute_overlap(lf)
    kin = obasis.compute_kinetic(lf)
    na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf) # mandatory pseudo_numbers??
    er = obasis.compute_electron_repulsion(lf)
    
    if user_input_defaults['calculation'] == 'hf' and user_input_defaults['restriction'] == 'restricted':
        initialize_restricted_hf(user_input_defaults, mol, lf, olp, kin, na, er, obasis)
    if user_input_defaults['calculation'] == 'hf' and user_input_defaults['restriction'] == 'unrestricted':
        initialize_unrestricted_hf(user_input_defaults, mol, lf, olp, kin, na, er, obasis)
    if user_input_defaults['calculation'] == 'dft' and user_input_defaults['restriction'] == 'restricted':
        lib = initialize_functionals(user_input_defaults)
        initialize_restricted_dft(user_input_defaults, mol, lib, lf, olp, kin, na, er, obasis)
    if user_input_defaults['calculation'] == 'dft' and user_input_defaults['restriction'] == 'unrestricted':
        lib = initialize_functionals(user_input_defaults)
        initialize_unrestricted_dft(user_input_defaults, mol, lib, lf, olp, kin, na, er, obasis)

# Starts up stuff needed for restricted hartree-fock
def initialize_restricted_hf(user_input_defaults, mol, lf, olp, kin, na, er, obasis):
    
    from horton import guess_core_hamiltonian, compute_nucnuc, ROneBodyTerm, RDirectTerm, RExchangeTerm, AufbauOccModel, REffHam

    exp_alpha = lf.create_expansion()
    guess_core_hamiltonian(olp, kin, na, exp_alpha)
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    terms = [
            ROneBodyTerm(kin, 'kin'),
            RDirectTerm(er, 'hartree'),
            RExchangeTerm(er, 'x_hf'),
            ROneBodyTerm(na, 'ne'),
        ]
    ham = REffHam(terms, external)
    occ_model = AufbauOccModel(mol.nalpha)
    occ_model.assign(exp_alpha)
    dm_alpha = exp_alpha.to_dm()
    scf_solver= initialize_scf_solver(user_input_defaults, mol)
    if scf_solver.kind == 'exp':
        scf_solver(ham, lf, olp, occ_model, exp_alpha)
    else:
        scf_solver(ham, lf, olp, occ_model, dm_alpha)


# Starts up stuff for unrestricted hartree-fock
def initialize_unrestricted_hf(user_input_defaults, mol, lf, olp, kin, na, er, obasis):

    from horton import guess_core_hamiltonian, compute_nucnuc, UOneBodyTerm, UDirectTerm, UExchangeTerm, AufbauOccModel, UEffHam

    exp_beta = lf.create_expansion()
    exp_alpha = lf.create_expansion()
    guess_core_hamiltonian(olp, kin, na, exp_alpha, exp_beta)
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    terms = [
        UOneBodyTerm(kin, 'kin'),
        UDirectTerm(er, 'hartree'),
        UExchangeTerm(er, 'x_hf'),
        UOneBodyTerm(na, 'ne'),
    ]
    ham = UEffHam(terms, external)
    occ_model = AufbauOccModel(mol.nalpha, mol.nbeta)
    occ_model.assign(exp_alpha, exp_beta)
    dm_alpha = exp_alpha.to_dm()
    dm_beta= exp_beta.to_dm()
    scf_solver= initialize_scf_solver(user_input_defaults, mol)
    if scf_solver.kind == 'exp':
        scf_solver(ham, lf, olp, occ_model, exp_alpha, exp_beta)
    else:
        scf_solver(ham, lf, olp, occ_model, dm_alpha, dm_beta)


# Starts up stuff for restricted dft
def initialize_restricted_dft(user_input_defaults, mol, lib, lf, olp, kin, na, er, obasis):

    from horton import guess_core_hamiltonian, compute_nucnuc, ROneBodyTerm, RDirectTerm, RExchangeTerm, RGridGroup, AufbauOccModel, REffHam
    
    exp_alpha = lf.create_expansion()
    guess_core_hamiltonian(olp, kin, na, exp_alpha)
    grid = initialize_grid(user_input_defaults, mol)
    
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    terms = [
        ROneBodyTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RGridGroup(obasis, grid, lib),
        ROneBodyTerm(na, 'ne')
    ]
    if exchange_term(user_input_defaults, mol, er):
        terms.append(exchange_term(user_input_defaults, mol, er)) 
    ham = REffHam(terms, external)
    occ_model = AufbauOccModel(mol.nalpha)
    occ_model.assign(exp_alpha)
    dm_alpha = exp_alpha.to_dm()
        
    scf_solver = initialize_scf_solver(user_input_defaults, mol)
    if scf_solver.kind == 'exp':
        scf_solver(ham, lf, olp, occ_model, exp_alpha)
    else:
        scf_solver(ham, lf, olp, occ_model, dm_alpha)


#Starts up stuff for unrestricted dft
def initialize_unrestricted_dft(user_input_defaults, mol, lib, lf, olp, kin, na, er, obasis):

    from horton import guess_core_hamiltonian, compute_nucnuc, UOneBodyTerm, UDirectTerm, UExchangeTerm, UGridGroup, AufbauOccModel, UEffHam

    exp_beta = lf.create_expansion()
    exp_alpha = lf.create_expansion()
    guess_core_hamiltonian(olp, kin, na, exp_alpha, exp_beta)
    grid = initialize_grid(user_input_defaults, mol)
        
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}

    terms = [
        UOneBodyTerm(kin, 'kin'),
        UDirectTerm(er, 'hartree'),
        UGridGroup(obasis, grid, lib),
        UOneBodyTerm(na, 'ne'),
    ]
    if exchange_term(user_input_defaults, mol, er):
        terms.append(exchange_term(user_input_defaults, mol, er))
    ham = UEffHam(terms, external)

    occ_model = AufbauOccModel(mol.nalpha,mol.nbeta)
    occ_model.assign(exp_alpha, exp_beta) 
    dm_alpha = exp_alpha.to_dm()
    dm_beta= exp_beta.to_dm()
        
    scf_solver = initialize_scf_solver(user_input_defaults, mol)

    if scf_solver.kind == 'exp': 
        scf_solver(ham, lf, olp, occ_model, exp_alpha, exp_beta)
    else:
        scf_solver(ham, lf, olp, occ_model, dm_alpha, dm_beta)

