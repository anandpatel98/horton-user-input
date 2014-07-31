from user_input_revised import get_filename, read_file, check_value, get_attributes, check_relativistic, select_defaults, set_verbosity, initialize_functionals, get_nalpha_nbeta, initialize_molecule, initialize_grid, initialize_scf_solver, exchange_term, initialize_horton, initialize_restricted_hf, initialize_unrestricted_hf, initialize_restricted_dft, initialize_unrestricted_dft

# Initializes the appropriate functions
user_input= read_file(filename = get_filename())
attributes = get_attributes(user_input)
check_relativistic(user_input)
user_input_defaults= select_defaults(user_input, diatomic = attributes['diatomic'], even_electrons = attributes['even_electrons'])
set_verbosity(user_input_defaults)
if user_input['calculation'] == 'dft':
    initialize_functionals(user_input_defaults)
mol = initialize_molecule(user_input_defaults)
initialize_horton(user_input_defaults, mol)
