# to do for this module: the other unit tests must be written, along with testing horton to get the same number to a particular precision. check that the numbers match up. test files should be h2 rather than h2o to save time

import subprocess

# Checks that the default system for hartree fock is working, and that horton is initializing properly
def test_hf_plain():
    process= subprocess.Popen(['python', 'horton_initialization.py', 'test_hf_plain.txt'], stdout= subprocess.PIPE)
    output = process.communicate()[0]
    output_list = output.split()
    i = 0
    for entry in output_list:
        if output_list[i] == 'total':
            energy = float(output_list[i+1])
        i += 1
# Statement that checks energy against the 'standard' value should go here

# Tests that the default system for dft is working properly, and that horton is initializing properly
def test_dft_plain():
    process= subprocess.Popen(['python', 'horton_initialization.py', 'test_dft_plain.txt'], stdout= subprocess.PIPE)
    output = process.communicate()[0]
    output_list = output.split()
    i = 0
    for entry in output_list:
        if output_list[i] == 'total':
            energy = float(output_list[i+1])
        i =+ 1

# this function still has a few bugs in it- parsing the dictionary output isn't as easy as it seems!
def test_user_input_parsing():
        import ast 
        process= subprocess.Popen(['python', 'user_input_parsing.py', 'test_user_input_parsing.txt', 'test'], stdout= subprocess.PIPE)
        output = process.communicate()[0]
        output_list = output.split()
        i = 0
        user_input = {}
        values_list=[]
        for entry in output_list:
            if output_list[i] == 'user_input':
                output_list_dict = ''.join(output_list[i+1:])
                #user_input = ast.literal_eval(output_list_dict)  converts string to dict, then assigns it
            i += 1
        print output_list_dict
        for value in user_input.values():
            values_list.append(value)
        assert values_list == ['dft', 'scf', '10e-7', '300', 'fine', 'becke', '6-31g*', 'restricted', 'doublet', '0', 'low']
    
def test_atom_attributes():
    process= subprocess.Popen(['python', 'atom_attributes.py', 'test_attributes.txt', 'test'], stdout= subprocess.PIPE)
    output = process.communicate()[0]
    output_list = output.split()
    i = 0
    for entry in output_list:
        if entry == 'diatomic':
            diatomic = output_list[i+1]
        if output_list[i] == 'total_electrons':
            total_electrons = output_list[i+1]
        if output_list[i] == 'even_electrons':
            even_electrons = output_list[i+1]
        i += 1
    assert diatomic
    assert total_electrons == '184'
    assert even_electrons

def test_input_conversion():
    process= subprocess.Popen(['python', 'input_conversion.py', 'test_conversion.txt', 'test'], stdout= subprocess.PIPE)
    output = process.communicate()[0]
    output_list = output.split()
    lib = False
    i = 0
    print output_list
    for entry in output_list:
        if output_list[i] == 'nalpha':
            nalpha = output_list[i+1]
        if output_list[i] == 'nbeta':
            nbeta = output_list[i+1]
        if output_list[i] == 'lib_initialized':
            lib = True
        i += 1
    assert nalpha == '5'
    assert nbeta == '5' # this should really be 4, use the built in nalpha/nbeta calculator in horton
    assert lib

# give a set of keywords and values that are converted by this module, and check to ssee that they are changed appropriately.
# also check that particular values are left untouched by the module.


test_input_conversion()

