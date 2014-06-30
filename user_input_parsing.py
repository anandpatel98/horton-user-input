#!/usr/bin/env python

# This module takes a user input file and parses the keywords and their corresponding entries into a dictionary variable called user_input. Includes sanity checks to ensure the user doesn't enter invalid keywords. ***Note that only valid keywords are checked in this module, but the input is not checked.

import sys
import numpy as np

filename = sys.argv[1]

test = False
try:
    if sys.argv[2] == 'test':
        test = True
except:
    pass

user_input= {}
atom_list=[]
coordinate_list=[]
line_count = 0

keywords= ["calculation", "scf", "correlation", "exchange", "functional", "grid", "basis", "pseudopotential", "1dm", "2dm", "wfn", "verbosity", "charge", "coordinates", "restriction", "iterations", "energy", "gradient", "spin", "integration"]

atoms = ['h', 'he', 'li', 'be', 'b', 'c', 'n', 'o', 'f', 'ne', 'na', 'mg', 'al', 'si', 'p', 's', 'cl', 'ar', 'k', 'ca', 'sc', 'ti', 'v', 'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr', 'rb', 'sr', 'y', 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in', 'sn', 'sb', 'te', 'i', 'xe', 'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu', 'hf', 'ta', 'w', 're', 'os', 'ir', 'pt', 'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn', 'fr', 'ra', 'ac', 'th', 'pa', 'u', 'np', 'pu', 'am', 'cm', 'bk', 'cf', 'es', 'fm', 'md', 'no', 'lr', 'rf', 'db', 'sg', 'bh', 'hs', 'mt', 'ds', 'rg', 'cn', 'uut', 'fl', 'uup', 'lv', 'uus', 'uuo']

# Adds entries into a dictionary, with the keyword as a key. Does not take coordinates into the dictionary.
with open(filename) as handle:
    for line in handle:
	    line_count += 1
	    if line == '\n':
	        continue
	    tempinput = line.split()
	    key= str(tempinput[0]).lower()
	    value = ' '.join(tempinput[1:]).lower()
	    if key in user_input.keys(): # checks for repeated keywords (excluding coordinates)
	        raise KeyError, 'One or more of the keywords has been repeated on line %s' % (line_count) 
	        exit()
	    if not key in keywords:
	        print "The keyword on line %s is not reconized." % (line_count)
	        exit()
	    if "coordinates" in key:
	        break
	    if len(value) == 0:
	        print "The keyword on line %s does not have a corresponding value" % (line_count)
	        exit()
	    user_input[str(key)] = value
# Imports coordinates (includes check to make sure only the coordinates are last)
    for line in handle:
        line_count += 1
        if line == '\n':
	        continue
        tempinput = line.split()
        try:
            key= str(tempinput[0]).lower()
        except:
            continue
        value= tempinput [1:]
        if len(value) != 3:
	        print "Please check line %s of your input file. An error was found in the coordinate matrix." % (line_count)
	        exit()       
        if len(value) == 0:
	        print "The keyword on line %s does not have a corresponding value" % (line_count)
	        exit()
        find_coordinate= "coordinates" in key
	# Error checking portion, are coordinates last?, is the coordinate keyword repeated?, did you make up atoms?      	
        if key in user_input.keys(): # checks for repeated keywords (excluding coordinates)
            raise KeyError, "One or more of the keywords has been repeated on line %s" % (line_count) 
            exit()
        atom_check= key in atoms
        if atom_check:
	        atom_array_temp = key # ensures that the key variable is an array
	        atom_list.append(atom_array_temp)
	        coordinate_array_temp = np.array(value, dtype=float)
	        coordinate_list.append(coordinate_array_temp)
        else: 
            print "Please check line %s of your input file. An error was found in the coordinate matrix." % (line_count)
            exit()

if 'calculation' not in user_input.keys():
    raise KeyError, 'The calculation keyword must be specified to run a calculation.'
    
if user_input["calculation"] == "hf":
    if "functional" in user_input or "exchange" in user_input or "correlation" in user_input:
        print "The functional, exchange or correlation keywords cannot be sepcified for hf calculations"
        exit()
	
# Stacks up lists into an array, then adds the completed arrays to the dictionary
atom_array= np.vstack(atom_list)
coordinate_array = np.vstack(coordinate_list)
user_input["atom_array"] = atom_array
user_input["coordinate_array"] = coordinate_array

# This is only run during unit testing. It is supposed to strip whitespace from dictionary output for ease of parsing.

if test:
    def removew(d):
        for k, v in d.iteritems():
            if isinstance(v, dict):
                removew(v)
            else:
                pass
    removew(user_input)
    print 'user_input %s' %(user_input)









