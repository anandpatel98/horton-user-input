
# Looks at the atom matrix in order to determine different attributes of the molecule (used later in default selection). This module determines if a molecule is diatomic, or if any atoms present need to be treated relativistically.

import user_input_parsing as uip
import numpy as np

total_electrons= 0
atom_electrons = {'h' : 1, 'he' : 2, 'li' : 3, 'be' : 4 , 'b' : 5 , 'c' : 6, 'n' : 7 , 'o' : 8 , 'f' : 9, 'ne': 10, 'na':11, 'mg':12, 'al':13, 'si':14, 'p':15, 's':16, 'cl':17, 'ar':18, 'k':19, 'ca':20, 'sc':21, 'ti':22, 'v':23, 'cr':24, 'mn':25, 'fe':26, 'co':27, 'ni':28, 'cu': 29, 'zn':30, 'ga':31, 'ge':32, 'as':33, 'se':34, 'br':35, 'kr':36, 'rb':37, 'sr':38, 'y':39, 'zr':40, 'nb':41, 'mo':42, 'tc':43, 'ru':44, 'rh':45, 'pd':46, 'ag':47, 'cd':48, 'in':49, 'sn':50, 'sb':51, 'te':52, 'i':53, 'xe':54, 'cs':55, 'ba':56, 'la':57, 'ce':58, 'pr':59, 'nd':60, 'pm':61, 'sm':62, 'eu':63, 'gd':64, 'tb':65, 'dy':66, 'ho':67, 'er':68, 'tm':69, 'yb':70, 'lu':71, 'hf':72, 'ta':73, 'w':74, 're':75, 'os':76, 'ir':77, 'pt':78, 'au':79, 'hg':80, 'tl':81, 'pb':82, 'bi':83, 'po':84, 'at':85, 'rn':86, 'fr':87, 'ra':88, 'ac':89, 'th':90, 'pa':91, 'u':92, 'np':93, 'pu':94, 'am':95, 'cm':96, 'bk':97, 'cf':98, 'es':99, 'fm':100, 'md':101, 'no':102, 'lr':103, 'rf':104, 'db':105, 'sg':106, 'bh':107, 'hs':108, 'mt':109, 'ds':110, 'rg':111, 'cn':112, 'uut':113, 'fl':114, 'uup':115, 'lv':116, 'uus':117, 'uuo':118}

if len(uip.atom_array)== 2:
    if uip.atom_array[0] == uip.atom_array[1]:
        diatomic = True
    else:
        diatomic = False
else:
    diatomic = False

for atom in np.nditer(uip.atom_array):
    electrons_temp = atom_electrons[str(atom)]
    total_electrons += electrons_temp
# even_electrons stores if the electron number is even= true or odd= false.

if total_electrons%2==0:
    even_electrons = True
else:
    even_electrons = False

if uip.test:
    print 'diatomic %s' %(diatomic)
    print 'even_electrons %s' %(even_electrons)
    print 'total_electrons %s' %(total_electrons)
    