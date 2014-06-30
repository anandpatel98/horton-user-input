import user_input_parsing as uip
import atom_attributes as aa
import numpy as np


# Selects defaults based off the molecule and type of calculation the user specified

# This list holds all the keywords that are selected by default (i.e. keywords that the user didn't specify)
unspecified_keywords= []

# Holds all of the supported functionals in horton
functional_list = ['lda_x', 'lda_c_wigner', 'lda_c_rpa', 'lda_c_hl', 'lda_c_gl', 'lda_c_xalpha', 'lda_c_vwn', 'lda_c_vwn_rpa', 'lda_c_pz', 'lda_c_pz_mod', 'lda_c_ob_pz', 'lda_c_pw', 'lda_c_pw_mod', 'lda_c_ob_pw', 'lda_c_2d_amgb', 'lda_c_2d_prm', 'lda_c_vbh', 'lda_c_1d_csc', 'lda_x_2d', 'lda_xc_teter93', 'lda_x_1d', 'lda_c_ml1', 'lda_c_ml2', 'lda_c_gombas', 'lda_c_pw_rpa', 'lda_c_1d_loos', 'lda_c_rc04', 'lda_c_vwn_1', 'lda_c_vwn_2', 'lda_c_vwn_3', 'lda_c_vwn_4', 'lda_k_tf', 'lda_k_lp', 'gga_c_op_xalpha', 'gga_c_op_g96', 'gga_c_op_pbe', 'gga_c_op_b88', 'gga_c_ft97', 'gga_c_spbe', 'gga_x_ssb_sw', 'gga_x_ssb', 'gga_x_ssb_d', 'gga_xc_hcth_407p', 'gga_xc_hcth_p76', 'gga_xc_hcth_p14', 'gga_xc_b97_gga1', 'gga_xc_hcth_a', 'gga_x_bpccac', 'gga_c_revtca', 'gga_c_tca', 'gga_x_pbe', 'gga_x_pbe_r', 'gga_x_b86', 'gga_x_herman', 'gga_x_b86_mgc', 'gga_x_b88', 'gga_x_g96', 'gga_x_pw86', 'gga_x_pw91', 'gga_x_optx', 'gga_x_dk87_r1', 'gga_x_dk87_r2', 'gga_x_lg93', 'gga_x_ft97_a', 'gga_x_ft97_b', 'gga_x_pbe_sol', 'gga_x_rpbe', 'gga_x_wc', 'gga_x_mpw91', 'gga_x_am05', 'gga_x_pbea', 'gga_x_mpbe', 'gga_x_xpbe', 'gga_x_2d_b86_mgc', 'gga_x_bayesian', 'gga_x_pbe_jsjr', 'gga_x_2d_b88', 'gga_x_2d_b86', 'gga_x_2d_pbe', 'gga_c_pbe', 'gga_c_lyp', 'gga_c_p86', 'gga_c_pbe_sol', 'gga_c_pw91', 'gga_c_am05', 'gga_c_xpbe', 'gga_c_lm', 'gga_c_pbe_jrgx', 'gga_x_optb88_vdw', 'gga_x_pbek1_vdw', 'gga_x_optpbe_vdw', 'gga_x_rge2', 'gga_c_rge2', 'gga_x_rpw86', 'gga_x_kt1', 'gga_xc_kt2', 'gga_c_wl', 'gga_c_wi', 'gga_x_mb88', 'gga_x_sogga', 'gga_x_sogga11', 'gga_c_sogga11', 'gga_c_wi0', 'gga_xc_th1', 'gga_xc_th2', 'gga_xc_th3', 'gga_xc_th4', 'gga_x_c09x', 'gga_c_sogga11_x', 'To be used with hyb_gga_x_SOGGA11-X', 'gga_x_lb', 'gga_xc_hcth_93', 'gga_xc_hcth_120', 'gga_xc_hcth_147', 'gga_xc_hcth_407', 'gga_xc_edf1', 'gga_xc_xlyp', 'gga_xc_b97', 'gga_xc_b97_1', 'gga_xc_b97_2', 'gga_xc_b97_d', 'gga_xc_b97_k', 'gga_xc_b97_3', 'gga_xc_pbe1w', 'gga_xc_mpwlyp1w', 'gga_xc_pbelyp1w', 'gga_xc_sb98_1a', 'gga_xc_sb98_1b', 'gga_xc_sb98_1c', 'gga_xc_sb98_2a', 'gga_xc_sb98_2b', 'gga_xc_sb98_2c', 'gga_x_lbm', 'gga_x_ol2', 'gga_x_apbe', 'gga_k_apbe', 'gga_c_apbe', 'gga_k_tw1', 'gga_k_tw2', 'gga_k_tw3', 'gga_k_tw4', 'gga_x_htbs', 'gga_x_airy', 'gga_x_lag', 'gga_xc_mohlyp', 'gga_xc_mohlyp2', 'gga_xc_th_fl', 'gga_xc_th_fc', 'gga_xc_th_fcfo', 'gga_xc_th_fco', 'gga_c_optc', 'gga_k_vw', 'gga_k_ge2', 'gga_k_golden', 'gga_k_yt65', 'gga_k_baltin', 'gga_k_lieb', 'gga_k_absr1', 'gga_k_absr2', 'gga_k_gr', 'gga_k_ludena', 'gga_k_gp85', 'gga_k_pearson', 'gga_k_ol1', 'gga_k_ol2', 'gga_k_fr_b88', 'gga_k_fr_pw86', 'gga_k_dk', 'gga_k_perdew', 'gga_k_vsk', 'gga_k_vjks', 'gga_k_ernzerhof', 'gga_k_lc94', 'gga_k_llp', 'gga_k_thakkar', 'gga_x_wpbeh', 'gga_x_hjs_pbe', 'gga_x_hjs_pbe_sol', 'gga_x_hjs_b88', 'gga_x_hjs_b97x', 'gga_x_ityh', 'hyb_gga_xc_b3pw91', 'hyb_gga_xc_b3lyp', 'hyb_gga_xc_b3p86', 'hyb_gga_xc_o3lyp', 'hyb_gga_xc_mpw1k', 'hyb_gga_xc_pbeh', 'hyb_gga_xc_b97', 'hyb_gga_xc_b97_1', 'hyb_gga_xc_b97_2', 'hyb_gga_xc_x3lyp', 'hyb_gga_xc_b1wc', 'hyb_gga_xc_b97_k', 'hyb_gga_xc_b97_3', 'hyb_gga_xc_mpw3pw', 'hyb_gga_xc_b1lyp', 'hyb_gga_xc_b1pw91', 'hyb_gga_xc_mpw1pw', 'hyb_gga_xc_mpw3lyp', 'hyb_gga_xc_sb98_1a', 'hyb_gga_xc_sb98_1b', 'hyb_gga_xc_sb98_1c', 'hyb_gga_xc_sb98_2a', 'hyb_gga_xc_sb98_2b', 'hyb_gga_xc_sb98_2c', 'hyb_gga_x_sogga11_x', 'hyb_gga_xc_hse03', 'hyb_gga_xc_hse06', 'hyb_gga_xc_hjs_pbe', 'hyb_gga_xc_hjs_pbe_sol', 'hyb_gga_xc_hjs_b88', 'hyb_gga_xc_hjs_b97x', 'hyb_gga_xc_cam_b3lyp', 'hyb_gga_xc_tuned_cam_b3lyp', 'hyb_gga_xc_bhandh', 'hyb_gga_xc_bhandhlyp', 'hyb_gga_xc_mb3lyp_rc04']

# Variables that hold stuff needed to check for relativistic effects
relativistic_atoms = ['co', 'ni', 'cu', 'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr', 'rb', 'sr', 'y', 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in', 'sn', 'sb', 'te', 'i', 'xe', 'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb', 'lu', 'hf', 'ta', 'w', 're', 'os', 'ir', 'pt', 'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn', 'fr', 'ra', 'ac', 'th', 'pa', 'u', 'np', 'pu', 'am', 'cm', 'bk', 'cf', 'es', 'fm', 'md', 'no', 'lr', 'rf', 'db', 'sg', 'bh', 'hs', 'mt', 'ds', 'rg', 'cn', 'uut', 'fl', 'uup', 'lv', 'uus', 'uuo']
relativistic_print = []
relativistic_effect= False


# default profiles for DFT and HF, can be modified through use of if statements.

dft_profile= {"calculation" : uip.user_input["calculation"], "basis": "cc-pvtz", "scf": "ediis2", "iterations" : 300, "energy": 10e-7, "functional": 'gga_x_b88 gga_c_p86', "grid": 'medium', "verbosity" : "medium", "charge" : 0, "spin" : "singlet", "restriction": None, "coordinate_array": None, "atom_array": None, "integration" : "becke", 'functional_header' : None}  #unify functional
hf_profile= {"calculation" : uip.user_input["calculation"], "basis": "cc-pvtz", "scf": "ediis2", "iterations" : 300, "energy": 10e-8, "functional": None , "grid": 'medium' , "verbosity" : "medium", "charge" : 0, "spin" : "singlet", "restriction": None, "coordinate_array": None, "atom_array": None }
user_input_defaults=None


# Check for relativistic atoms (in the atom array) -> should this be under atom attributes?

for x in uip.atom_array:
    if x in relativistic_atoms:
        relativistic_print.append(''.join(x))
        relativistic_output = ', '.join(relativistic_print)
        relativistic_effect = True

if relativistic_effect:
    print "Warning, %s must be treated relativistically to get accurate output" %(relativistic_output)

            
if uip.user_input["calculation"] == "dft":

    # overwrites any 'blank' spots in the user input with the appropriate defaults. Note that it leaves the original user input intact.
    for key in dft_profile:
        if key in uip.user_input:
            dft_profile[key]= uip.user_input[key]
        if key not in uip.user_input:
            unspecified_keywords.append(key)

            
    # If statements that modify the default profiles. If the user has specified something, then the if statements shouldn't overwrite any of their options (what the 'and' conditions are for). 

    if aa.diatomic and 'basis' in unspecified_keywords:
        dft_profile["basis"] = 'aug-' + dft_profile["basis"]
    if  aa.even_electrons and "restriction" in unspecified_keywords:
        dft_profile["restriction"] = "restricted"
        print "Warning, type of calculation changed to %s" %(dft_profile["restriction"])
    if not aa.even_electrons and "restriction"  in unspecified_keywords:
        dft_profile["restriction"] = "unrestricted"
        print "Warning, type of calculation changed to %s" %(dft_profile["restriction"])
    if not "spin" in uip.user_input:
        dft_profile["spin"]=  "doublet"
        print "Warning, the default spin multiplicity of a(n) %s will be used" %(dft_profile["spin"]) # change the warning

    user_input_defaults = dft_profile
        
                           
if uip.user_input["calculation"] == "hf":

# overwrites any 'blank' spots in the user input with the appropriate defaults 
    for key in hf_profile:
        if key in uip.user_input:
            hf_profile[key]= uip.user_input[key]
        if key not in uip.user_input:
            unspecified_keywords.append(key) 
            
    if aa.diatomic and 'basis' in unspecified_keywords:
        hf_profile["basis"] = 'aug-' + hf_profile["basis"]
    if aa.even_electrons and "restriction" not in uip.user_input:
        hf_profile["restriction"] = "restricted"
        print "Warning, type of restriction changed to %s" %(hf_profile["restriction"])
    if not aa.even_electrons and "restriction" not in uip.user_input:
        uip.user_input["restriction"] = "unrestricted"
        print "Warning, type of calculation changed to %s" %(hf_profile["restriction"])
    if not "spin" in uip.user_input:
        hf_profile["spin"]= "doublet"
        print "Warning, the default spin multiplicity of a(n) %s will be used" %(hf_profile["spin"])
    user_input_defaults = hf_profile
  
if uip.test:
    print 'default_input_selection %s' %(user_input_defaults)






