# Extra code that has been saved in case it is needed

for x in user_input:
    if x in ["calculation", "scf", "energy", "correlation", "exchange", "functional", "basis", "1dm", "2dm", "wfn", "verbosity", "charge", "spin", "restriction"]:
        v= ''.join(user_input[x])
        print v + "hi"
    else:
        continue
    if x == 'calculation':
        pass
    if x == 'scf':
        pass
    if x == 'correlation':
        pass
    if x == 'exchange':
        pass
    if x== 'functional':
        pass
    if x== 'basis':
        pass
    if x== '1dm':
        pass
    if x== '2dm':
        pass
    if x== 'wfn':
        pass
    if x== 'verbosity':
        pass
    if x== 'charge':
        pass
    if x== 'spin':
        pass
    if x== 'restriction':
        if not v == "restricted" or v== "unrestricted":
            print "Error goes here"
            exit()
