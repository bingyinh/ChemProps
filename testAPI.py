from nmChemPropsAPI import nmChemPropsAPI as nmCPAPI
nm = nmCPAPI("testnmid")
# case 1: uSMILES reported, uSMILES does not exist in ChemProps, uSMILES does not exist in unknowns
result = nm.searchPolymers({'ChemicalName': 'donotexist', 'Abbreviation':'dne', 'uSMILES': 'ccccccccccccccccccccccccccccccccccccccc7'})
# case 2: uSMILES reported, uSMILES does not exist in ChemProps, uSMILES exists in unknowns, duplicate nmid
result = nm.searchPolymers({'ChemicalName': 'donotexist2', 'Abbreviation':'dne2', 'uSMILES': 'ccccccccccccccccccccccccccccccccccccccc7'})
# case 3: uSMILES reported, uSMILES does not exist in ChemProps, uSMILES exists in unknowns, new nmid
nm = nmCPAPI("testnmid1")
result = nm.searchPolymers({'ChemicalName': 'donotexist', 'TradeName': 'trade', 'uSMILES': 'ccccccccccccccccccccccccccccccccccccccc7'})
# case 4: uSMILES not reported, ChemicalName does not exist in ChemProps, ChemicalName does not exist in unknowns
result = nm.searchPolymers({'ChemicalName': 'donotexist3'})
# case 5: uSMILES not reported, ChemicalName does not exist in ChemProps, ChemicalName exists in unknowns, duplicate nmid
result = nm.searchPolymers({'ChemicalName': 'donotexist3'})
# case 6: uSMILES not reported, ChemicalName does not exist in ChemProps, ChemicalName exists in unknowns, newnmid
nm = nmCPAPI("testnmid2")
result = nm.searchPolymers({'ChemicalName': 'donotexist3', 'TradeName': 'trade3', 'Abbreviation': 'dne3'})
# case 7: find a winning match, reported Abbreviation and TradeName not in ChemProps
result = nm.searchPolymers({'ChemicalName': 'Polyurethane', 'Abbreviation': 'ignoreME', 'TradeName': 'ignoreME!'})
# case 8: find a winning match, regular case
result = nm.searchPolymers({'ChemicalName': 'epoxy', 'Abbreviation':'DGEBA'})
