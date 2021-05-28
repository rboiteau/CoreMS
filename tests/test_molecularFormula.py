__author__ = "Yuri E. Corilo"
__date__ = "Jul 22, 2019"

import sys
sys.path.append(".")

import pytest
from corems.molecular_formula.factory.MolecularFormulaFactory import MolecularFormula     
from corems.encapsulation.constant import Labels
from copy import deepcopy

def test_molecular_formula():
    
    '''test the MolecularFormula class and the calculation of isotopologues'''
    
    formula_dict = {'C':10, 'H':0, 'O':10,'Cl':2, Labels.ion_type: 'Radical'}
    
    ion_charge = 1 
    formula_obj = MolecularFormula(formula_dict, ion_charge)
    print("ion_type", formula_obj.ion_type)
    assert round(formula_obj.mz_calc,2) == round(349.886303060457,2)
    
    min_abundance, current_abundance = 1,1 
    #print(min_abundance, current_abundance)
    isotopologues = list(formula_obj.isotopologues(0.01, current_abundance, 500))
    
    assert round(isotopologues[0].mz_calc,2) == round(351.883352980637,2)
    assert round(isotopologues[0].prob_ratio,2) == round(0.6399334750069298,2)
    assert isotopologues[0].string == 'C10 O10 Cl1 37Cl1'
    
    formula_obj.ion_type = 'RADICAL'
    formula_obj.kmd
    formula_obj.kendrick_mass
    formula_obj.knm
    formula_obj.atoms_qnt('C')
    formula_obj.class_label
    formula_obj.atoms_symbol('13C')
    
    formula_str = 'C10 H21 N1'
    formula_obj = MolecularFormula(formula_str, ion_charge)
    '''
    for isotopologue_obj in formula_obj.isotopologues(0.01, current_abundance):
        
        print("formula:", isotopologue_obj.string, 
              "mz_calc:", isotopologue_obj.mz_calc,
              "prob_ratio:", isotopologue_obj.prob_ratio)
      '''

if __name__ == "__main__":
      test_molecular_formula()
   

    