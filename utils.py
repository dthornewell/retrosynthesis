import sqlite3
import random
from typing import List, Tuple
from dataclasses import dataclass

import pandas as pd
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdkit import DataStructs
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from template_extractor import extract_from_reaction
from rdenzyme import rdchiralRun
import numpy as np
from scorer import SCScorer


class Retrosim:
    def __init__(self, reference_data_path='RHEA_atom_mapped_timepoint_7_success.pkl', retrobiocat_path='data/retrobiocat_database.pkl'):
        """Initialize RDEnzyme for similarity-based retrosynthesis analysis.
        
        For RDEnzyme, the template cache for each analysis is maintained in jx_cache.
        For RetroBioCat, it is maintained in template_cache
        
        Parameters:
        reference_data_path (str): Path to reference reaction database pickle file
        
        """
        self.jx_cache = {}
        self.reference_data = None
        self.load_reference_data(reference_data_path)

        # for retrobiocat
        self.template_set = pd.read_pickle(retrobiocat_path)
        self.template_cache = {}
        
        # for delta scscore
        self.model = SCScorer()
        self.model.restore('data/model.ckpt-10654.as_numpy.json.gz')


    def load_reference_data(self, file_path):
        """Load and process reference reaction database."""
        self.reference_data = pd.read_pickle(file_path)
        self.reference_data['prod_fp'] = [self.calculate_fingerprint(smi) for smi in self.reference_data['prod_smiles']]                           
        
    @staticmethod
    def calculate_fingerprint(smiles):
        """Calculate Morgan fingerprint for a given SMILES string."""
        fp = AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smiles), 2, useChirality=True, useFeatures=True)
        return fp
             
    def deltascscore(self, prod, react):
        react_list = react.split(".")
        prod_list = prod.split(".")

        react_scores = []
        for smi in react_list:
            can_smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
            (smi, sco) = self.model.get_score_from_smi(smi)
            react_scores.append(sco)
        react_score = max(react_scores)

        prod_scores = []
        for smi in prod_list:
            can_smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
            (smi, sco) = self.model.get_score_from_smi(smi)
            prod_scores.append(sco)
        prod_score = max(prod_scores)

        return (prod_score - react_score).item()
        
    def RetroBioCat(self, prod_smiles):

        can_prod_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(prod_smiles), canonical=True)
        prod = rdchiralReactants(can_prod_smiles)

        # results storage list
        results = []

        # loop through the template set
        for idx, name, rxn_smarts, rxn_type in self.template_set.itertuples():
            if rxn_smarts in self.template_cache:
                rxn = self.template_cache[rxn_smarts]
            else:
                # If not in cache, create and store it
                rxn = rdchiralReaction(rxn_smarts)
                self.template_cache[rxn_smarts] = rxn

            # apply the template
            outcomes = rdchiralRun(rxn, prod, combine_enantiomers=False)
            
            for precursors in outcomes:
                score = self.deltascscore(can_prod_smiles, precursors)
                results.append((name, precursors, score))

        return results     
    
    def single_step_retro(self, target_molecule, max_precursors=50, debug=True, retrobiocat=False):
        """
        Perform single-step retrosynthesis analysis on the target molecule.
        
        Parameters:
        molecule_idx: index of molecule in csv file
        datasub_test: csv file with target molecule and target rxn
        datasub: reference pkl file (with Rhea ID rxns)
        max_precursors: get top 50 most similar rxns
        debug: just to help characterize problems
          
        """
        
        can_target_molecule = canonicalize_smiles(target_molecule)
        product_smiles = [can_target_molecule]
        
        #loads product SMILES into RDKit object
        ex = Chem.MolFromSmiles(target_molecule) 
        #loads product SMILES into RDChiral object
        rct = rdchiralReactants(target_molecule)

        
        if debug:
            print(f"Analyzing product: {product_smiles[0]}") 
        
       
        # Calculate similarities
        fp = self.calculate_fingerprint(product_smiles[0])
        sims = DataStructs.BulkDiceSimilarity(fp, [fp_ for fp_ in self.reference_data['prod_fp']])
        
        # Sort similarity metric in the reverse order, from most to least similar
        js = np.argsort(sims)[::-1]
        probs = {}
        rhea_id = {} # Store the best Rhea ID for each precursor
        rhea_history = {}
        delta_scscore = {}
        
        # Look into each similar rxn in js
        for ji, j in enumerate(js[:max_precursors]):
            jx = self.reference_data.index[j]
            current_rhea_id = self.reference_data['id'][jx]
            
            if debug and ji < 5:
                print(f"\nPrecedent {ji+1}")
                print(f"Similarity score: {sims[j]}")
                print(f"Reference reaction: {self.reference_data['rxn_smiles'][jx]}")
            
            if jx in self.jx_cache:
                (rxn, template, rcts_ref_fp) = self.jx_cache[jx] 
            else:
                try:
                    rxn_smiles = self.reference_data['rxn_smiles'][jx]

                    if isinstance(rxn_smiles, list):
                        rxn_smiles = rxn_smiles[0]
                    elif "']" in rxn_smiles:
                        rxn_smiles = rxn_smiles[2:-2]

                    rct_0, rea_0, prd_0 = rxn_smiles.split(' ')[0].split('>')
                    
                    # Extract template
                    reaction = {'reactants': rct_0,'products': prd_0,'_id': self.reference_data['id'][jx]}
                    template = extract_from_reaction(reaction)

                    #Load into rdChiralReaction
                    rxn = rdchiralReaction(template['reaction_smarts'])

                    #get the reactants to compute reactant fingerprint
                    prec_rxn = self._get_precursor_goal(rxn_smiles)
                    
                    # get rcts reference fingerprint
                    rcts_ref_fp = self.calculate_fingerprint(prec_rxn)

                    #Save for future use
                    self.jx_cache[jx] = (rxn, template, rcts_ref_fp)

                    if debug and ji < 5:
                        print(f"Template: {template['reaction_smarts']}")
                except:
                    continue
                
            try:
                # Run retrosynthesis
                outcomes = rdchiralRun(rxn, rct, combine_enantiomers=False)
            except Exception as e:
                print(e)
                outcomes = []

            if debug and ji < 5:
                print(f"Number of outcomes: {len(outcomes)}")

            # Process outcomes
            for precursors in outcomes:
                precursors = canonicalize_smiles(precursors)
                precursors_fp = self.calculate_fingerprint(precursors)
                precursors_sim = DataStructs.BulkDiceSimilarity(precursors_fp, [rcts_ref_fp])[0]

                overall_score = precursors_sim * sims[j]

                if precursors not in delta_scscore:
                    score = self.deltascscore(can_target_molecule, precursors)
                    delta_scscore[precursors] = score
                
                if precursors not in rhea_history:
                    rhea_history[precursors] = []

                rhea_history[precursors].append({'rhea_id': current_rhea_id, 'score': overall_score})

                # If this precursor structure was already found through a different template/reaction
                if precursors in probs:
                    probs[precursors] = max(probs[precursors], overall_score)
                else:
                    probs[precursors] = overall_score
                    rhea_id[precursors] = current_rhea_id

                if debug and ji < 5:
                    print(f"Found precursor: {precursors}")
                    print(f"Score: {overall_score}")
                
        
        # Rank results
        ranked_output = []
        
        
        # Sort RHEA histories by score for each precursor
        for precursor in rhea_history:
            rhea_history[precursor].sort(key=lambda x: x['score'], reverse=True)
        
        for r, (prec, prob) in enumerate(sorted(probs.items(), key=lambda x:x[1], reverse=True)):
            ranked_output.append((
                "RHEA", # Source
                int(rhea_id[prec]),  # RHEA ID
                prec,  # precursor SMILES
                prob, # probability
            ))

        if retrobiocat:
            retrobiocat_output = self.RetroBioCat(product_smiles[0])
            for name, prec, score in retrobiocat_output:
                if prec not in ranked_output:
                    ranked_output.append(("RetroBioCat", name, prec, 0))  
                    if prec not in delta_scscore:
                        delta_scscore[prec] = score

        scores = [delta_scscore[x[2]] for x in ranked_output]
        
        #return found_rank, product_smiles, ranked_output, prec_goal
        return ranked_output, scores
    
    def _get_precursor_goal(self, rxn_smiles):
        """Extract and process precursor goal from reaction SMILES."""
        if isinstance(rxn_smiles, list):
            rxn_smiles = rxn_smiles[0]
        reactants = rxn_smiles.split('>')[0]
        prec_goal = Chem.MolFromSmiles(reactants)
        [a.ClearProp('molAtomMapNumber') for a in prec_goal.GetAtoms()]
        prec_goal = Chem.MolToSmiles(prec_goal, True)
        return Chem.MolToSmiles(Chem.MolFromSmiles(prec_goal), True)

@dataclass
class Reaction:
    def __init__(self, name, smarts, react_type, precursors:List[str]):
        self.name = name
        self.smarts = smarts
        self.react_type = react_type
        self.precursors = precursors

    def __eq__(self, other): 
        if not isinstance(other, Reaction):
            return False
        return self.name == other.name and self.smarts == other.smarts and self.react_type == other.react_type and self.precursors == other.precursors


def canonicalize_smiles(smi):
    """
    Canonicalize mol SMILES
    """
    try:
        canon_smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi), canonical=True)
    except:
        print(f"ERROR: Cannot Canonicalize {smi}")
        canon_smi = smi
    return canon_smi


def check_in(smile:str, cursor:sqlite3.Cursor, table_name:str) -> bool:
    """
    Check if a SMILES string exists in the database.
    """
    cursor.execute(f'SELECT 1 FROM {table_name} WHERE SMILES = ?', (smile,))
    result = cursor.fetchone()
    return bool(result)

def check_buyable(smile:str, buyables:sqlite3.Cursor, excluded:sqlite3.Cursor) -> bool:
    """
    Check if a SMILES string is buyable.
    """
    if check_in(smile, excluded, 'excluded'):
        return False
    return check_in(smile, buyables, 'buyable')

def weight_normalizer(weights:List) -> List:
    """
    Normalize a list of weights using sigmoid function
    """
    return [1/(1 + np.exp(-weight)) for weight in weights]
    


# Choice functions
def complete_random(choices:List):
    """
    Choose a random element from a list
    """
    if not choices:
        return None
    return random.choice(choices)

def weighted_random(choices:List, weights:List):
    """
    Choose a random element from a list with weights
    """
    if not choices:
        return None
    return random.choices(choices, weights=weights)[0]
