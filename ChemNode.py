import sqlite3
import pandas as pd
from ReactNode import ReactNode
from utils import *
from typing import List
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from rdenzyme import rdchiralRun


EXPLORATION_PARAM = 1.414
MAX_DEPTH = 6

class ChemNode:

    buyables:sqlite3.Cursor = None # Cursor for buyable lookup
    abundants:sqlite3.Cursor = None # Here to rule out simple chemicals (i.e. water, oxygen, ammonium)
    retrobiocat:pd.DataFrame = None # Retrobiocat templates
    analyzer:Retrosim = None # RdEnzyme analyzer

    def __init__(self, smiles:str, depth:int, parent_reaction:ReactNode, buyables:sqlite3.Cursor = None, 
                 templates:pd.DataFrame = None, retrosim:Retrosim = None, abundant:sqlite3.Cursor = None):
        # Chemical data
        self.smiles = smiles
        self.parent_reaction = parent_reaction
        self.depth = depth
        self.reactions:List[ReactNode] = []

        # MCTS variables
        self.possible_reactions:List[Reaction] = []
        self.weights:List[float] = []
        self.visits:int = 0
        self.score:float = 0.0

        # Needed for buyable lookup and reaction generation
        if buyables is not None:
            ChemNode.buyables = buyables
        if abundant is not None:
            ChemNode.abundants = abundant
        if templates is not None:
            ChemNode.retrobiocat = templates
        if retrosim is not None:
            ChemNode.analyzer = retrosim

        self.buyable:bool = check_buyable(smiles, ChemNode.buyables)
        self.solution:bool = self.buyable

        if (not self.solution): # If buyable, no need to generate reactions
            self.generate_reactions_rdenzyme()


    def copy(self) -> 'ChemNode':
        """
        Copy this node
        """
        new_node = ChemNode(self.smiles, self.depth, self.parent_reaction)
        new_node.reactions = []
        new_node.possible_reactions = []
        new_node.visits = self.visits
        new_node.score = self.score
        new_node.solution = self.solution
        return new_node
    
    def generate_reactions_retrobiocat(self):
        """
        Populate the possible reactions for this node using RetroBioCat dataset
        """
        prod = rdchiralReactants(self.smiles)
        for idx, name, rxn_smarts, rxn_type in self.retrobiocat.itertuples():
            rxn = rdchiralReaction(rxn_smarts)
            outcomes = rdchiralRun(rxn, prod, combine_enantiomers=False)
            for reaction in outcomes:
                reagents = [canonicalize_smiles(reagent) for reagent in reaction.split('.')]
                self.possible_reactions.append(Reaction(name, rxn_smarts, rxn_type, reagents))
        n = len(self.possible_reactions)
        self.weights = [1/n for _ in range(n)]

    def generate_reactions_rdenzyme(self):
        """
        Populate the possible reactions for this node using RdEnzyme
        """
        results, delta_scscore = ChemNode.analyzer.single_step_retro(self.smiles, max_precursors=10, debug=False, retrobiocat=True)
        self.weights = weight_normalizer(delta_scscore)
        for reaction in results:
            # Add abundance check here
            self.possible_reactions.append(Reaction(reaction[1], "", reaction[0], reaction[2].split('.')))
        
    def is_buyable(self) -> bool:
        """
        Check if this node is buyable
        """
        return self.buyable     

    def is_terminal(self):
        """
        Check if this node is a terminal node (solution or unmakeable)
        """
        return self.buyable or (not self.possible_reactions and not self.reactions)

    def is_fully_expanded(self):
        """
        Check if this node is fully expanded
        """
        return len(self.possible_reactions) == 0 and len(self.reactions) > 0

    def get_score(self) -> float:
        """
        Get the score of this node
        """
        if self.visits == 0:
            return 0
        return self.score / self.visits

    def get_MCTS_reaction(self) -> 'ReactNode':
        """
        Get the best reaction to expand based on MCTS selection function
        """
        if self.is_terminal():
            return None
        random.shuffle(self.reactions)
        return max(self.reactions, key = lambda x : x.get_mcts_value())
    
    def get_best_reaction(self) -> 'ChemNode':
        """
        Get the best child node
        """
        if self.is_terminal():
            return None
        random.shuffle(self.reactions) # Shuffle to avoid max first pick bias
        return max(self.reactions, key = lambda x : x.get_reaction_score())
    
    def get_random_reaction(self) -> Reaction:
        """
        Get a random reaction
        """
        if self.is_terminal():
            return None
        react = weighted_random(self.possible_reactions, self.weights)
        self.weights.remove(self.weights[self.possible_reactions.index(react)])
        self.possible_reactions.remove(react)
        return react