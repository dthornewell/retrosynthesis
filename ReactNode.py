from typing import List
import math

EXPLORATION_PARAM = 1.414

class ReactNode:
    """
    Basic class for a reaction node in the tree
    Allows reactions to have multiple reagents
    """
    def __init__(self, reaction_name:str, smarts:str, react_type:str, parent_chemical:'ChemNode', weight:float):
        """
        Initialize a ReactNode.

        :param reaction_name: Name of the reaction.
        :param smarts: SMARTS string for the reaction.
        :param react_type: Type of the reaction.
        :param parent_chemical: Parent chemical node.
        :param weight: Weight of the reaction.
        """
        self.reaction_name = reaction_name
        self.smarts = smarts
        self.react_type = react_type 
        self.parent_chemical = parent_chemical
        self.precursors:List['ChemNode'] = []  # Will hold ChemicalNodes representing precursors
        self.weight = weight

    def copy(self) -> 'ReactNode':
        """
        Copy a ReactNode
        """
        react = ReactNode(self.reaction_name, self.smarts, self.react_type, self.parent_chemical, self.weight)
        react.precursors = []
        return react

    def add_reagent(self, chemical_node:'ChemNode'):
        self.precursors.append(chemical_node)
    
    def get_mcts_value(self) -> float:
        """
        UCT value of this node
        """
        visits = 0
        score = 0
        # for reagent in self.precursors:
        #     visits += reagent.visits
        #     score += reagent.score
        
        solution = 0
        if self.parent_chemical.solution:
            solution = 0.1 * self.parent_chemical.depth # Penalize solution nodes

        for reagent in self.precursors:
            score += reagent.get_score()
            visits += reagent.visits
        if visits == 0:
            return float('inf')
        uct = EXPLORATION_PARAM * math.sqrt(math.log(self.parent_chemical.visits) / visits)
        return score + uct - solution
    


    def get_reaction_score(self):
        """
        Reaction score based on reagent successes
        """
        # visits = 0
        # score = 0
        # for reagent in self.precursors:
        #     visits += reagent.visits
        #     score += reagent.score
        # if visits == 0:
        #     return float('inf')
        # return score / visits

        return sum([precursor.get_score() for precursor in self.precursors]) / len(self.precursors)

# For type hinting
from ChemNode import ChemNode
