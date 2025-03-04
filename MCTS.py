from ReactNode import ReactNode
from ChemNode import ChemNode
from utils import *
from typing import List, Deque
from collections import deque
from multiprocessing import Process, Queue, Pool
from tqdm import tqdm
import sqlite3
import pandas as pd
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from rdenzyme import rdchiralRun

PAYOUT_DEPTH_LIMIT = 0
PAYOUT_DEAD_END = -0.01
PAYOUT_BUYABLE = 1

# Monte Carlo Tree Search
class MCTS:
    def __init__(self, root:ChemNode, max_depth:int = 8, exploration_param:float = 1.414, max_processes:int = 12):
        self.root = root
        self.simulate_stack:Deque[str] = deque()
        self.max_processes = max_processes
        self.max_depth = max_depth
        self.exploration_param = exploration_param
        self.buyables:sqlite3.Cursor = ChemNode.buyables
        self.abundants:sqlite3.Cursor = ChemNode.abundants
        self.template_set:pd.DataFrame = ChemNode.retrobiocat
        self.analyzer:Retrosim = ChemNode.analyzer


    # def select(self, root:ChemNode) -> 'ReactNode':
    #     temp = root
    #     while not temp.is_terminal() and temp.is_fully_expanded():
    #         temp = temp.get_MCTS_reaction()
    #     if not temp.is_fully_expanded():
    #         return self.expand(temp)
    #     elif temp.is_terminal():
    #         print("Terminal node reached")
    #         return None
    #     return temp     

    def multi_select(self, root:ChemNode):
        selected:List[ChemNode] = []
        select_stack:Deque[ChemNode] = deque()
        select_stack.append(root)

        while True:
            while select_stack:
                temp = select_stack.pop()

                if temp.is_terminal():
                    if temp.is_buyable():
                        if temp.visits == 0:
                            selected.append(temp)
                    else:
                        selected.append(temp)
                    continue
                elif not temp.is_fully_expanded():
                    expanded = self.expand(temp)
                    if not expanded:
                        print("Nothing expanded")
                    for child in expanded:
                        selected.append(child)
                    continue
                react = temp.get_MCTS_reaction()
                for precursor in react.precursors:
                    select_stack.append(precursor)

            if selected:
                return selected
            else:
                rand = random.choice(root.reactions)
                for precursor in rand.precursors:
                    select_stack.append(precursor)

    def expand(self, node:ChemNode) -> List[ChemNode]:
        """
        Expand a random child chem node 
        """        
        reaction = node.get_random_reaction()
        if reaction is None:
            return []
        react = ReactNode(reaction.name, reaction.smarts, reaction.react_type, node, 0.95)
        for reagent in reaction.precursors:
            react.add_reagent(ChemNode(reagent, node.depth + 1, react))
        node.reactions.append(react)
        return react.precursors
    

    def simulate(self, node:ChemNode):
        """
        Simulation/rollout of the node
        """
        smile = node.smiles
        depth = node.depth
        while True:
            if depth >= self.max_depth:
                return PAYOUT_DEPTH_LIMIT # Too far
            if check_buyable(smile, self.buyables):
                return PAYOUT_BUYABLE # buyable/abundant
            reactants = self.generate_random_reactions_rdenzyme(smile)
            if reactants == None:
                return PAYOUT_DEAD_END # Unmakeable (dead end)
            smile = complete_random(reactants)
            depth += 1
    

    def simulate_handler(self, node:ChemNode):
        smile = node.smiles
        depth = node.depth
        num_processes = 0
        self.simulate_stack.append(smile)
        
        while self.simulate_stack:
            if num_processes < self.max_processes:
                smile = self.simulate_stack.pop()
                p = Process(target=self.multi_simulate, args=(smile, depth))
                p.start()
                num_processes += 1
            else:
                p.join()
                num_processes -= 1


    def multi_simulate(self, smile:str, depth:int):
        while True:
            if depth >= self.max_depth:
                return PAYOUT_DEPTH_LIMIT # Too far
            if check_buyable(smile, self.cursor):
                return PAYOUT_BUYABLE # buyable
            reactants = self.generate_random_reactions_rdenzyme(smile)
            if reactants == None:
                return PAYOUT_DEAD_END # Unmakeable (dead end)
            smile = reactants[0]
            if len(reactants) > 1:
                for reagent in reactants[1:]:
                    self.simulate_stack.append(reagent) # Add all additional reagents to the stack
            depth += 1
        
    def backpropagate(self, node:ChemNode, reward:float):
        node.visits += 1
        node.score += reward
        if node.parent_reaction and node.parent_reaction.parent_chemical:
            if all([precursor.solution for precursor in node.parent_reaction.precursors]):
                node.parent_reaction.parent_chemical.solution = True
            self.backpropagate(node.parent_reaction.parent_chemical, reward * node.parent_reaction.weight)
        
    def MCTS(self, iterations:int = 1000):
        """
        MCTS algorithm
        """
        for _ in tqdm(range(iterations)):
            selected = self.multi_select(self.root)
            if not selected:
                continue
            for node in selected:
                reward = self.simulate(node)
                self.backpropagate(node, reward)
    

    def generate_random_reactions_retrobiocat(self, smile:str) -> List:
        """
        Get a random reaction from the retrobiocat database
        """
        prod = rdchiralReactants(smile)
        children = []
        for idx, name, rxn_smarts, rxn_type in self.template_set.itertuples():
            rxn = rdchiralReaction(rxn_smarts)
            outcomes = rdchiralRun(rxn, prod, combine_enantiomers=False)
            for reagents in outcomes:
                children.append(reagents.split('.'))

        return complete_random(children)
    
    def generate_random_reactions_rdenzyme(self, smile:str) -> str:
        """
        Get a random reaction from the rdenzyme database
        """
        children = []
        results, delta_scscore = ChemNode.analyzer.single_step_retro(smile, max_precursors=10, debug=False, retrobiocat=False)
        for reaction in results:
            for reagent in reaction[2].split('.'):
                children.append(reagent)

        if len(children) == 0:
            return None
        
        return children