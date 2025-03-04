from ChemNode import ChemNode
from typing import List
from itertools import product
from typing import List

# def prune_tree(node: ChemNode):
#     """
#     Prune the tree by removing all non-solution pathways.

#     :param node: The root ChemNode to prune.
#     :return: None (modifies the tree in-place).
#     """
#     # Iterate over a copy of the reactions to safely modify the original list
#     for reaction in node.reactions[:]:
#         # Check if all reagents in the reaction are solutions
#         if all(reagent.solution for reagent in reaction.precursors):
#             # Recursively prune each reagent
#             for reagent in reaction.precursors:
#                 prune_tree(reagent)
#         else:
#             # Remove non-solution reagents from the reaction
#             for reagent in reaction.precursors[:]:
#                 if not reagent.solution:
#                     reaction.precursors.remove(reagent)
            
#             # If the reaction has no reagents left, remove the reaction
#             if not reaction.precursors:
#                 node.reactions.remove(reaction)


# # Go down tree, copying along a given path to add to a single path list
# def generate_subtrees(node: ChemNode) -> ChemNode:
#     """
#     Generate a solution subtree from a given node.

#     :param node: The root ChemNode to generate a subtree from.
#     :return: The root of the generated subtree.
#     """
#     # Base case: if the node is a leaf
#     if node.solution and not node.reactions: # Leaf node
#         node.solution = False
#         return node.copy()
    
#     copy = node.copy()
#     # If the node is not a leaf, generate subtrees for all reactions
       
#     reaction = node.reactions[0]
#     react = reaction.copy()
#     for reagent in reaction.precursors:
#         if reagent.solution:
#             react.add_reagent(generate_subtrees(reagent))

#     if all(not reagent.solution for reagent in node.reactions[0].precursors):
#         node.reactions.pop(0)

#     if not node.reactions:
#         node.solution = False

    

#     copy.reactions.append(react)
#     return copy

# def generate_paths1(root: ChemNode) -> List[ChemNode]:
#     """
#     Generate all solution paths from a given node.

#     :param root: The root ChemNode to generate paths from.
#     :return: A list of all solution paths.
#     """
#     paths = []
#     while root.solution:
#         paths.append(generate_subtrees(root))
#     return paths

def path_explorer(root: ChemNode):
    """
    Prints the path given a subtree

    :param root: The root of the solution subtree.
    """
    stack = []
    stack.append(root)
    count = 0
    while stack:
        print(f"BRANCH {count}:\n-----------------")
        count += 1
        root = stack.pop()
        print(f'Chem: {root.smiles}')
        while root.reactions:
            print(f'Reaction: {root.reactions[0].reaction_name}')
            if len(root.reactions[0].precursors) > 1:
                print(f"SPLIT: Chem1: {root.reactions[0].precursors[0].smiles}, Chem2: {root.reactions[0].precursors[1].smiles}\n")
                stack.append(root.reactions[0].precursors[1])
            else:
                print(f'Chem: {root.reactions[0].precursors[0].smiles}')
            root = root.reactions[0].precursors[0]
        print("BUYABLE\n")

def generate_paths(node: ChemNode) -> List[ChemNode]:
    """
    Recursively generate all solution paths from a given ChemNode.
    Each ChemNode in the returned tree will have exactly one child reaction,
    and each reaction can have multiple precursor ChemNodes.
    """
    # Base case: if the node is a leaf solution (flagged and no reactions)
    if node.solution and not node.reactions:
        leaf = node.copy()
        # leaf.solution = False  # Optionally mark as processed
        return [leaf]
    
    all_paths = []
    # Iterate over each reaction option (OR branch) at the current node.
    for reaction in node.reactions:
        precursor_solution_lists = []
        valid_reaction = True
        # For an AND reaction, all precursors must provide at least one solution.
        for precursor in reaction.precursors:
            # If the precursor isn't flagged or doesn't lead to a full pathway, skip this reaction.
            if not precursor.solution:
                valid_reaction = False
                break
            sub_paths = generate_paths(precursor)
            if not sub_paths:
                valid_reaction = False
                break
            precursor_solution_lists.append(sub_paths)
        
        # Only continue if every precursor returned at least one solution path.
        if not valid_reaction:
            continue
        
        # Combine precursor solutions: each combination yields a complete branch.
        for combo in product(*precursor_solution_lists):
            node_copy = node.copy()
            reaction_copy = reaction.copy()
            # Attach each precursor solution from the combination as a reagent to the reaction.
            for precursor_path in combo:
                reaction_copy.add_reagent(precursor_path)
            # Now the ChemNode gets a single child reaction (the copied reaction)
            node_copy.reactions = [reaction_copy]
            all_paths.append(node_copy)
    
    # If no reaction yielded a full pathway, mark this node as not part of a solution.
    if not all_paths:
        # node.solution = False
        pass
    return all_paths