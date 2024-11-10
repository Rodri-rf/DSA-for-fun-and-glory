# Description: Implementation of an interval tree for genomics applications
from typing import Any
import random
from igraph import Graph, EdgeSeq
import plotly.graph_objects as go

interval_tree_genomics = Any

class interval_tree_genomics:

    id_counter = 0  # Class-level counter for unique IDs

    # Implementation as augmented AVL (self-balancing) tree
    # Augmented properties:
    # - max_position: the maximum position of the interval in the subtree rooted at the node
    # - balance_factor: the difference between the height of the left and right subtrees
    # The goal is to serve as a memoization table for sequence alignment algorithms.
    # It anwswers the question: "Have we alread alligned part of this sequence before?"
    # sources in case I loose my sanity: https://cs.stackexchange.com/questions/48861/balance-factor-changes-after-local-rotations-in-avl-tree
    
    def __init__(self, start, end, left: interval_tree_genomics=None, right: interval_tree_genomics=None, parent: interval_tree_genomics=None, sequence=None, ref_seq_id=None, targ_seq_id=None):
        self.node_id = interval_tree_genomics._id_counter
        interval_tree_genomics._id_counter += 1
        self.left_child = left
        self.right_child = right
        self.parent = parent
        self.sequence = sequence # DNA sequence 
        self.balance_factor = 0
        self.interval_start = start # let's use the start of the interval as the key
        self.interval_end = end
        self.ref_seq_id = ref_seq_id
        self.targ_seq_id = targ_seq_id
        self.max_position = end # this is the augmented property

    def insert(self, node: interval_tree_genomics):

        '''
        AVL insertion mantra:
            - Let n be the new node we just inserted
            - Let p the n's parent
            - Let m be the closest ancestor of p that is *not* balanced (i.e., has a balance factor != 0)
        then, only the balance factors between m and p (inclusive) will change after insertion
        '''

        # first insert as in a regular BST

        if self.interval_start >= node.interval_start:
            # insert at left subtree
            if self.left_child is not None:
                self.left_child.insert(node)
            else:
                node.parent = self
                self.left_child = node
                self._update_balance_factors_and_rebalance(insertion_at_right_subtree=False)

        elif self.interval_start < node.interval_start:
            # insert at right subtree
            if self.right_child is not None:
                self.right_child.insert(node)
            else:
                node.parent = self
                self.right_child = node
                self._update_balance_factors_and_rebalance(insertion_at_right_subtree=True)
        if self.left_child is None and self.right_child is None:
            return
        if self.right_child is None:
            self.max_position = max(self.max_position, self.left_child.max_position)
        elif self.left_child is None:
            self.max_position = max(self.max_position, self.right_child.max_position)
        else:
            self.max_position = max(self.left_child.max_position, self.right_child.max_position, self.max_position)

    
    def get_root(self):
        current = self
        while current.parent is not None:
            current = current.parent
        return current

    def update_root(self):
        root = self.get_root()
        if root.parent is not None:
            root.parent = None
        return root

    def _single_left_rotation(self):
        print("Single left rotation")
        new_root = self.right_child
        self.right_child = new_root.left_child
        if new_root.left_child is not None:
            new_root.left_child.parent = self
        new_root.parent = self.parent
        if self.parent is None:
            # This node is the root, so we need to update the root reference outside this method
            pass
        else:
            if self == self.parent.left_child:
                self.parent.left_child = new_root
            else:
                self.parent.right_child = new_root
        new_root.left_child = self
        self.parent = new_root

        self.balance_factor = 0
        new_root.balance_factor = 0

        return new_root.update_root()

    def _single_right_rotation(self):
        print("Single right rotation")
        new_root = self.left_child
        self.left_child = new_root.right_child
        if new_root.right_child is not None:
            new_root.right_child.parent = self
        new_root.parent = self.parent
        if self.parent is None:
            # This node is the root, so we need to update the root reference outside this method
            pass
        else:
            if self == self.parent.left_child:
                self.parent.left_child = new_root
            else:
                self.parent.right_child = new_root
        new_root.right_child = self
        self.parent = new_root

        self.balance_factor = 0
        new_root.balance_factor = 0

        return new_root.update_root()
    
    def _double_right_left_rotation(self):
        self.right_child._single_right_rotation()
        self._single_left_rotation()
    
    def _double_left_right_rotation(self):
        self.left_child._single_left_rotation()
        self._single_right_rotation()


    def _update_balance_factors_and_rebalance(self, insertion_at_right_subtree):
        if insertion_at_right_subtree:
            self.balance_factor +=1
        else:
            self.balance_factor -= 1
        
        if self.balance_factor == 0:
            # we gucci
            return
        
        elif abs(self.balance_factor) > 2:
            raise ValueError("AAAAAA ILLEGAL BALANCE FACTOR")
        
        elif self.balance_factor == 2 and self.right_child.balance_factor == 1:
            self._single_left_rotation()
            return
            

        elif self.balance_factor == 2 and self.right_child.balance_factor == -1:
            self._double_right_left_rotation()
            return
            

        elif self.balance_factor == -2 and self.left_child.balance_factor == -1:
            self._single_right_rotation()
            return
                
        
        elif self.balance_factor == -2 and self.left_child.balance_factor == 1:
            self._double_left_right_rotation()
            return
        
        
        if self.parent is not None:
            insertion_at_right_subtree = self.parent.right_child.node_id == self.node_id
            self.parent._update_balance_factors_and_rebalance(insertion_at_right_subtree=insertion_at_right_subtree)
        
    def delete(self, node):
        pass
        # probably not super important for now

    def __repr__(self):
        return f"({self.interval_start}, {self.interval_end})"
    
    def visualize(self):
        g = Graph(directed=True)
        g.vs["name"] = []  # Initialize the 'name' attribute for vertices

        def add_nodes_and_edges(graph, src_node):
            if src_node is None:
                return
            
            node_label = f"{src_node.interval_start}-{src_node.interval_end}\nBF={src_node.balance_factor}\nMax={src_node.max_position}"
            if node_label not in graph.vs["name"]:
                graph.add_vertex(name=node_label, dummy=False)

            if src_node.left_child is not None:
                child_label = f"{src_node.left_child.interval_start}-{src_node.left_child.interval_end}\nBF={src_node.left_child.balance_factor}\nMax={src_node.left_child.max_position}"
                if child_label not in graph.vs["name"]:
                    graph.add_vertex(name=child_label, dummy=False)
                graph.add_edge(node_label, child_label)
                add_nodes_and_edges(graph, src_node.left_child)
            else:
                # Add a dummy left child to ensure correct visualization
                dummy_label = f"dummy_left_{node_label}"
                graph.add_vertex(name=dummy_label, dummy=True)
                graph.add_edge(node_label, dummy_label)

            if src_node.right_child is not None:
                child_label = f"{src_node.right_child.interval_start}-{src_node.right_child.interval_end}\nBF={src_node.right_child.balance_factor}\nMax={src_node.right_child.max_position}"
                if child_label not in graph.vs["name"]:
                    graph.add_vertex(name=child_label, dummy=False)
                graph.add_edge(node_label, child_label)
                add_nodes_and_edges(graph, src_node.right_child)
            else:
                # Add a dummy right child to ensure correct visualization
                dummy_label = f"dummy_right_{node_label}"
                graph.add_vertex(name=dummy_label, dummy=True)
                graph.add_edge(node_label, dummy_label)

        add_nodes_and_edges(g, self)

        layout = g.layout_reingold_tilford()
        Y = [layout[k][1] for k in range(len(layout.coords))]  # y-coordinates of nodes
        M = max(Y)
        X = [layout[k][0] for k in range(len(layout.coords))]  # x-coordinates of nodes
        E = [e.tuple for e in g.es]
        labels = [g.vs[k]["name"] for k in range(len(g.vs))]  # names of the vertices

        # Create Edges
        edge_x = []
        edge_y = []
        for edge in E:
            if not g.vs[edge[1]]["dummy"]:  # Skip edges going to dummy nodes
                x0, y0 = X[edge[0]], Y[edge[0]]
                x1, y1 = X[edge[1]], Y[edge[1]]
                edge_x.extend([x0, x1, None])
                edge_y.extend([y0, y1, None])

        # Create Nodes
        node_x = []
        node_y = []
        node_labels = []
        for i in range(len(X)):
            if not g.vs[i]["dummy"]:
                x, y = X[i], Y[i]
                node_x.append(x)
                node_y.append(y)
                node_labels.append(labels[i])

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode='lines', 
                                line=dict(color='rgb(210,210,210)', width=3), hoverinfo='none',
                                showlegend=False))
        fig.add_trace(go.Scatter(x=node_x, y=node_y, mode='markers', name='Nodes',
                                marker=dict(symbol='circle-dot', size=30, color='rgb(50,220,210)', 
                                line=dict(color='darkturquoise', width=1)), 
                                text=node_labels, hoverinfo='text', opacity=0.8))

        fig.update_layout(title='AVL Tree Visualization', showlegend=False, 
                        xaxis=dict(showline=False, zeroline=False), 
                        yaxis=dict(showline=False, zeroline=False), 
                        plot_bgcolor='white')
        
        # Invert the y-axis to ensure the tree is not upside down
        fig.update_yaxes(autorange="reversed")
        
        fig.show()


        
 # --------------------------------    

if __name__ == "__main__":
    # lets test the AVL tree
    # we will insert intervals covering from 0 to 1000:
    prev = None

    for i in range(50):
        start = random.randint(0, 1000)
        end = random.randint(start, 1000)
        print(f"Inserting interval {start} to {end}")
        node = interval_tree_genomics(start, end)
        if i == 0:
            tree = node
        else:
            tree.insert(node)
    tree.visualize()
        

