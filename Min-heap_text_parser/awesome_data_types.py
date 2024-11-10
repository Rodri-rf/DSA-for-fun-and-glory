# Let us start with the necessary classes for a k-arry tree implementation of
# a priority queue. It needs to satisfy the following properties:
# - 1: Min-heap property: the key of a node is less than or equal to the keys of its children
# - 2: As we are adding items to the tree from a an unsorted list, we want the closest ancestor 
# to each node to be the most recently added node with a smaller key, subject to the min-heap property
import re
from typing import TypeAlias
import pandas as pd
import igraph as ig
from igraph import Graph, EdgeSeq
import plotly.graph_objects as go
import time
import random

strain_name: TypeAlias = str
pcr_cycles: TypeAlias = str
segment: TypeAlias = str
variants_df: TypeAlias = pd.DataFrame
partitioned_variants: TypeAlias = dict[str, int]

class Node:
    def __init__(self, key, name):
        self.key = key
        self.name = name.strip()
        self.content = []
        self.children = []
        self.parent = None
        self.colour = 'white'
        self.id = random.randint(1, 1000000) # unique id for each node
    
    def __repr__(self) -> str:
        return f"Node({self.key}, '{self.name}')"

    def add_child(self, child):
        self.children.append(child)
        child.parent = self # make sure the child knows who its parent is

    def add_parent(self, parent):
        self.parent = parent

    def add_content(self, content):
        self.content.append(content)
    
    def __str__(self):
        return str(self.name)

    
class stronger_min_heap(Node):

    def __init__(self, key, name) -> None:
        # heap inherits from Node, that means it has all node attributes
        self.root = self
        self.current = self
        self.size = 1
        Node.__init__(self, key, name)

    def insert(self, node, curr):
        # we want to insert from the position of the last node added (curr)
        # we assume that the heap starts off staistfying properties 1 and 2.

        if self.root is None:
            self.root = node
            self.current = node
            self.size += 1
            return
        
        else:
            # case 1: key is more than the key of the current node, the current node has a new child!
            if node.key > curr.key:
                curr.add_child(node)
                self.current = node
                self.size += 1
                
            
            # case 2: key is exactly equal to the key of the current node, the current has a new sibling!
            elif node.key == curr.key:
                # what is a sibling if not a child of your parent?
                curr.parent.add_child(node)
                self.current = node
                self.size += 1
                
            
            # case 3: key is greater than the key of the current node, the current node has a new ancestor!
            elif node.key < curr.key:
                curr.parent.insert(node, curr.parent) 
        return 

    def get_node(self, name, key=None):
        # BFS
        name = name.strip()
        queue = []
        for child in self.children:
                queue.append(child)
        while queue:
            node = queue.pop(0)
            if node.name == name:
                return node
            for child in node.children:
                queue.append(child)
        raise ValueError("Node not found, perhaps you made a typo?")
    
    def show_tree(self):
        # print the tree in a nice format
        def print_tree(node, indent=""):
            print(indent + str(node))
            for child in node.children:
                print_tree(child, indent + "\t")
        print_tree(self.root)
        

    def regex_search(self, pattern):
        # given a regex pattern, return the first node whose name matches the pattern
        # BFS
        queue = []
        for child in self.children:
                queue.append(child)
        while queue:
            node = queue.pop(0)
            if re.match(pattern, node.name):
                return node
            for child in node.children:
                queue.append(child)
        raise ValueError("Node not found, perhaps you made a typo?")
    
    def nicest_visualization(self):
        # use igraph to set up a Tree and then plot it with plotly (traces)
        # To ensure the name of every vertex is unique, we will add the epoch time to the name of the vertex
        g = Graph(directed=True)
        def add_node_and_edges(graph, src_node):
            src_node.colour = 'grey'
            node_id = src_node.name + ', ' + str(src_node.key) + ', ' + f"node-id:{src_node.id}"
            graph.add_vertex(node_id)
            queue = []
            queue.append(src_node) # ENQUEUE
            while queue:
                curr = queue.pop(0) # DEQUEUE
                curr_name = curr.name + ', ' + str(curr.key) + ', ' + f"node-id:{curr.id}"
                for child in curr.children:
                    if child.colour == 'white': # if it has not been discovered yet
                        child.colour = 'grey'
                        name = child.name + ', ' + str(child.key) + ', ' + f"node-id:{child.id}"
                        graph.add_vertex(name)
                        graph.add_edge(curr_name, name)
                        queue.append(child) # ENQUEUE
                        # we dont really need to compute the distance beteen the nodes, we just need to add the nodes and edges
                curr.colour = 'black'

                    
        add_node_and_edges(g, self.root)
        # print number of edges and nodes in g 
        print(f"Number of edges: {g.ecount()}")
        print(f"Number of nodes: {g.vcount()}")
        fig = go.Figure()
        # Plot as a tree. Root is at the top, and children are below
        lay = g.layout_reingold_tilford(root=[0]) # layout determines the position of the nodes. 
        position = {k: lay[k] for k in range(len(lay.coords))} # position of each node
        Y = [lay[k][1] for k in range(len(lay.coords))] # y-coordinates of nodes
        M = max(Y)
        L = len(Y)
        X = [lay[k][0] for k in range(L)] # x-coordinates of nodes
        E = [e.tuple for e in g.es]
        V = [v['name'] for v in g.vs] # names of the vertices
        L = len(X)
        Y = [M - y for y in Y]
        labels = [g.vs[k]['name'] for k in range(L)] # 1 label per coordinate
        # Create Edges
        edge_x = []
        edge_y = []
        for edge in E: # E is a list of tuples, each tuple is an edge
            x0, y0 = X[edge[0]], Y[edge[0]]
            x1, y1 = X[edge[1]], Y[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
        # Create Nodes
        node_x = []
        node_y = []
        for i in range(L): # L is the number of nodes (strictly speaking, it is the number of coordinates in the layout)
            x, y = X[i], Y[i]
            node_x.append(x)
            node_y.append(y)
        fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode='lines', 
                                 line=dict(color='rgb(210,210,210)', width=3), hoverinfo='none',
                                 showlegend=False))
        fig.add_trace(go.Scatter(x=node_x, y=node_y, mode='markers', name='bla',
                                    marker=dict(symbol='circle-dot', size=30, color='rgb(50,220,210)', 
                                    line=dict(color='darkturquoise', width=1)), 
                                    text=labels, hoverinfo='text', opacity=0.8))
        # Fun fact: plotly parses strings as HTML, so we can use HTML tags to format the text. Show names of nodes
        fig.update_layout(title=f'Reingold-Tilford layout: <br><sup>Tree built from {self.root.name}</sup>', title_x=0.5, showlegend=False, 
                          hovermode='closest', margin=dict(b=20,l=5,r=5,t=40), height=700, xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                          yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
        fig.show()
        
        

class variant:
    # this is an attempt to nicely wrapp up all the potentially relevant information
    # about a variant in a genome. The conscious decision of going with an 
    # object-oriented approach, as opposed to something like a pandas dataframe,
    # is to make the code more extendable and dynamic.

    # We need to know,
    # - the postion of the variant (and read depth at that postion)
    # - the reference base at that position
    # - the alternative base(s) at that position, along with the frequency of each mutant base
    # - the type of variant (SNP, insertion, deletion, etc.)
    # - the effect of the variant (ORF, UTR, intron, differnt AA, etc.)

    def __init__(self, position, read_depth , ref_base, alt_bases:dict, insertions=None, deletions=None, variant_type=None, effect=None, total_freq=None, info=None, segment=None, pcr_cycles=None):
        self.position = position
        self.read_depth = read_depth
        self.ref_base = ref_base
        self.alt_bases = alt_bases
        self.variant_type = variant_type # SNP, insertion, deletion
        self.effect = effect
        self.total_freq = total_freq
        self.insertions = insertions
        self.deletions = deletions
        self.info = info # any other information that might be relevant. e.g. gene name, strain, etc.
        self.segment = segment # viral DNA is divided into segments, this is the segment the variant belongs to
        self.pcr = pcr_cycles

    def __repr__(self):
        return f"Variant at position {self.position}, with reference base {self.ref_base} and alternative bases {self.alt_bases}"
    

if __name__ == '__main__':
    # Test implementation
    # create a tree
    root = stronger_min_heap(key=0, name="root")
    root.current = root # This is kinda awkward, but it is necessary for the insert method to work

    from random import randint
    for i in range(60):
        curr = root.current # this is the latest node added
        node = stronger_min_heap(key=randint(1,4), name=f"node_{i%3}")
        root.current.insert(node, root.current)
        node.root = curr
        root.current = node

    root.show_tree()
    root.nicest_visualization()



