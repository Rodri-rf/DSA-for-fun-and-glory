import re
import parsing_and_tokenizing.awesome_data_types as adt
import plotly.graph_objects as go
import networkx as nx
import pandas as pd
from typing import Union


class node:
    def __init__(self, name):
        self.name = name
        self.children = []
        self.parents = []
        self.content = []
        self.visited = False
    
    def add_child(self, child):
        self.children.append(child)

    def add_parent(self, parent):
        self.parents.append(parent)
    
    def add_content(self, content):
        self.content.append(content)
    
    def __str__(self):
        return self.name
    
def get_main_sections(filename):
    with open(filename) as f:
        lines = f.readlines()
    head = node("head")
    current = head
    for line in lines:
        if re.match(r"\\title", line):
            head.add_content(line)

        elif re.match(r"\\section", line):
            # the name of each node is the title of the section. 
            # We format it to remove the \section command and the brackets
            new_node = node(re.sub(r"\\(sub)*section\{([\w\s]+?)\}(.*)", r"\2", line))
            new_node.add_parent(current)
            current.add_child(new_node)

        elif current.children != []: 
            current.children[-1].add_content(line)
        
        if re.match(r"\\end{document}", line):
            break
    return head

def create_tex_tree(filename):
   # Create a tree of the latex file, with the following properties:
    # - Each node is a section of the latex file.
    # - Each subsection of a given section is a child of that section.
    # - Each node has a list of the content of the section.
    # - the tree implements a min heap, where the keys represent whether a given node
    # is a section, subsection, subsubsection, etc.
    with open(filename) as f:
        lines = f.readlines()
    head = adt.stronger_min_heap(key=0, name="head")
    head.current = head
    head.root = head
    for line in lines:
        if re.match(r"\\title", line):
            head.add_content(line)
            head.name = re.sub(r"\\title\{(.*)\}", r"\1", line)
            
        elif re.match(r"\\(sub)*section", line):
            # the name of each node is the title of the section. 
            # We format it to remove the \section command and the brackets
            key = len(re.findall(r"sub", line)) + 1
            curr = head.current # this is the latest node added
            new_node = adt.stronger_min_heap(key=key, 
                                             name=re.sub(r"\\(sub)*section\{([\w\s\W]+?)\}(.*)", r"\2", line))
            head.current.insert(new_node, curr)
            head.size += 1
            new_node.root = curr
            head.current = new_node
        else:
            head.current.add_content(line)

        if re.match(r"\\end{document}", line):
            break
    print(f"Tree created successfully!    size: {head.size}    root: {head.root}    current: {head.current}")
    return head

def graph_tex_tree(head):
    # Create a directed graph from the tree data. Note that each element in tree_data is a tuple (parent, current_node)
    G = nx.DiGraph()
    for node in head.children:
        G.add_edge(head.name, node)
    for child in head.children:
        for node in child.content:
            G.add_edge(child.name, node)
        for grandchild in child.children:
            for node in grandchild.content:
                G.add_edge(grandchild.name, node)
    return G

def show_graph(G):
    # Create a plotly figure to display the graph
    pos = nx.spring_layout(G)
    edge_x = []
    edge_y = []

    # Create edges
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
    
    # Create nodes
    node_x = []
    node_y = []
    node_text = []
    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_text.append(node)

    # Create figure
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode='lines', line=dict(width=0.5, color='black')))
    fig.add_trace(go.Scatter(x=node_x, y=node_y, mode='markers', marker=dict(size=12, color='blue')))
    fig.update_layout(showlegend=False)
    fig.show()

# ----------------------------------- Genetic Stability Specific -----------------------------------------

def parse_content_variants(mess, name=None, segment=None, strain=None) -> Union[pd.DataFrame, list]:
    # example:
    # - pos:711 depth:421401 ref:G A:1.56454% C:0.00308495% G:98.4302% T:0.00213573% Indels:-A:6, +A:7 { ORF:NA ref:Gly Gly(=):98.399% Glu(6=):1.58161% Other:0.01939% }
    variants = pd.DataFrame(columns=["position", "depth", "ref", "SNPs_variant", "Insertions", "Deletions", "total freq", "other info", 'name'])
    variant_objects = []
    
    for line in mess:

        if not line.startswith("pos:"):
            continue # skip the line if it does not contain the variant information

        alternative_bases = {"A" : 0, "C" : 0, "G" : 0, "T" : 0}
        other_info = []
        variant_types = set()
        insertions = {}
        deletions = {}
        flag = False

        for word in line.split():
            if word.startswith("pos:"):
                position = int(word.split(":")[1])
            elif word.startswith("depth:"):
                depth = int(word.split(":")[1])
            elif word.startswith("ref:") and not flag:
                ref = word.split(":")[1]
            elif re.match(r"[A, G, C, T]:(.)", word):
                base = word[0]
                # for the frequency we need to match all ASCII characters from ":" to the "%" sign
                alternative_bases[base] = float(re.search(r":(.*?)(?=\\)", word).group(1))

            elif word.startswith('Indels') or "+" in word or "-" in word:
                # + means insertions, - means deletions
                if '+' in word:
                    match = re.match(r".*\+([A, G, C, T]+):(\d+)", word)
                    insertion = match.group(1)
                    occurrences = int(match.group(2))
                    insertions[insertion] = occurrences
                    variant_types.add("Insertion")
                elif '-' in word:
                    match = re.match(r".*\-([A, G, C, T]+):(\d+)", word)
                    deletion = match.group(1)
                    occurrences = int(match.group(2))
                    deletions[deletion] = occurrences
                    variant_types.add("Deletion")

            elif flag:
                other_info.append(word)

            elif word == "\{":
                # this is the part of the string that contains the effect of the variant (ORFs, etc.)
                flag = True
            elif word == "}/":
                flag = False

        if any(value != 0 for value in alternative_bases.values()): # we define replacements (SNPs) as positions with an alternative base
            variant_types.add("SNP")
            snps_freq = sum([value for key, value in alternative_bases.items() if key != ref]) 
            insertions_freq = sum([value / depth for key, value in insertions.items()]) # implicit representation invariant assumed: depth is not zero
            deletions_freq = sum([value / depth for key, value in deletions.items()]) # implicit representation invariant assumed: depth is not zero
            total_freq = snps_freq + insertions_freq + deletions_freq # total frequency of variants at a given position is the sum of the frequencies of SNPs, insertions and deletions

        variants = variants._append({"position": position, "depth": depth, "ref": ref, "SNPs_variant": alternative_bases,
                                        "Insertions": insertions, "Deletions": deletions, "total freq": total_freq, 
                                        "other info": other_info, 'name': name}, ignore_index=True)
        var = adt.variant(position=position, read_depth=depth, 
                          ref_base=ref, # the base that is present in the reference genome at this position
                          alt_bases=alternative_bases, # dictionary where the keys are the alternative bases and the values are the frequencies of those bases
                          variant_type=variant_types, # set of the types of variants present at this position (SNP, Insertion, Deletion, etc.)
                          total_freq=total_freq, # the total frequency of all the variants at this position (calculated above)
                          effect=other_info, # the effect of the variant (ORF, etc.)
                          segment=segment, # the segment of the influenza virus to which the variant belongs (NA, HA, etc.)
                          info=strain, # the strain of the influenza virus to which the variant belongs
                          insertions=insertions, # dictionary where the keys are the inserted bases and the values are the number of occurrences of that base
                          deletions=deletions) # dictionary where the keys are the deleted bases and the values are the number of occurrences of that base
        variant_objects.append(var)
    return variants, variant_objects # return a pandas dataframe and a list of variant objects              





        