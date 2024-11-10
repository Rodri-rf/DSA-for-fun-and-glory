import parsing_and_tokenizing.awesome_data_types as adt
from parsing_and_tokenizing.latex_parser_rodrigos_version import *
from analysis_tools import parse_content_variants, partition_variants_count
from analysis_tools import *
from os import listdir
import re
import statistics as stats
import pandas as pd
from tqdm import tqdm
from typing import TypeAlias
import matplotlib.pyplot as plt
import importlib

# Useful background info: https://www.nature.com/articles/s41564-019-0513-7
strain_name: TypeAlias = str
pcr_cycles: TypeAlias = str
segment: TypeAlias = str
variants_df: TypeAlias = pd.DataFrame
partitioned_variants: TypeAlias = dict[str, int]
avg_read_depth: TypeAlias = float
variant_objects: TypeAlias = list[adt.variant]
godzilla_dict: TypeAlias = dict[strain_name: dict[pcr_cycles: dict[segment: list[variants_df, partitioned_variants, avg_read_depth, variant_objects]]]] 

def construct_summary_dict(strains: list[str], plot_tree=False, add_all_per_position=True) -> godzilla_dict:
    '''
    This function constructs a dictionary that summarizes the information contained in the latex files of the strains.
    It starts by creating a tree structure of the latex file, and then it extracts the variant analysis information.

    The dictionary is structured as follows:
    - The keys are the names of the strains.
    - The values are dictionaries with the following structure:
        - The keys are the number of PCR cycles.
        - The values are dictionaries with the following structure:
            - The keys are the segments of the influenza virus.
            - The values are lists with the following structure:
                - The first element is a pandas dataframe with the variants information.
                - The second element is a dictionary with the partitioned variants count.
                - The third element is the average read depth.
                - The fourth element is a list of variant objects.
    '''
    dic_summary = {}
    for strain in tqdm(strains, desc="Processing strains..."):
        strain_name = strain.split("\\")[1]
        files = [strain + "\\" + f for f in listdir(strain) if f.endswith(".tex")] # get all the .tex files in the directory, ignore the rest.
        dic_summary[strain_name] = {}

        for file_path in sorted(files, key=lambda f: re.search(r".*(_[1-9]{2}_).*", f).group(1)[1:-1]): # sort by number of PCR cycles
            pcr_cycles = re.search(r".*(_[1-9]{2}_).*", file_path).group(1)[1:-1]
            file_tree = create_tex_tree(file_path)
            if plot_tree:
                file_tree.show_tree()
                file_tree.nicest_visualization() # shows the tree structure of the latex file 

            variant_analysis = file_tree.get_node("Variant Analysis") # get the node that contains the variant analysis information
            na = variant_analysis.regex_search("NA.*")
            variants_na, variant_objects_na  = parse_content_variants(mess=na.content, name=pcr_cycles, segment="NA", strain=strain_name)
            # set the pcr_cycles attribute of the variant objects to the number of PCR cycles
            for var in variant_objects_na:
                var.pcr_cycles = pcr_cycles 
            ha = variant_analysis.regex_search("HA.*")
            variants_ha, variant_objects_ha = parse_content_variants(mess=ha.content, name=pcr_cycles, segment="HA", strain=strain_name)
            for var in variant_objects_ha:
                var.pcr_cycles = pcr_cycles

            dic_summary[strain_name][pcr_cycles] = {"NA": [variants_na, partition_variants_count(variants=variant_objects_na, partitions={'0.5-1':0, '1-2':0, '2-5':0, '5-10':0, '10+':0}, add_all_per_position=add_all_per_position), 
                                                           stats.mean(variants_na["depth"]), variant_objects_na],
                                                    "HA": [variants_ha, partition_variants_count(variants=variant_objects_ha, partitions={'0.5-1':0, '1-2':0, '2-5':0, '5-10':0, '10+':0}, add_all_per_position=add_all_per_position), 
                                                           stats.mean(variants_ha["depth"]), variant_objects_ha]} # Godzilla dictionary
    return dic_summary

def construct_summary_df(summary_dict: godzilla_dict, output_dir: str, summary_df: pd.DataFrame, export_csv=True) -> variants_df:
    for strain, info in summary_dict.items():
        for pcr_cycles, segments in info.items():
            for segment, variants in segments.items():
                # add revelant information to the summary dataframe as a new row.
                tot_indel_freq = 0
                total_insertion_freq = 0
                total_deletion_freq = 0
                tot_indels = 0
                for var in variants[3]:
                    #for a given batch, the total indel frequency is the sum of the frequency of deletions and the frequency of insertions for every variant in the batch
                    total_insertion_freq += sum([val / var.read_depth for key, val in var.insertions.items()])
                    total_deletion_freq += sum([val / var.read_depth for key, val in var.deletions.items()])
                    tot_indel_freq =+ total_insertion_freq + total_deletion_freq
                    tot_indels += sum([val for key, val in var.insertions.items()])
                    tot_indels += sum([val for key, val in var.deletions.items()])

                summary_df = summary_df._append({'Influenza Strain': strain, 'PCR cycles': pcr_cycles, 'Segment': segment, 'Average Read Depth': variants[2], 
                                                 'Total Indels': tot_indels, 'Total Indel Frequency' :tot_indel_freq,
                                                'Variants:0.5%-1%': variants[1]['0.5-1'], 'Variants: 1%-2%': variants[1]['1-2'], 'Variants: 2%-5%': variants[1]['2-5'],
                                                'Variants: 5%-10%': variants[1]['5-10'], 'Variants: 10%+': variants[1]['10+']}, ignore_index=True)
    if export_csv:
        summary_df.to_csv(output_dir + "\\" + "summary.csv", index=False)
    return summary_df




    
    





                
