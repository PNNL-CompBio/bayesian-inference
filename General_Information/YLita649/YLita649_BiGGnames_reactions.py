#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 15:08:56 2020

@author: cies677
"""

#########################################################################
# DON'T RUN THIS CODE ALL AT ONCE. RUN THE FIRST FIVE BLOCKS, FOLLOW THE
# INSTRUCTIONS, THEN PROCEED WITH THE REST OF THE CODE.
#
# Set your working directory to the `bayesian-inference' folder or update
# the OMICS_FILE and JSON_MODEL pathnames.
#########################################################################


#################################
# YLita649_BiGGnames_reactions.py
#
# Molly Mersinger (EBSD)
# molly.mersinger@pnnl.gov
#
# Danielle Ciesielski (EBSD)
# danielle.ciesielski@pnnl.gov
#
# Last update: 02 October 2020
#################################


## Import pandas read in excel raw data files and to store & manage dataframes.
## Import csv to read in and write out data for external processing.
## Import cobra to read in YLita649 model data and match genes to reactions.
import pandas, csv, cobra, os

## Identify experimental data file location (OMICS_FILE, OMICS_SHEET) and 
## which columns of the excel will be unnecessary; identify model data file 
## location (JSON_MODEL). 
path='/home/mollymersinger/bayesian-inference/'
os.chdir(path)
OMICS_FILE ='Omics_Data/23_strain_data/Proteomics/Yarrowia_protein_relative_mass.xlsx'
OMICS_SHEET = 'Mass percentage'
EXTRA_OMICS_COLUMNS = [2, 3, 4, 5, 6, 7]
JSON_MODEL = 'General_Information/YLita649/YLita649.itaconate.json'


#####
# Note that there are sometimes unusable rows in the raw data that can be
# identified by missing data in the 'Fasta headers' column. For simplicity, 
# rows that don't contain relevent data are considered unnecessary. Which 
# columns are unnecessary are determine ahead of time and identified in 
# EXTRA_OMICS_COLUMNS.
#####

## Import & clean raw data as Pandas dataframe `df'
df = pandas.read_excel(OMICS_FILE, sheet_name = OMICS_SHEET, header = 0)
df = df[df['Fasta headers'].map(type) == str] # remove unusable rows
df = df.drop(df.columns[EXTRA_OMICS_COLUMNS], axis = 1) # remove unnecessary columns


#####
# Note that the Protein IDs aren't always unique. To work around this, the ID
# used will always be the UniProt ID found in the Fasta header.
#####

## Extract the UniProt ID from the 'Fasta headers' and record them in a list.
uniprot_id_list = []
for uniprot_id in df['Fasta headers']:
    uniprot_id_list.append([uniprot_id[3:9]])

## Create a csv file to be read into the UniProt.org Retrieve/ID mapping tool
file = open('uniprot_id_list.csv', 'w+', newline = '')
with file:
    write = csv.writer(file)
    write.writerows(uniprot_id_list)


#####
# Go to https://www.uniprot.org/uploadlists/ and `upload your own file' in
# step 1 by locating the `uniprot_id_list.csv' file now located in the same
# location as this program. In step 2, choose `From: UniProtKB AC/ID' and 
# `To: Gene name', then submit. When the next page loads, click the download
# icon above the table and select `Mapping Table.' Copy and paste this data
# into an excel document. Change the information in the first row to `Maps.'
# Save this file as `mapping_table.csv' then proceed with the remaining code.
#####

## Read in the mapping table from the UniProt.org tool as a Pandas dataframe. 
## Note the UniProt IDs have been alphabetized and no longer map to the main 
## dataframe `df'.
map_df = pandas.read_csv('mapping_table.csv', sep = ',')
## Split the UniProt IDs and Gene names into separate columns in the dataframe.
map_df[['UniProt ID', 'Gene name']] = map_df.Maps.str.split('  ', expand = True)
## Remove the underscore from the Gene name so they match the Cobra model names.
map_df['Gene name'] = map_df['Gene name'].str.replace('_', '')
## Remove the raw UniProt.org data from the cleaned data.
map_df = map_df.drop(map_df.columns[[0]], axis = 1)

## Load in the reaction model.
model = cobra.io.load_json_model(JSON_MODEL)

## Create a list of all the gene names in the model as strings and use the 
## list to cross reference the gene names from the UniProt.org tool. Genes 
## that are in both sets are stored in a new dataframe `overlap_df.'
model_genes = [str(gene) for gene in list(model.genes)]
overlap_df = map_df.loc[map_df['Gene name'].isin(model_genes)]

## Create a list of the reactions associated with each gene in `overlap_df' 
## using their BiGG IDs.
rxn_lists = []
## Read each gene name through CobraPy to retrieve a Cobra object connecting
## the genes to their associated reactions and save them as a list objects.
for gene in overlap_df['Gene name']:
    cobra_gene = model.genes.get_by_id(gene)
    reactions = list(cobra_gene.reactions)
    ## Retrieve the BiGG ID for each reaction associated with the gene and 
    ## store them as a list.
    rxn_list_entry = []
    for rxn in reactions:
        ## Retrieve the reaction name from the Cobra object as a string.
        cobra_rxn = str(rxn)
        ## Retrieve the BiGG ID from the reaction name string.
        bigg_id = cobra_rxn.split(':')[0]
        ## Add the BiGG ID to the gene's list of reactions.
        rxn_list_entry.append(bigg_id)
    ## Add each gene's associated reactions to the list of associated reactions.   
    rxn_lists.append(rxn_list_entry)
## Append the BiGG IDs to the overlapping genes in `overlap_df'
overlap_df.insert(loc = 2, column = 'Reactions', value = rxn_lists)
overlap_df.to_csv('overlap.csv')
## Create three new rows in the main dataframe `df' to track the alternate IDs.
## `UniProt ID' is the unique Protein ID from the 'Fasta headers.'
df.insert(2, 'UniProt ID', value = pandas.Series([], dtype = str))
## 'Gene name' is the gene name in the style of the Cobra objects.
df.insert(3, 'Gene name', value = pandas.Series([], dtype = str))
## `Reactions' are a list of BiGG IDs retrieved from the Cobra model
df.insert(4, 'Reactions', value = pandas.Series([], dtype = str))


#####
# Because invalid rows were removed from the main dataframe `df,' the UniProt 
# IDs in the different dataframes must be matched to ensure the reaction IDs 
# are assigned appropriately in `df.'
#####

## Create series of the UniProt IDs in `map_df' and `overlap_df.' 
map_df_ids = map_df['UniProt ID']
overlap_df_ids = overlap_df['UniProt ID']
## For each row in the main dataframe, extract the UniProt ID from the Fasta 
## header as `match_str' and record the row index as an integer.
for header in df['Fasta headers']:
    match_str = header[3:9]
    df_idx = int(df[df['Fasta headers'] == header].index.values)
    ## Find the index of the UniProt ID matching `match_str' in `map_df.'
    map_idx = map_df_ids[map_df_ids == match_str].index[0]
    ## Assign the UniProt ID from `map_df' to the corresponding gene in `df.'
    df.loc[df_idx, 'UniProt ID'] = map_df.iloc[map_idx]['UniProt ID']
    ## Assign the gene name from `map_df' to the corresponding gene in `df.'
    df.loc[df_idx, 'Gene name'] = map_df.iloc[map_idx]['Gene name']
    ## If the gene is in `overlap_df,' find its index.
    if not overlap_df_ids[overlap_df_ids == match_str].empty:
        overlap_df_idx = overlap_df_ids[overlap_df_ids == match_str].index[0]
        ## Assign the associated reactions to the gene in `df' as a comma 
        ## separated string.
        df.loc[df_idx, 'Reactions'] = ', '.join(
            overlap_df.loc[overlap_df_idx,'Reactions']
            )        
df.to_csv('df.csv')
#####
# As a check, compare the first Protein ID in the first column to the UniProt
# ID in the finalized dataframe `df' to make sure they all match. 
#####

## Retrive the first six digits of each entry in the `Protein ID' column and 
## record them in a list.
raw_data_id = []
for id in df['Protein IDs']:
    raw_data_id.append(id[0:6])
## Retrieve the `UniProt ID ' series from `df' and save it as a list.
derived_data_id = list(df['UniProt ID'])
## Compare the entries of each list and add one to the counter for each error
## observed. If the count is zero at the end, there are no discrepancies.
counter = 0
for i in range(len(raw_data_id)):
    if raw_data_id[i] != derived_data_id[i]:
        counter += 1
