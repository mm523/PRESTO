import numpy as np
import matplotlib.pyplot as pl
import mdtraj as md
import matplotlib as mpl
#from nglview.player import TrajectoryPlayer
import os
from Bio import pairwise2
import pandas
import itertools
from matplotlib import rc
from matplotlib import rcParams

# Woohoo my own homemade scripts!
import tcr_structParse as tcrParse

# Define MHC class 1 or 2
mhc_class = 1

# Load in a TCR3D metadata file to get out PDB names
if mhc_class == 1:
    input_df = pandas.read_csv('tcr3d_classImeta.csv')
elif mhc_class == 2:
    input_df = pandas.read_csv('tcr3d_classIImeta.csv')

# Only pick out human TCRs
human_df = input_df[input_df['Species'] == 'Human']
print(len(human_df))
print('Is when this script will end')

# Run through all of the PDBs in this database. Pretty slow overall
for dataset in np.arange(len(human_df)):
    # Pull out PDB name
    pdb = human_df['PDB<BR>ID'].values[dataset]
    # Pull out TRAV name for later reference
    trav = human_df['TRAV<BR>gene'].values[dataset]
    # Pull out TRBV name for later reference
    trbv = human_df['TRBV<BR>gene'].values[dataset]
    # Really only use this next one for classII
    mhcID_pre = human_df['MHC<BR>Name'].values[dataset]
    
    # Only need this for class II
    # DQ/DR/DP need different alignment sequences
    if pandas.isna(mhcID_pre):
        # Default to HLA-DQ and see if it works
        mhcID = 'HLA-DQ'
    else:
        # This won't capture the full class I allele
        # but we don't care really
        mhcID = mhcID_pre[0:6]
    
    # Unclear why but this PDB crashes things, check it out
    if pdb == '6R0E':
        continue
    # Oh hell yea we can load in PDBs from URL
    # It's pretty fast too, from what I can tell
    struct = md.load_pdb('http://www.rcsb.org/pdb/files/'+pdb+'.pdb')
    # Get a table of the atoms in the pdb
    table,bonds = struct.topology.to_dataframe()

    # Pull out the proper chains and chain names from the pdb
    xxx = tcrParse.get_chains(struct,mhc_class = mhc_class,mhcID = mhcID)
    
    # Double check that our function worked, and then separate out the chains
    if len(xxx) == 0:
        print('No chains?')
        continue
    else:
        if mhc_class == 1:
            alpha_chain,alpha_nameF,beta_chain,beta_nameF,mhc_chain = xxx
        elif mhc_class == 2:
            alpha_chain,alpha_nameF,beta_chain,beta_nameF,mhc_alpha_chain,mhc_beta_chain = xxx
            
    # Go back and double check if we have TRAV/TRBV mismatches
    if alpha_nameF != trav:
        print('ERROR: TRAV mismatch PDB: '+pdb)
        print('Should be: ' +trav)
        print('IDed as: '+alpha_nameF)
        continue
    if beta_nameF != trbv:
        print('ERROR: TRBV mismatch PDB: '+pdb)
        print('Should be: ' +trbv)
        print('IDed as: '+beta_nameF)
        continue

    # Really annoying there are so many inputs to this script, but leave it be for now
    # These scripts automatically pull out sidechain-sidechain interactions with some
    # custom restrictions. Check these scripts for details
    if mhc_class == 1:
        alpha_df = tcrParse.calc_process_dist(struct,alpha_chain,mhc_chain,alpha_nameF,beta_nameF,table,ab='alpha',dist_cutoff = 0.35)
        beta_df = tcrParse.calc_process_dist(struct,beta_chain,mhc_chain,alpha_nameF,beta_nameF,table,ab='beta',dist_cutoff = 0.35)
    elif mhc_class == 2:
        alpha_output = tcrParse.calc_process_classIIdist(struct,alpha_chain,mhc_alpha_chain,mhc_beta_chain,
                                                         alpha_nameF,beta_nameF,table,ab='alpha',dist_cutoff = 0.35,mhcID = mhcID)
        beta_output = tcrParse.calc_process_classIIdist(struct,beta_chain,mhc_alpha_chain,mhc_beta_chain,
                                                        alpha_nameF,beta_nameF,table,ab='beta',dist_cutoff = 0.35,mhcID = mhcID)
        
        # Cancel the loop if we get a BadReg error
        # BadReg means there are a bunch of weird breaks and mismatches in the PDB
        # Hopefully can rectify these issues later on
        if alpha_output == 'BadReg':
            continue
        if beta_output == 'BadReg':
            continue
        
        # Populate our distance dataframe from the output of the structural analysis
        if len(alpha_output) != 0:
            tcrA_alpha_df = alpha_output[0]; tcrA_beta_df = alpha_output[1]
            if len(alpha_output[0]) != 0:
                meta1_df = pandas.DataFrame(['mhc_alpha']*len(alpha_output[0]))
                meta1_df.columns = ['mhc_chain']
                alpha_pre1 = pandas.concat([alpha_output[0],meta1_df],axis=1)
                alpha_pre1.columns = [['tcrRes', 'tcrNum', 'tcrName', 'mhcRes', 'mhcNum', 'mhcName', 'distance', 'loop','mhc_chain']]
            else:
                alpha_pre1 = pandas.DataFrame([])

            if len(alpha_output[1]) != 0:
                meta2_df = pandas.DataFrame(['mhc_beta']*len(alpha_output[1]))
                meta2_df.columns = ['mhc_chain']
                alpha_pre2 = pandas.concat([alpha_output[1],meta2_df],axis=1)
                alpha_pre2.columns = [['tcrRes', 'tcrNum', 'tcrName', 'mhcRes', 'mhcNum', 'mhcName', 'distance', 'loop','mhc_chain']]
            else:
                alpha_pre2 = pandas.DataFrame([])
            
            alpha_df = pandas.concat([alpha_pre1,alpha_pre2])
        else:
            alpha_df = []
        
        # Same deal as above, but now do it for the beta chain-MHC interactions
        if len(beta_output) != 0:
            tcrB_alpha_df = beta_output[0]; tcrB_beta_df = beta_output[1]
            if len(beta_output[0]) != 0:
                meta1_df = pandas.DataFrame(['mhc_alpha']*len(beta_output[0]))
                meta1_df.columns = ['mhc_chain']
                beta_pre1 = pandas.concat([beta_output[0],meta1_df],axis=1)
                beta_pre1.columns = [['tcrRes', 'tcrNum', 'tcrName', 'mhcRes', 'mhcNum', 'mhcName', 'distance', 'loop','mhc_chain']]
            else:
                beta_pre1 = pandas.DataFrame([])

            if len(beta_output[1]) != 0:
                meta2_df = pandas.DataFrame(['mhc_beta']*len(beta_output[1]))
                meta2_df.columns = ['mhc_chain']
                beta_pre2 = pandas.concat([beta_output[1],meta2_df],axis=1)
                beta_pre2.columns = [['tcrRes', 'tcrNum', 'tcrName', 'mhcRes', 'mhcNum', 'mhcName', 'distance', 'loop','mhc_chain']]
            else:
                beta_pre2 = pandas.DataFrame([])
                
            beta_df = pandas.concat([beta_pre1,beta_pre2])
        else:
            beta_df = []
        
        if len(alpha_output) == 0 and len(beta_output) == 0:
            continue

    # Most of these that "don't have germline contacts"
    # are actually just poorly deposited PDBs        
    if len(alpha_df) == 0 and len(beta_df) == 0:
        print('No Germline Contacts PDB: '+pdb)
        continue
    elif len(alpha_df) == 0:
        final_df = beta_df
    elif type(alpha_df) == str:
        # I need to double check, but I think I've made this a redundant failsafe?
        # It might be already covered further up in the code?
        if alpha_df == 'BadReg':
            print('ERROR: TCR DistCalc Fail PDB: '+pdb)
            continue
        elif beta_df == 'BadSel':
            print('ERROR: TCR DistCalc Fail PDB: '+pdb)
            continue
    elif len(beta_df) == 0:
        final_df = alpha_df
    elif type(beta_df) == str:
        if beta_df == 'BadReg':
            print('ERROR: TCR DistCalc Fail PDB: '+pdb)
            continue
        elif beta_df == 'BadSel':
            print('ERROR: TCR DistCalc Fail PDB: '+pdb)
            continue
    else:
        final_df = pandas.concat([alpha_df,beta_df])
    
    # Save every one of these processed distances to a csv
    # This script is slow enough you'd rather not need to run
    # it multiple times
    final_df.to_csv(pdb+'_dists.csv',index=False)

########## End of script! ########