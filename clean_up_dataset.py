# use this script to clean up STCRDab export .tsv, to only one chain combination per each structure
# also extract sequences for each of your structures
# uses custom functions, which it looks for in a functions/ folder within repo

import pandas as pd
from Bio import PDB
from functions.myfunctions import *
import functions.sequencefunctions as sf
import functions.structurefunctions as stf
from datetime import date
import numpy as np
import pymol

structures_dir = "/path/to/stcrdab/download/imgt/"
structures_dir1 = "/path/where/I/save/PDB/from/pymol.fetch/"
structures_dir2 = "/path/to/stcrdab/download/raw/"

datalist = pd.read_csv("data/2022-05-31_tcrdb_export_cleaned.csv")
vdjdb_info = pd.read_csv("data/2022-06-27_vdjdb_structure_info.csv")

## a little manual dataset cleanup - I have curated these manually, because they behaved weird/are supercomplexes
idx = datalist.loc[(datalist.pdb == "5yxu") & (datalist.Achain == "F")].index.values.tolist()
idx1 = datalist.loc[(datalist.pdb == "6rp9") & (datalist.Achain == "K")].index.values.tolist()
idx2 = datalist.loc[(datalist.pdb == "4c56")].index.values.tolist()
idx3 = datalist.loc[(datalist.pdb == "4ei6")].index.values.tolist() # weird, one chain is missing one res (can't figure out why)
idx4 = datalist.loc[(datalist.pdb == "2p1y")].index.values.tolist() # singlechain, but chains are weird in datalist
idx5 = datalist.loc[(datalist.pdb == "2wbj") & (datalist.Achain == "G")].index.values.tolist()
idx6 = datalist.loc[(datalist.pdb == "4qrp") & (datalist.Achain == "K")].index.values.tolist()

idx0 = datalist.loc[datalist.pdb == "1kb5"].index.values.tolist()
datalist.at[idx0[0], "antigen_chain"] = np.nan
datalist_clean = datalist.drop(index= idx+idx1+idx2+idx3+idx4+idx5+idx6)

p = PDB.PDBParser()

pymol.cmd.set('fetch_path', structures_dir1, quiet=0)
pymol.cmd.set("cif_keepinmemory")
datalist_unique_pdbs = pd.DataFrame() # keep the most complete of each for avg
seqs = {}
# removed_pdbs = []
singlechain = ["2p1y"]
missing_residues_inAorB = []
no_vdjdb_val = []
rows_to_remove = []
pdbnotes = {}
i = 1
j = 1

for num, pdb in enumerate(list(set(datalist_clean.pdb))):

    # get information about pdb and handle special cases
    achain = list(datalist_clean.loc[datalist_clean.pdb == pdb]["Achain"])
    bchain = list(datalist_clean.loc[datalist_clean.pdb == pdb]["Bchain"])
    echain = list(datalist_clean.loc[datalist_clean.pdb == pdb]["antigen_chain"])

    if pdb in ["4z7w", "4lcc", "2gj6"]:
        # 4z7w has the ag encoded by two chains because it's a peptide + sugar (I think...)
        # 4lcc is mait TCR, it has C|C has epitope chain (but that's MR1, because no peptide)
        echain = [x.split("|")[0].strip() for x in echain]
    
    if [x.upper() for x in achain] == [x.upper() for x in bchain]:
        print("Single chain TCR, ", pdb, ", total count: ", j)
        j += 1
        singlechain.append(pdb)
        continue

    chain_pairs = list(zip(achain, bchain, echain))

    s = pdb + ".pdb"
    print(num, pdb)

    pymol.cmd.fetch(pdb)
    
    manual = stf.get_missing_residues_from_header(structures_dir2, pdb)
    if manual is not None:
        print("checking out missing residues")
        manual.SSSEQI = manual.SSSEQI.astype("str")
        # automatic = stf.bioPDB_missing_residues(structures_dir2, pdb)
        automatic = stf.pymol_missing_residues(pdb)
        automatic.SSSEQI = automatic.SSSEQI.astype("str")
        try:
            assert manual.equals(automatic)
        except:
            print(pdb, " - missing residues lists do not correspond")
            problematic_chains = set(manual.loc[~manual.SSSEQI.eq(automatic.SSSEQI)].C)
            print(manual)
            print(automatic)
            achain = [x for x in achain if x not in problematic_chains]
            bchain = [x for x in bchain if x not in problematic_chains]
            chain_pairs = [x for x in chain_pairs if x[0] in achain and x[1] in bchain]
            if len(achain) > 0 and len(bchain)>0 and len(chain_pairs)>0:
                print(pdb, " has missing residues, but not in all A or B")
                print("allowed pairs: ", chain_pairs)
                rows_to_remove += datalist_clean.loc[(datalist_clean.pdb == pdb) & ((datalist_clean.Achain.isin(problematic_chains)) | (datalist_clean.Bchain.isin(problematic_chains)))].index.values.tolist()
                pass
            else:
                print(pdb, " has problematic missing residues, skipped")
                missing_residues_inAorB.append(pdb)
                continue
    else:
        automatic = stf.pymol_missing_residues(pdb)
        assert manual == automatic, "ERROR: Pymol finds missing residues, " + pdb

    
    ## expected sequence (pymol includes the missing residues)
    # if you only extract sequence from the renumbered pdb, it might miss unresolved residues
    
    pymolseqs = pymol.cmd.get_fastastr(pdb).split(">")
    pymol.cmd.delete(pdb)
    pymolseqs = [x.replace("\n", "") for x in pymolseqs if x!=""]
    pymolseqs_d = {x[5:6]:x[8:] for x in pymolseqs}

    if pdb == "6vm8":
        pymolseqs_d["D"] = pymolseqs_d["D"].strip("?") #not present in the fasta on PDB

    # renumber pymol sequences to find sequence of v region only

    pymol_imgt_a = {a:sf.get_vregion_seq(pymolseqs_d[a]) for a in achain}
    pymol_imgt_a_gaps = {a:sf.get_vregion_seq_with_gaps(pymolseqs_d[a]) for a in achain}
    pymol_imgt_b = {b:sf.get_vregion_seq(pymolseqs_d[b]) for b in bchain}
    pymol_imgt_b_gaps = {b:sf.get_vregion_seq_with_gaps(pymolseqs_d[b]) for b in bchain}

    assert len(set(pymol_imgt_a.values())) == 1, pymol_imgt_a # sanity check that all pymol sequences are the same
    assert len(set(pymol_imgt_b.values())) == 1, pymol_imgt_a_gaps # sanity check that all pymol sequences are the same
    assert len(set(pymol_imgt_a_gaps.values())) == 1, pymol_imgt_b # sanity check that all pymol sequences are the same
    assert len(set(pymol_imgt_b_gaps.values())) == 1, pymol_imgt_b_gaps # sanity check that all pymol sequences are the same

    ## check cdr3s are what I expect them to be, since I can get that info from vdjdb

    if pdb in set(vdjdb_info["structure.id"]):
        cdr3a = vdjdb_info.loc[(vdjdb_info["structure.id"] == pdb) & (vdjdb_info.Gene == "TRA"), "CDR3"].tolist()
        cdr3b = vdjdb_info.loc[(vdjdb_info["structure.id"] == pdb) & (vdjdb_info.Gene == "TRB"), "CDR3"].tolist()

        # print(cdr3a)

        assert len(set(cdr3a)) == 1
        assert len(set(cdr3b)) == 1

        assert cdr3a[0] == sf.get_cdr3_seq(pymolseqs_d[achain[0]])
        assert cdr3b[0] == sf.get_cdr3_seq(pymolseqs_d[bchain[0]])
    else:
        print("cdr3 sequence cannot be validated as not in vdjdb")
        no_vdjdb_val.append(pdb)

    ## look at imgt pdb file 
    ## make a note of which chains have the correct amino acids in the v region
    struc = p.get_structure(pdb, structures_dir + s)
    sa = {a:sf.get_chain_seq(struc[0][a]) for a in achain}
    sb = {b:sf.get_chain_seq(struc[0][b]) for b in bchain}

    # handle epitope (less concerned about this)

    pymol_e = {}
    se = {}

    for e in echain:
        if pd.isna(e):
            pymol_e[e] = ""
            se[e] = ""
        else:
            se[e] = sf.get_chain_seq(struc[0][e])
            pymol_e[e] = pymolseqs_d[e]

    assert len(set(pymol_e.values())) == 1, pymol_e
   
    sa_ann = {a:(sa[a], sa[a] in pymol_imgt_a[a] or pymol_imgt_a[a] in sa[a]) for a in sa.keys()}
    sb_ann = {b:(sb[b], sb[b] in pymol_imgt_b[b] or pymol_imgt_b[b] in sb[b]) for b in sb.keys()}
    se_ann = {e:(se[e], se[e] == pymol_e[e]) for e in se.keys()}
    pdbnotes[pdb] = {}
    pdbnotes[pdb] = {(a,b,e):(sa_ann[a][1], sb_ann[b][1], se_ann[e][1]) for (a, b, e) in chain_pairs}

    # write down if in complex or not

    mhc = set(datalist_clean.loc[datalist_clean.pdb == pdb]["mhc_type"].isna())
    assert len(mhc) == 1

    complex = 1-int(list(mhc)[0])

    # save V region sequence to dictionary

    seqs[pdb] = {"alpha_aa_imgt":list(set(pymol_imgt_a.values()))[0], "alpha_aa_imgt_withGaps":list(set(pymol_imgt_a_gaps.values()))[0], 
                "beta_aa_imgt":list(set(pymol_imgt_b.values()))[0], "beta_aa_imgt_withGaps":list(set(pymol_imgt_b_gaps.values()))[0], "epitope_aa":list(set(pymol_e.values()))[0]}

datalist_clean = datalist_clean.drop(index = rows_to_remove)
sequences = pd.DataFrame.from_dict(seqs, orient = "index")
data_with_seq = pd.merge(datalist_clean, sequences, left_on="pdb", right_index=True)
data_with_seq_nodupl = data_with_seq.drop_duplicates(subset = ["alpha_aa_imgt", "beta_aa_imgt", "epitope_aa"])
best_combo = {}
best_pdb_set = pd.DataFrame()
for pdb in pdbnotes.keys():
    info = pdbnotes[pdb]
    print(info)
    for chain_mix in info:
        # print(chain_mix)
        if info[chain_mix][0]==True and info[chain_mix][1]==True and info[chain_mix][2] == True:
            best_combo[pdb] = {chain_mix:info[chain_mix]}
            break
    if pdb not in best_combo.keys():
        for chain_mix in info:
            if info[chain_mix][0]==True and info[chain_mix][1]==True:
                best_combo[pdb] = {chain_mix:info[chain_mix]}
                break
    if pdb not in best_combo.keys():
        for chain_mix in info:
            if info[chain_mix][0]==True or info[chain_mix][1]==True:
                best_combo[pdb] = {chain_mix:info[chain_mix]}
                break
    if pdb not in best_combo.keys():
        best_combo[pdb] = {chain_mix:info[list(info.keys())[0]]}

    print(best_combo[pdb])
    achain = list(best_combo[pdb].keys())[0][0]
    # print(achain)     
    df_line = data_with_seq.loc[(data_with_seq.pdb == pdb) & (data_with_seq.Achain == achain)].reset_index(drop = True)
    # print(df_line)
    assert df_line.shape[0] == 1
    best_pdb_set = pd.concat([best_pdb_set, df_line])
