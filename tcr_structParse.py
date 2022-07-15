# Leave credit at the top here because at some point in development I made a dumb mistake and deleted 
# a non-backed up version of this python file... Oof
# uncompyle6 version 3.8.0
# Python bytecode 2.7 (62211)
# Decompiled from: Python 2.7.18 |Anaconda, Inc.| (default, Apr 23 2020, 17:44:47) 
# [GCC 4.2.1 Compatible Clang 4.0.1 (tags/RELEASE_401/final)]

# Embedded file name: tcr_structParse.py
# Compiled at: 2022-05-02 12:35:38
import numpy as np, matplotlib.pyplot as pl, matplotlib as mpl, mdtraj as md, os
from Bio import pairwise2
import pandas, itertools

def convert_3Let(inp):
    first = True
    three_let = ['ALA', 'GLY', 'ARG', 'LYS', 'ASP', 'GLU', 'ASN', 'GLN', 'MET', 'CYS', 'PHE', 'TYR', 'THR', 'TRP', 'PRO', 'SER', 'LEU', 'VAL', 'HIS', 'ILE']
    sin_let = ['A', 'G', 'R', 'K', 'D', 'E', 'N', 'Q', 'M', 'C', 'F', 'Y', 'T', 'W', 'P', 'S', 'L', 'V', 'H', 'I']
    sin_final = []
    for i in inp:
        hold = []
        for scan in np.arange(len(three_let)):
            if i.lower() == three_let[scan].lower():
                hold = sin_let[scan]
                break

        if len(hold) == 0:
            continue
        if first:
            sin_final = hold
            first = False
        else:
            sin_final = np.hstack((sin_final, hold))

    if len(sin_final) == 0:
        return ()
    return sin_final

# Have HLA-DP as a standard
def get_chains(struct, mhc_class = 1,mhcID='HLA-DP'):
    trav_genes = pandas.read_csv('trav_human_full.csv')
    trbv_genes = pandas.read_csv('trbv_human_full.csv')
    if mhc_class == 1:
        mhc_seq = 'GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQKMEPRAPWIEQEGPEYWDQETRNMKAHSQTDRANLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAVHAAEQRRVYLEGRCVDGLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWE'
    elif mhc_class == 2:
        if mhcID == 'HLA-DP':
            mhc_alpha = 'KADHVSTYAAFVQTHRPTGEFMFEFDEDEMFYVDLDKKETVWHLEEFGQAFSFEAQGGLANIAILNNNLNTLIQRSNHT'
        elif mhcID == 'HLA-DQ':
            mhc_alpha = 'ADHVASCGVNLYQFYGPSGQYTHEFDGDEEFYVDLERKETAWRWPEFSKFGGFDPQGALRNMAVAKHNLNIMIKRYNST'
        elif mhcID == 'HLA-DR':
            mhc_alpha = 'KEEHVIIQAEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKANLEIMTKRSNYT'
        mhc_beta = 'MVLLTSVVQGRATPENYVYQGRQECYAFNGTQRFLERYIYNREEFVRFDSDVGEFRAVTELGRPDEDYWNSQKDLLEEKRAVPDRVCRHNYELDEAVTLQ'

    table, bonds = struct.topology.to_dataframe()
    totChain = table['chainID'].drop_duplicates().values[(-1)]
    alpha_chain = -1
    beta_chain = -1
    if mhc_class == 1:
        mhc_chain = -1
    elif mhc_class == 2:
        mhc_alpha_chain = -1
        mhc_beta_chain = -1
    alphaMax = 0
    betaMax = 0
    for chainNum in np.arange(totChain + 1):
        alpha_score = 0
        beta_score = 0
        chain_seq = [ residue for residue in struct.topology.chain(chainNum).residues ]
        seq_only = [ str(i)[:3] for i in chain_seq ]
        sinLet_seq = convert_3Let(seq_only)
        final_sequence = ('').join(sinLet_seq)
        if mhc_class == 1:
            aa = pairwise2.align.globalms(mhc_seq, final_sequence, 0.5, -0.1, -5, -0.5)
        elif mhc_class == 2:
            # Need to change the scoring a bit for class II, I think because of shorter sequence
            # in each individual matching
            aa_alpha = pairwise2.align.globalms(mhc_alpha, final_sequence, 1.5, -0.1, -5, -0.5)
            aa_beta = pairwise2.align.globalms(mhc_beta, final_sequence, 1.5, -0.1, -5, -0.5)
            if aa_alpha == [] or aa_beta == []:
                aa = []
            else:
                if aa_alpha[0][2] > aa_beta[0][2]:
                    aa = aa_alpha
                    classII_chain = 'alpha'
                else:
                    aa = aa_beta
                    classII_chain = 'beta'
        if aa == []:
            score = 0
        else:
            score = aa[0][2]
        if score > 0:
            if mhc_class == 1:
                mhc_chain = chainNum
            else:
                if classII_chain == 'alpha':
                    mhc_alpha_chain = chainNum
                elif classII_chain == 'beta':
                    mhc_beta_chain = chainNum
            continue
        for tcrA in trav_genes.values:
            temp_name = tcrA[0]
            temp_seq = tcrA[1]
            aa = pairwise2.align.localms(temp_seq, final_sequence, 0.5, -0.1, -5, -0.5)
            if aa == []:
                score = 0
            else:
                score = aa[0][2]
            if score > alpha_score:
                alpha_score = score
            if score > alphaMax:
                alpha_nameF = temp_name
                alphaMax = score

        for tcrB in trbv_genes.values:
            temp_name = tcrB[0]
            temp_seq = tcrB[1]
            aa = pairwise2.align.localms(temp_seq, final_sequence, 0.5, -0.1, -5, -0.5)
            if aa == []:
                score = 0
            else:
                score = aa[0][2]
            if score > beta_score:
                beta_score = score
            if score > betaMax:
                beta_nameF = temp_name
                betaMax = score

        if alpha_score > 20 and alpha_score >= alphaMax:
            alpha_chain = chainNum
        if beta_score > 20 and beta_score >= betaMax:
            beta_chain = chainNum
        if mhc_class == 1:
            if alpha_chain != -1 and beta_chain != -1 and mhc_chain != -1:
                break
        elif mhc_class == 2:
            if alpha_chain != -1 and beta_chain != -1 and mhc_alpha_chain != -1 and mhc_beta_chain != -1:
                break

    if alpha_chain == -1:
        print('Cannot find TRAV match!')
        return ()
    if beta_chain == -1:
        print('Cannot find TRBV match!')
        return ()
    if mhc_class == 1:
        if mhc_chain == -1:
            print('Cannot find MHC match!')
            return ()
        return (alpha_chain, alpha_nameF, beta_chain, beta_nameF, mhc_chain)
    elif mhc_class == 2:
        if mhc_alpha_chain == -1:
            print('Cannot find MHCalpha match!')
            return ()
        if mhc_beta_chain == -1:
            print('Cannot find MHCbeta match!')
        return (alpha_chain, alpha_nameF, beta_chain, beta_nameF, mhc_alpha_chain,mhc_beta_chain)

def calc_process_dist(struct, tcr_chain, mhc_chain, alpha_nameF, beta_nameF, table, ab='alpha', dist_cutoff=0.35):

    # Try to edit this to make sure our numbers match from the get-go
    alpha1 = [55,56,59,60,63,66,67,70,74,77,80]
    pep_contact = [3,5,7,20,22,24,43,57,61,64,65,68,71,72,75,78,79,82,93,95,97,112,114,141,145,150,154,157,169]
    alpha2 = [143,144,147,148,149,152,153,156,160,161,164,165,167,168]

    hlaA_0101 = 'SHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQKMEPRAPWIEQEGPEYWDQETRNMKAHSQTDRANLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAVHAAEQRRVYLEGRCVDGLRRYLENGKETL'

    trav_cdrs = pandas.read_csv('trav_human_cdrs.csv')
    trbv_cdrs = pandas.read_csv('trbv_human_cdrs.csv')
    trav_12seq = trav_cdrs[(trav_cdrs['gene'] == alpha_nameF)][['cdr1', 'cdr2']].values[0]
    trbv_12seq = trbv_cdrs[(trbv_cdrs['gene'] == beta_nameF)][['cdr1', 'cdr2']].values[0]
    tcr_sub = [ residue for residue in struct.topology.chain(tcr_chain).residues ]
    seq_only1 = [ str(i)[:3] for i in tcr_sub ]
    num_only1 = [ str(i)[3:] for i in tcr_sub ]
    sinLet_seq1 = convert_3Let(seq_only1)
    tcr_sequence = ('').join(sinLet_seq1)
    if ab == 'alpha':
        tcr_cdr1Start = tcr_sequence.find(trav_12seq[0])
        tcr_cdr2Start = tcr_sequence.find(trav_12seq[1])
        if tcr_cdr1Start == -1:
            tcr_cdr1Start = tcr_sequence.find(trav_12seq[0][0:3])
        if tcr_cdr1Start == -1:
            aa = pairwise2.align.localms(trav_12seq[0], tcr_sequence, 0.5, -0.1, -5, -0.5)
            tcr_cdr1Start = aa[0][0].find(trav_12seq[0])
        if tcr_cdr2Start == -1:
            tcr_cdr2Start = tcr_sequence.find(trav_12seq[1][0:3])
        if tcr_cdr2Start == -1:
            aa = pairwise2.align.localms(trav_12seq[1], tcr_sequence, 0.5, -0.1, -5, -0.5)
            tcr_cdr2Start = aa[0][0].find(trav_12seq[1])
        tcr_cdr1End = tcr_cdr1Start + len(trav_12seq[0])
        tcr_cdr2End = tcr_cdr2Start + len(trav_12seq[1])
    else:
        if ab == 'beta':
            tcr_cdr1Start = tcr_sequence.find(trbv_12seq[0])
            tcr_cdr2Start = tcr_sequence.find(trbv_12seq[1])
            if tcr_cdr1Start == -1:
                tcr_cdr1Start = tcr_sequence.find(trbv_12seq[0][0:3])
            if tcr_cdr2Start == -1:
                tcr_cdr2Start = tcr_sequence.find(trbv_12seq[1][0:3])
            tcr_cdr1End = tcr_cdr1Start + len(trbv_12seq[0])
            tcr_cdr2End = tcr_cdr2Start + len(trbv_12seq[1])
    mhc_sub = [ residue for residue in struct.topology.chain(mhc_chain).residues ]
    seq_only3 = [ str(i)[:3] for i in mhc_sub ]
    num_only3 = [ str(i)[3:] for i in mhc_sub ]
    sinLet_seq3 = convert_3Let(seq_only3)
    mhc_sequence = ('').join(sinLet_seq3)

    aa = pairwise2.align.localms(hlaA_0101,mhc_sequence,1.0, -0.1, -5,-0.5)
    num_shift = aa[0][0].find('SHSMRYF')
    resid_shift = 1 - int(num_only3[0])
    
    cdr1_len = tcr_cdr1End - tcr_cdr1Start
    cdr2_len = tcr_cdr2End - tcr_cdr2Start
    tcrdiff = 0
    tcrdiffCount = 0
    for i in np.arange(cdr1_len):
        mdsel = struct.topology.select('chainid == ' + str(tcr_chain) + ' and (residue ' + num_only1[(tcr_cdr1Start + i)] + ') and name CA')
        CAs = table[(table['name'] == 'CA')]
        chainCA = CAs[(CAs['chainID'] == tcr_chain)]
        tabsel = chainCA[(chainCA['resSeq'] == int(num_only1[(tcr_cdr1Start + i)]))]['serial'].values
        if len(mdsel) > 1:
            mdsel = mdsel[0]
        if len(tabsel) > 1:
            tabsel = tabsel[0]
        newdiff = mdsel - tabsel
        if newdiff != tcrdiff:
            tcrdiffCount += 1
            tcrdiff = newdiff

    for i in np.arange(cdr2_len):
        mdsel = struct.topology.select('chainid == ' + str(tcr_chain) + ' and (residue ' + num_only1[(tcr_cdr2Start + i)] + ') and name CA')
        CAs = table[(table['name'] == 'CA')]
        chainCA = CAs[(CAs['chainID'] == tcr_chain)]
        tabsel = chainCA[(chainCA['resSeq'] == int(num_only1[(tcr_cdr2Start + i)]))]['serial'].values
        if len(mdsel) > 1:
            mdsel = mdsel[0]
        if len(tabsel) > 1:
            tabsel = tabsel[0]
        newdiff = mdsel - tabsel
        if newdiff != tcrdiff:
            tcrdiffCount += 1
            tcrdiff = newdiff

    mhcdiff = 0
    mhcdiffCount = 0
    for i in np.arange(len(mhc_sequence)):
        mdsel = struct.topology.select('chainid == ' + str(mhc_chain) + ' and (residue ' + num_only3[i] + ') and name CA')
        CAs = table[(table['name'] == 'CA')]
        chainCA = CAs[(CAs['chainID'] == mhc_chain)]
        tabsel = chainCA[(chainCA['resSeq'] == int(num_only3[i]))]['serial'].values
        if len(mdsel) > 1:
            mdsel = mdsel[0]
        if len(tabsel) > 1:
            tabsel = tabsel[0]
        newdiff = mdsel - tabsel
        if newdiff != mhcdiff:
            mhcdiffCount += 1
            mhcdiff = newdiff
    if mhcdiffCount != 1:
        print('Mid-sel register shift in MHC')
        return('BadReg')

    if tcrdiffCount != 1:
        print('Mid-sel register shift in TCR')
        return('BadReg')

    cdr1 = struct.topology.select('chainid == ' + str(tcr_chain) + ' and (residue ' + num_only1[tcr_cdr1Start] + ' to ' + num_only1[tcr_cdr1End] + ') and sidechain')
    cdr2 = struct.topology.select('chainid == ' + str(tcr_chain) + ' and (residue ' + num_only1[tcr_cdr2Start] + ' to ' + num_only1[tcr_cdr2End] + ') and sidechain')
    
    mhc_struct = struct.topology.select('chainid == ' + str(mhc_chain) + ' and sidechain')
    cdr1_mhc = list(itertools.product(cdr1, mhc_struct))
    cdr2_mhc = list(itertools.product(cdr2, mhc_struct))

    if len(cdr1) > 0:
        cdr1_dists = md.compute_distances(struct, atom_pairs=cdr1_mhc, periodic=False)
    else:
        return 'BadSel'
    if len(cdr2) > 0:
        cdr2_dists = md.compute_distances(struct, atom_pairs=cdr2_mhc, periodic=False)
    else:
        return 'BadSel'
    min_dist = dist_cutoff
    cdr1_pos = []
    cdr1_seldist = []
    for i in np.arange(len(cdr1_dists[0])):
        if cdr1_dists[0][i] < min_dist:
            cdr1_pos = cdr1_pos + [i]
            cdr1_seldist = cdr1_seldist + [cdr1_dists[0][i]]

    cdr2_pos = []
    cdr2_seldist = []
    for i in np.arange(len(cdr2_dists[0])):
        if cdr2_dists[0][i] < min_dist:
            cdr2_pos = cdr2_pos + [i]
            cdr2_seldist = cdr2_seldist + [cdr2_dists[0][i]]

    first = True
    for k in np.arange(len(cdr1_pos)):
        pairs = cdr1_pos[k]
        tcr_index = cdr1_mhc[pairs][0] - tcrdiff[0]
        mhc_index = cdr1_mhc[pairs][1] - mhcdiff[0]
        tcr_ID = table[(table['serial'] == tcr_index)][['resName', 'resSeq', 'name']].values[0]
        mhc_ID = table[(table['serial'] == mhc_index)][['resName', 'resSeq', 'name']].values[0]

        # Try to edit the 'resSeq'. Make it match the alpha1/2 numbering we've had
        ###################################
        for i in np.arange(len(num_only3)):
            if int(num_only3[i]) == mhc_ID[1]:
                start_seq = int(i)
        hold_seq = mhc_sequence[start_seq:start_seq+5]
        if hold_seq[0] != convert_3Let([mhc_ID[0]]):
            print('whoops')
        seq_loc = aa[0][1].find(hold_seq)

        for num in np.arange(len(alpha1)):
            if alpha1[num] == seq_loc-num_shift:
                mhc_ID[1] = alpha1[num]

        for num in np.arange(len(alpha2)):
            if alpha2[num] == seq_loc-num_shift:
                mhc_ID[1] = alpha2[num]
        ##################################

        if first:
            pre_df = np.hstack((tcr_ID, mhc_ID, cdr1_seldist[k], 'cdr1' + ab[0]))
            first = False
        else:
            pre_df = np.vstack((pre_df, np.hstack((tcr_ID, mhc_ID, cdr1_seldist[k], 'cdr1' + ab[0]))))

    for k in np.arange(len(cdr2_pos)):
        pairs = cdr2_pos[k]
        tcr_index = cdr2_mhc[pairs][0] - tcrdiff[0]
        mhc_index = cdr2_mhc[pairs][1] - mhcdiff[0]
        tcr_ID = table[(table['serial'] == tcr_index)][['resName', 'resSeq', 'name']].values[0]
        mhc_ID = table[(table['serial'] == mhc_index)][['resName', 'resSeq', 'name']].values[0]

        # Try to edit the 'resSeq'. Make it match the alpha1/2 numbering we've had
        ###################################
        for i in np.arange(len(num_only3)):
            if int(num_only3[i]) == mhc_ID[1]:
                start_seq = int(i)
        hold_seq = mhc_sequence[start_seq:start_seq+5]
        if hold_seq[0] != convert_3Let([mhc_ID[0]]):
            print('whoops')
        seq_loc = aa[0][1].find(hold_seq)

        for num in np.arange(len(alpha1)):
            if alpha1[num] == seq_loc-num_shift:
                mhc_ID[1] = alpha1[num]

        for num in np.arange(len(alpha2)):
            if alpha2[num] == seq_loc-num_shift:
                mhc_ID[1] = alpha2[num]
        ##################################

        if first:
            pre_df = np.hstack((tcr_ID, mhc_ID, cdr2_seldist[k], 'cdr2' + ab[0]))
            first = False
        else:
            pre_df = np.vstack((pre_df, np.hstack((tcr_ID, mhc_ID, cdr2_seldist[k], 'cdr2' + ab[0]))))

    if first == True:
        return ()
    final_df = pandas.DataFrame(pre_df)
    if np.shape(final_df)[1] != 8:
        final_df = np.transpose(final_df)
    final_df.columns = [['tcrRes', 'tcrNum', 'tcrName', 'mhcRes', 'mhcNum', 'mhcName', 'distance', 'loop']]
    return final_df
# okay decompiling tcr_structParse.pyc

def calc_process_classIIdist(struct, tcr_chain, mhc_alpha_chain, mhc_beta_chain, alpha_nameF, beta_nameF,
table, ab='alpha', dist_cutoff=0.35,mhcID = 'HLA-DP'):
    trav_cdrs = pandas.read_csv('trav_human_cdrs.csv')
    trbv_cdrs = pandas.read_csv('trbv_human_cdrs.csv')
    trav_12seq = trav_cdrs[(trav_cdrs['gene'] == alpha_nameF)][['cdr1', 'cdr2']].values[0]
    trbv_12seq = trbv_cdrs[(trbv_cdrs['gene'] == beta_nameF)][['cdr1', 'cdr2']].values[0]
    tcr_sub = [ residue for residue in struct.topology.chain(tcr_chain).residues ]
    seq_only1 = [ str(i)[:3] for i in tcr_sub ]
    num_only1 = [ str(i)[3:] for i in tcr_sub ]
    sinLet_seq1 = convert_3Let(seq_only1)
    tcr_sequence = ('').join(sinLet_seq1)
    if ab == 'alpha':
        tcr_cdr1Start = tcr_sequence.find(trav_12seq[0])
        tcr_cdr2Start = tcr_sequence.find(trav_12seq[1])
        if tcr_cdr1Start == -1:
            tcr_cdr1Start = tcr_sequence.find(trav_12seq[0][0:3])
        if tcr_cdr1Start == -1:
            aa = pairwise2.align.localms(trav_12seq[0], tcr_sequence, 0.5, -0.1, -5, -0.5)
            tcr_cdr1Start = aa[0][0].find(trav_12seq[0])
        if tcr_cdr2Start == -1:
            tcr_cdr2Start = tcr_sequence.find(trav_12seq[1][0:3])
        if tcr_cdr2Start == -1:
            aa = pairwise2.align.localms(trav_12seq[1], tcr_sequence, 0.5, -0.1, -5, -0.5)
            tcr_cdr2Start = aa[0][0].find(trav_12seq[1])
        tcr_cdr1End = tcr_cdr1Start + len(trav_12seq[0])
        tcr_cdr2End = tcr_cdr2Start + len(trav_12seq[1])
    else:
        if ab == 'beta':
            tcr_cdr1Start = tcr_sequence.find(trbv_12seq[0])
            tcr_cdr2Start = tcr_sequence.find(trbv_12seq[1])
            if tcr_cdr1Start == -1:
                tcr_cdr1Start = tcr_sequence.find(trbv_12seq[0][0:3])
            if tcr_cdr1Start == -1:
                aa = pairwise2.align.localms(trbv_12seq[0], tcr_sequence, 0.5, -0.1, -5, -0.5)
                tcr_cdr1Start = aa[0][0].find(trbv_12seq[0])
            if tcr_cdr2Start == -1:
                tcr_cdr2Start = tcr_sequence.find(trbv_12seq[1][0:3])
            if tcr_cdr2Start == -1:
                aa = pairwise2.align.localms(trbv_12seq[1], tcr_sequence, 0.5, -0.1, -5, -0.5)
                tcr_cdr2Start = aa[0][0].find(trbv_12seq[1])
            tcr_cdr1End = tcr_cdr1Start + len(trbv_12seq[0])
            tcr_cdr2End = tcr_cdr2Start + len(trbv_12seq[1])

    mhc_alpha_sub = [ residue for residue in struct.topology.chain(mhc_alpha_chain).residues ]
    mhc_beta_sub = [ residue for residue in struct.topology.chain(mhc_beta_chain).residues ]
    mhcalpha_seqonly3 = [ str(i)[:3] for i in mhc_alpha_sub ];mhcbeta_seqonly3 = [ str(i)[:3] for i in mhc_beta_sub ]
    mhcalpha_numonly3 = [ str(i)[3:] for i in mhc_alpha_sub ];mhcbeta_numonly3 = [ str(i)[3:] for i in mhc_beta_sub ]
    alphasinLet_seq3 = convert_3Let(mhcalpha_seqonly3); betasinLet_seq3 = convert_3Let(mhcbeta_seqonly3)
    mhcalpha_sequence = ('').join(alphasinLet_seq3); mhcbeta_sequence = ('').join(betasinLet_seq3)

    if mhcID == 'HLA-DP':
        mhc_alpha = 'KADHVSTYAAFVQTHRPTGEFMFEFDEDEMFYVDLDKKETVWHLEEFGQAFSFEAQGGLANIAILNNNLNTLIQRSNHT'
    elif mhcID == 'HLA-DQ':
        mhc_alpha = 'ADHVASCGVNLYQFYGPSGQYTHEFDGDEEFYVDLERKETAWRWPEFSKFGGFDPQGALRNMAVAKHNLNIMIKRYNST'
    elif mhcID == 'HLA-DR':
        mhc_alpha = 'KEEHVIIQAEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKANLEIMTKRSNYT'
    mhc_beta = 'MVLLTSVVQGRATPENYVYQGRQECYAFNGTQRFLERYIYNREEFVRFDSDVGEFRAVTELGRPDEDYWNSQKDLLEEKRAVPDRVCRHNYELDEAVTLQ'

    alpha = [51, 53, 55, 56, 58, 59, 62, 63, 66, 69, 70, 73, 74]
    beta = [63, 66, 67, 70, 71, 73, 74, 76, 77, 80, 83, 84, 87, 88, 91, 92]

    aa = pairwise2.align.localms(mhc_alpha,mhcalpha_sequence,1.0, -0.1, -5,-0.5)
    num_shift_alpha = aa[0][0].find(mhc_alpha[0:5])
    resid_shift_alpha = 1 - int(mhcalpha_numonly3[0])

    bb = pairwise2.align.localms(mhc_beta,mhcbeta_sequence,1.0, -0.1, -5,-0.5)
    num_shift_beta = bb[0][0].find(mhc_beta[0:5])
    resid_shift_beta = 1 - int(mhcbeta_numonly3[0])
    
    cdr1_len = tcr_cdr1End - tcr_cdr1Start
    cdr2_len = tcr_cdr2End - tcr_cdr2Start
    tcrdiff = 0
    tcrdiffCount = 0
    for i in np.arange(cdr1_len):
        mdsel = struct.topology.select('chainid == ' + str(tcr_chain) + ' and (residue ' + num_only1[(tcr_cdr1Start + i)] + ') and name CA')
        CAs = table[(table['name'] == 'CA')]
        chainCA = CAs[(CAs['chainID'] == tcr_chain)]
        tabsel = chainCA[(chainCA['resSeq'] == int(num_only1[(tcr_cdr1Start + i)]))]['serial'].values
        if len(mdsel) > 1:
            mdsel = mdsel[0]
        if len(tabsel) > 1:
            tabsel = tabsel[0]
        newdiff = mdsel - tabsel
        if newdiff != tcrdiff:
            tcrdiffCount += 1
            tcrdiff = newdiff

    for i in np.arange(cdr2_len):
        mdsel = struct.topology.select('chainid == ' + str(tcr_chain) + ' and (residue ' + num_only1[(tcr_cdr2Start + i)] + ') and name CA')
        CAs = table[(table['name'] == 'CA')]
        chainCA = CAs[(CAs['chainID'] == tcr_chain)]
        tabsel = chainCA[(chainCA['resSeq'] == int(num_only1[(tcr_cdr2Start + i)]))]['serial'].values
        if len(mdsel) > 1:
            mdsel = mdsel[0]
        if len(tabsel) > 1:
            tabsel = tabsel[0]
        newdiff = mdsel - tabsel
        if newdiff != tcrdiff:
            tcrdiffCount += 1
            tcrdiff = newdiff


    alphadiff = 0; betadiff = 0
    alphadiffCount = 0; betadiffCount = 0
    for i in np.arange(len(mhcalpha_sequence)):
        # Why TF are any of these pdbs numbered with a -1 resid? unclear
        if mhcalpha_numonly3[i].find('-') != -1:
            continue
        mdsel = struct.topology.select('chainid == ' + str(mhc_alpha_chain) + ' and (residue ' + mhcalpha_numonly3[i] + ') and name CA')
        CAs = table[(table['name'] == 'CA')]
        chainCA = CAs[(CAs['chainID'] == mhc_alpha_chain)]
        tabsel = chainCA[(chainCA['resSeq'] == int(mhcalpha_numonly3[i]))]['serial'].values
        if len(mdsel) > 1:
            mdsel = mdsel[0]
        if len(tabsel) > 1:
            tabsel = tabsel[0]
        newdiff = mdsel - tabsel
        if newdiff != alphadiff:
            alphadiffCount += 1
            alphadiff = newdiff
    if alphadiffCount != 1:
        print('Mid-sel register shift in MHCalpha')
        return('BadReg')

    for i in np.arange(len(mhcbeta_sequence)):
        # Why TF are any of these pdbs numbered with a -1 resid? unclear
        if mhcbeta_numonly3[i].find('-') != -1:
            continue
        mdsel = struct.topology.select('chainid == ' + str(mhc_beta_chain) + ' and (residue ' + mhcbeta_numonly3[i] + ') and name CA')
        CAs = table[(table['name'] == 'CA')]
        chainCA = CAs[(CAs['chainID'] == mhc_beta_chain)]
        tabsel = chainCA[(chainCA['resSeq'] == int(mhcbeta_numonly3[i]))]['serial'].values
        if len(mdsel) > 1:
            mdsel = mdsel[0]
        if len(tabsel) > 1:
            tabsel = tabsel[0]
        newdiff = mdsel - tabsel
        if newdiff != betadiff:
            betadiffCount += 1
            betadiff = newdiff
    if betadiffCount != 1:
        print('Mid-sel register shift in MHCbeta')
        return('BadReg')

    if tcrdiffCount != 1:
        print('Mid-sel register shift in TCR')
        return('BadReg')

    cdr1 = struct.topology.select('chainid == ' + str(tcr_chain) + ' and (residue ' + num_only1[tcr_cdr1Start] + ' to ' + num_only1[tcr_cdr1End] + ') and sidechain')
    cdr2 = struct.topology.select('chainid == ' + str(tcr_chain) + ' and (residue ' + num_only1[tcr_cdr2Start] + ' to ' + num_only1[tcr_cdr2End] + ') and sidechain')
    

    mhc_alpha_struct = struct.topology.select('chainid == ' + str(mhc_alpha_chain) + ' and sidechain')
    mhc_beta_struct = struct.topology.select('chainid == ' + str(mhc_beta_chain) + ' and sidechain')
    cdr1_mhcalpha = list(itertools.product(cdr1, mhc_alpha_struct))
    cdr1_mhcbeta = list(itertools.product(cdr1, mhc_beta_struct))
    cdr2_mhcalpha = list(itertools.product(cdr2, mhc_alpha_struct))
    cdr2_mhcbeta = list(itertools.product(cdr2, mhc_beta_struct))
    
    if len(cdr1) > 0:
        cdr1_dists_alpha = md.compute_distances(struct, atom_pairs=cdr1_mhcalpha, periodic=False)
        cdr1_dists_beta = md.compute_distances(struct, atom_pairs=cdr1_mhcbeta, periodic=False)
    else:
        return 'BadSel'
    if len(cdr2) > 0:
        cdr2_dists_alpha = md.compute_distances(struct, atom_pairs=cdr2_mhcalpha, periodic=False)
        cdr2_dists_beta = md.compute_distances(struct, atom_pairs=cdr2_mhcbeta, periodic=False)
    else:
        return 'BadSel'

    min_dist = dist_cutoff
    cdr1_pos_alpha = []; cdr1_seldist_alpha = []
    cdr1_pos_beta = []; cdr1_seldist_beta = []

    for i in np.arange(len(cdr1_dists_alpha[0])):
        if cdr1_dists_alpha[0][i] < min_dist:
            cdr1_pos_alpha = cdr1_pos_alpha + [i]
            cdr1_seldist_alpha = cdr1_seldist_alpha + [cdr1_dists_alpha[0][i]]

    for i in np.arange(len(cdr1_dists_beta[0])):
        if cdr1_dists_beta[0][i] < min_dist:
            cdr1_pos_beta = cdr1_pos_beta + [i]
            cdr1_seldist_beta = cdr1_seldist_beta + [cdr1_dists_beta[0][i]]

    cdr2_pos_alpha = []; cdr2_seldist_alpha = []
    cdr2_pos_beta = []; cdr2_seldist_beta = []

    for i in np.arange(len(cdr2_dists_alpha[0])):
        if cdr2_dists_alpha[0][i] < min_dist:
            cdr2_pos_alpha = cdr2_pos_alpha + [i]
            cdr2_seldist_alpha = cdr2_seldist_alpha + [cdr2_dists_alpha[0][i]]

    for i in np.arange(len(cdr2_dists_beta[0])):
        if cdr2_dists_beta[0][i] < min_dist:
            cdr2_pos_beta = cdr2_pos_beta + [i]
            cdr2_seldist_beta = cdr2_seldist_beta + [cdr2_dists_beta[0][i]]

    first_alpha = True
    for k in np.arange(len(cdr1_pos_alpha)):
        pairs = cdr1_pos_alpha[k]
        tcr_index = cdr1_mhcalpha[pairs][0] - tcrdiff[0]
        mhc_index = cdr1_mhcalpha[pairs][1] - alphadiff[0]
        tcr_ID = table[(table['serial'] == tcr_index)][['resName', 'resSeq', 'name']].values[0]
        mhc_ID = table[(table['serial'] == mhc_index)][['resName', 'resSeq', 'name']].values[0]

          # Try to edit the 'resSeq'. Make it match the alpha1/2 numbering we've had
        ###################################
        for i in np.arange(len(mhcalpha_numonly3)):
            if int(mhcalpha_numonly3[i]) == mhc_ID[1]:
                start_seq = int(i)
        hold_seq = mhcalpha_sequence[start_seq:start_seq+5]
        findit = aa[0][1].find(hold_seq)
        find_alpha=aa[0][0][findit:findit+5]
        for num in alpha:
            stand_alpha=mhc_alpha[num:num+5]
            if find_alpha == stand_alpha:
                mhc_ID[1] = num
                break
            else:
                mhc_ID[1] = 'DROP'
        ##################################
        if mhc_ID[1] == 'DROP':
            continue

        if first_alpha:
            pre_df_alpha = np.hstack((tcr_ID, mhc_ID, cdr1_seldist_alpha[k], 'cdr1' + ab[0]))
            first_alpha = False
        else:
            pre_df_alpha = np.vstack((pre_df_alpha, np.hstack((tcr_ID, mhc_ID, cdr1_seldist_alpha[k], 'cdr1' + ab[0]))))
    
    first_beta = True
    for k in np.arange(len(cdr1_pos_beta)):
        pairs = cdr1_pos_beta[k]
        tcr_index = cdr1_mhcbeta[pairs][0] - tcrdiff[0]
        mhc_index = cdr1_mhcbeta[pairs][1] - betadiff[0]
        tcr_ID = table[(table['serial'] == tcr_index)][['resName', 'resSeq', 'name']].values[0]
        mhc_ID = table[(table['serial'] == mhc_index)][['resName', 'resSeq', 'name']].values[0]

         # Try to edit the 'resSeq'. Make it match the alpha1/2 numbering we've had
        ###################################
        for i in np.arange(len(mhcbeta_numonly3)):
            if int(mhcbeta_numonly3[i]) == mhc_ID[1]:
                start_seq = int(i)
        hold_seq = mhcbeta_sequence[start_seq:start_seq+5]
        findit = bb[0][1].find(hold_seq)
        find_beta=bb[0][0][findit:findit+5]
        for num in beta:
            stand_beta=mhc_beta[num:num+5]
            if find_beta == stand_beta:
                mhc_ID[1] = num
                break
            else:
                mhc_ID[1] = 'DROP'
        if mhc_ID[1] == 'DROP':
            continue
        ##################################

        if first_beta:
            pre_df_beta = np.hstack((tcr_ID, mhc_ID, cdr1_seldist_beta[k], 'cdr1' + ab[0]))
            first_beta = False
        else:
            pre_df_beta = np.vstack((pre_df_beta, np.hstack((tcr_ID, mhc_ID, cdr1_seldist_beta[k], 'cdr1' + ab[0]))))

    for k in np.arange(len(cdr2_pos_alpha)):
        pairs = cdr2_pos_alpha[k]
        tcr_index = cdr2_mhcalpha[pairs][0] - tcrdiff[0]
        mhc_index = cdr2_mhcalpha[pairs][1] - alphadiff[0]
        tcr_ID = table[(table['serial'] == tcr_index)][['resName', 'resSeq', 'name']].values[0]
        mhc_ID = table[(table['serial'] == mhc_index)][['resName', 'resSeq', 'name']].values[0]

         # Try to edit the 'resSeq'. Make it match the alpha1/2 numbering we've had
        ###################################
        for i in np.arange(len(mhcalpha_numonly3)):
            if int(mhcalpha_numonly3[i]) == mhc_ID[1]:
                start_seq = int(i)
        hold_seq = mhcalpha_sequence[start_seq:start_seq+5]
        findit = aa[0][1].find(hold_seq)
        find_alpha=aa[0][0][findit:findit+5]
        for num in alpha:
            stand_alpha=mhc_alpha[num:num+5]
            if find_alpha == stand_alpha:
                mhc_ID[1] = num
                break
            else:
                mhc_ID[1] = 'DROP'
        ##################################
        if mhc_ID[1] == 'DROP':
            continue

        if first_alpha:
            pre_df_alpha = np.hstack((tcr_ID, mhc_ID, cdr2_seldist_alpha[k], 'cdr2' + ab[0]))
            first_alpha = False
        else:
            pre_df_alpha = np.vstack((pre_df_alpha, np.hstack((tcr_ID, mhc_ID, cdr2_seldist_alpha[k], 'cdr2' + ab[0]))))

    for k in np.arange(len(cdr2_pos_beta)):
        pairs = cdr2_pos_beta[k]
        tcr_index = cdr2_mhcbeta[pairs][0] - tcrdiff[0]
        mhc_index = cdr2_mhcbeta[pairs][1] - betadiff[0]
        tcr_ID = table[(table['serial'] == tcr_index)][['resName', 'resSeq', 'name']].values[0]
        mhc_ID = table[(table['serial'] == mhc_index)][['resName', 'resSeq', 'name']].values[0]

         # Try to edit the 'resSeq'. Make it match the alpha1/2 numbering we've had
        ###################################
        for i in np.arange(len(mhcbeta_numonly3)):
            if int(mhcbeta_numonly3[i]) == mhc_ID[1]:
                start_seq = int(i)
        hold_seq = mhcbeta_sequence[start_seq:start_seq+5]
        findit = bb[0][1].find(hold_seq)
        find_beta=bb[0][0][findit:findit+5]
        for num in beta:
            stand_beta=mhc_beta[num:num+5]
            if find_beta == stand_beta:
                mhc_ID[1] = num
                break
            else:
                mhc_ID[1] = 'DROP'
        if mhc_ID[1] == 'DROP':
            continue
        ##################################

        if first_beta:
            pre_df_beta = np.hstack((tcr_ID, mhc_ID, cdr2_seldist_beta[k], 'cdr2' + ab[0]))
            first_beta = False
        else:
            pre_df_beta = np.vstack((pre_df_beta, np.hstack((tcr_ID, mhc_ID, cdr2_seldist_beta[k], 'cdr2' + ab[0]))))

    if first_alpha == True and first_beta == True:
        print('no found contacts')
        return ()

    if first_alpha:
        final_df_alpha = []
    else:
        final_df_alpha = pandas.DataFrame(pre_df_alpha)
        if np.shape(final_df_alpha)[1] != 8:
            final_df_alpha = np.transpose(final_df_alpha)
        final_df_alpha.columns = [['tcrRes', 'tcrNum', 'tcrName', 'mhcRes', 'mhcNum', 'mhcName', 'distance', 'loop']]
    
    if first_beta:
        final_df_beta = []
    else:
        final_df_beta = pandas.DataFrame(pre_df_beta)
        if np.shape(final_df_beta)[1] != 8:
            final_df_beta = np.transpose(final_df_beta)
        final_df_beta.columns = [['tcrRes', 'tcrNum', 'tcrName', 'mhcRes', 'mhcNum', 'mhcName', 'distance', 'loop']]
    
    return(final_df_alpha,final_df_beta)
