import re
import gzip
import pickle
import datetime
import warnings
import numpy as np
import pandas as pd
from collections import OrderedDict
from pprint import pprint
from copy import deepcopy
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
#from Bio.Alphabet import generic_dna
#from Bio.Alphabet.IUPAC import unambiguous_rna
from collections import Counter
warnings.filterwarnings('ignore')
from .common import *

def get_features(features, no_maxi=True, mRNA=None):
    for mO_name, feature_list in features.items():
        if no_maxi and mO_name == 'Maxicircle':
            continue
        for feature in feature_list:
            if mRNA is None or feature['mRNA_name'] == mRNA:
                yield feature, mO_name

def complement(seq):
    conv = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N', '-':'-'}
    return ''.join([conv[i] for i in seq])

def pairs(mRNA, gRNA):
    pairings = {'GC':'|', 'CG':'|', 'AU':'|', 'UA':'|', 'GU':':', 'UG':':'}
    return ''.join([pairings[mb+gb] if mb+gb in pairings else '-' if '-' in mb+gb else '.' for mb, gb in zip(mRNA, gRNA)])

def get_maxicircle(filename, file_format='fasta'):
    return SeqIO.read(filename, file_format)

def get_minicircles(filename, file_format='fasta'):
    return SeqIO.to_dict(SeqIO.parse(filename, file_format))

def get_mRNAs(insertion_file, deletion_file):
    #edits:   ACGTt
    #DNA_seq: ACGT
    #seq:     ACGUu
#
    deletions = {}
    with open(deletion_file) as f:
        for line in f:
            if line.startswith('>'):
                mRNA_name = line[1:-1]
            else:
                deletions[mRNA_name] = line[:-1]
#
    insertions = {}
    with open(insertion_file) as f:
        for line in f:
            if line.startswith('>'):
                mRNA_name = line[1:-1]
            else:
                insertions[mRNA_name] = line[:-1]
#
    mRNAs = {}
    for name, edits in sorted(insertions.items()):
        #determine open reading frame as giving the the longest protein sequence between 2 stop codons
        mRNA = Seq(edits.upper().replace('T', 'U'))#,unambiguous_rna)
        max_lengths = [max([len(j) for j in re.split('\*', str(mRNA[orf:].translate(table=4)))]) for orf in range(3)]
        orf = max_lengths.index(max(max_lengths))
#
        mRNAs[name] = {'orf'      :orf,
                       'seq'      :edits.replace('t', 'u').replace('T', 'U'),
                       'edits'    :edits,
                       'length'   :len(edits),
                       'deletions':deletions[name],
                       'DNA_seq'  :edits.replace('u','').replace('U', 'T'),
                       'translate':mRNA[orf:].translate(table=4)
                       }
    return mRNAs

def identify_CSBs(minicircles,CSB_regexes):
    CSB1_regex = re.compile(CSB_regexes['CSB1'])
    CSB2_regex = re.compile(CSB_regexes['CSB2'])
    CSB3_regex = re.compile(CSB_regexes['CSB3'])
#
    CSB1 = {}
    CSB2 = {}
    CSB3 = {}
#
    # move CSB3 to start of sequence
    for mO_name, minicircle_record in minicircles.items():
        seq = minicircle_record.seq
        match = CSB3_regex.search(str(seq))
        if match:
            start = match.start(0)
            minicircle_record.seq = seq[start:]+seq[:start]
        else:
            print('CSB 3 not found: {}'.format(mO_name))
#
    pos = Counter()
    for mO_name, minicircle_record in minicircles.items():
        seq = str(minicircle_record.seq)
        match = CSB1_regex.search(seq)
        if match:
            CSB1[mO_name] = {'start':match.start(0), 'end':match.end(0)}
            pos[len(seq)-match.start(0)] += 1
        else:
            print('no CSB1: {}'.format(mO_name))
#
        match = CSB2_regex.search(seq)
        if match:
            CSB2[mO_name] = {'start':match.start(0), 'end':match.end(0)}
        else:
            print('no CSB2: {}'.format(mO_name))
#
        match = CSB3_regex.search(seq)
        if match:
            CSB3[mO_name] = {'start':match.start(0), 'end':match.end(0)}
        else:
            print('no CSB3: {}'.format(mO_name))
#
    print(pos)
    print('#CSB1: {}'.format(len(CSB1)))
    print('#CSB2: {}'.format(len(CSB2)))
    print('#CSB3: {}'.format(len(CSB3)))
    return CSB1, CSB2, CSB3

def get_gRNAs(filename, minicircles, maxicircle, mRNAs): #need to twick the hq_gRNA file: just apply no filter
    gRNAs = {}
    with open(filename) as f:
        total_len=len(next(f).rstrip().split())
        for alignment in f:
            r = alignment.rstrip()
            if len(r.split()[1:])!=total_len:
                print('not match',total_len, len(r.split()))
            else:
            #[mO_name, mRNA_name, strand, region, position, length, score] = r.split() #original
                [mO_name, strand, length, circle_start, circle_end, mRNA_name, product, mRNA_start, mRNA_end, mRNA_seq, gRNA_seq, pairing, mismatches, anchor_length] = r.split()[1:]
                #mRNA_seq = next(f).lstrip().rstrip()
                #pairing  = next(f).lstrip().rstrip()
                #gRNA_seq = next(f).lstrip().rstrip()
                #next(f)
    #
                gDNA_seq = gRNA_seq.replace('U', 'T').replace('-', '')[::-1]
                #length = len(gDNA_seq)
    #
                mRNA_record = mRNAs[mRNA_name]
                mDNA_seq = mRNA_seq.replace('U', 'T').replace('u', 'T').replace('-', '')
                #mRNA_pos = mRNA_record['DNA_seq'].index(mRNA_seq)
    #
                #new = mRNA_record['edits'][mRNA_pos:mRNA_pos+len(mDNA_seq)].replace('t', 'u').replace('T', 'U')
                new=mRNA_seq
                for match in re.finditer('-', mRNA_seq):
                    pos = match.start(0)
                    new = new[:pos]+'-'+new[pos:]
    #
                aligned_gRNA = {}
                #if strand == 'template':
                #    if 'Maxi' in mO_name:
                #        aligned_gRNA['gRNA_start'] = str(maxicircle.seq).index(reverse_complement(gDNA_seq))
                #    else:
                #        aligned_gRNA['gRNA_start'] = str(minicircles[mO_name].seq).index(reverse_complement(gDNA_seq))
                #else:
                #    if 'Maxi' in mO_name:
                #        aligned_gRNA['gRNA_start'] = str(maxicircle.seq).index(gDNA_seq)
                #    else:
                #        try:
                #            aligned_gRNA['gRNA_start'] = str(minicircles[mO_name].seq).index(gDNA_seq)
                #        except:
                #            print(str(minicircles[mO_name].seq))
                #            print(gRNA_seq)
                #            print(gDNA_seq)
                #            print([mO_name, mRNA_name, strand, position, length])
                #            exit()
                #aligned_gRNA['gRNA_end']      = aligned_gRNA['gRNA_start']+length
                aligned_gRNA['gRNA_start']    = int(circle_start)
                aligned_gRNA['gRNA_end']      = int(circle_end)
                aligned_gRNA['mO_name']       = mO_name
                aligned_gRNA['mRNA_name']     = mRNA_name
                aligned_gRNA['mRNA_start']    = int(mRNA_start)
                aligned_gRNA['mRNA_end']      = int(mRNA_end)
                aligned_gRNA['mRNA_seq']      = mRNA_seq
                aligned_gRNA['length']        = length
                aligned_gRNA['strand']        = strand
                aligned_gRNA['gapped_length'] = len(gRNA_seq)
                aligned_gRNA['gRNA_seq']      = gRNA_seq
                aligned_gRNA['pairing']       = pairing
                aligned_gRNA['mismatches']    = aligned_gRNA['pairing'].count('.')
                aligned_gRNA['method']        = 'bias'
                aligned_gRNA['type']          = 'canonical'
                aligned_gRNA['name']          = '' # name is set in remove_false_positives due to changes in position
                if mO_name not in gRNAs:
                    gRNAs[mO_name] = []
                gRNAs[mO_name].append(aligned_gRNA)
    return gRNAs

def format_gRNAs_by_mRNA(gRNAs_by_mO, mRNAs):
    # return a copy of all gRNAs sorted by mRNA that have meet all the filter criteria
#
    gRNAs_by_mRNA = dict([(i, []) for i in mRNAs])
    for gRNA, mO_name in get_features(gRNAs_by_mO, no_maxi=False):
        #print(mO_name)#
        try:
            gRNAs_by_mRNA[gRNA['mRNA_name']].append(gRNA)
        except KeyError:
            gRNAs_by_mRNA[gRNA['mRNA_name']+'_v1'].append(gRNA)
            gRNAs_by_mRNA[gRNA['mRNA_name']+'_v2'].append(gRNA)
#
    # sort by gapped length and mRNA_end
    for mRNA_name in gRNAs_by_mRNA:
        g = sorted(gRNAs_by_mRNA[mRNA_name], key=itemgetter('gapped_length'))
        gRNAs_by_mRNA[mRNA_name] = sorted(g, key=itemgetter('mRNA_end'), reverse=True)
    return gRNAs_by_mRNA

def identify_anchours(aligned_gRNAs, mRNAs):
    # For each gRNA we idenitfy its maximum anchour as the longest, 5’-most 
    # stretch of Watson-Crick basepairs. An anchour is considered an extender 
    # if there are any insertions or deletions in its first 6 basepairs, an 
    # initiator if there are no insertions or deletions in its first 6 basepairs.
#
    anchour6 = re.compile('\|{6,}')
    anchour5 = re.compile('\|{5}:\|{1,}')
    anchour4 = re.compile('\|{4}:\|{2,}')
    anchour3 = re.compile('\|{3}:\|{3,}')
    anchour2 = re.compile('\|{2}:\|{4,}')
    anchour1 = re.compile('\|{1}:\|{5,}')
#
    # anchour_regex = re.compile('\|{{{},}}'.format(4))
    guiding_regex = re.compile('[\|!:\.]+')
    gRNAs_by_mRNA = format_gRNAs_by_mRNA(aligned_gRNAs, mRNAs) #aligned_gRNAs is the output of get_gRNAs
    a_char = {1:'_', 2:':', 3:'.', 4:'_'}
#
    for mRNA_name, mRNA in mRNAs.items():
        # track prior editing and anchours on the mRNA
        mRNA['anchour_count'] = np.zeros(len(mRNA['seq']))
        mRNA['anchour'] = np.zeros(len(mRNA['seq']))
        mRNA['edited']  = np.zeros(len(mRNA['seq']))
#
        g = sorted(gRNAs_by_mRNA[mRNA_name], key=itemgetter('gapped_length'))
        g = sorted(g, key=itemgetter('mRNA_end'), reverse=True)
        print(mRNA_name,'number of gRNAs',len(g)) #number of gRNAs
#       
        for gRNA in g:
            # find position of minimal anchour on gRNA and the mRNA sequence of this anchour
            match = anchour6.match(gRNA['pairing'][::-1])
            if match is None:
                match = anchour5.match(gRNA['pairing'][::-1])
                if match is None:
                    match = anchour4.match(gRNA['pairing'][::-1])
                    if match is None:
                    #    print('no anchour found')
                    # else:
                    #     print(gRNA['mO_name'])
                         match = anchour3.match(gRNA['pairing'][::-1])
                         if match is None:
                             match = anchour2.match(gRNA['pairing'][::-1])
                             if match is None:
                                 match = anchour1.match(gRNA['pairing'][::-1])
                                 #if match is None:
                                 #    print('no anchour found')
                                 #else:
                                 #    print(gRNA['mO_name'])
            # find position of anchour in alignment pairing
            if match is None:
                #print('no anchour found')
                a_start = 0
                a_end   = 0
            else:
                #print(gRNA['mO_name'])
                a_start = match.start(0)
                a_end   = match.end(0)
            # maximal anchour on mRNA
            max_a_slice = slice(gRNA['mRNA_end']-a_end, gRNA['mRNA_end']-a_start)
            # minimal anchour on mRNA (ie min_anchour_length)
            min_a_slice = slice(gRNA['mRNA_end']-a_start-6, gRNA['mRNA_end']-a_start)
            # get mRNA sequence of min and max anchours
            min_seq  = mRNA['seq'][min_a_slice]
            max_seq  = mRNA['seq'][max_a_slice]
            min_prior_edits = mRNA['edited'][min_a_slice]
            max_prior_edits = mRNA['edited'][max_a_slice]
#
            # determine the type of anchour of this gRNA
            if 'u' not in min_seq:
                if np.sum(min_prior_edits) == 0:
                    a_type = 'initiator'
                else:
                    a_type = 'extenderB' # based on min_anchour_length could be an initiator or an extender
            else:
                for m, e in zip(min_seq, min_prior_edits):
                    if m == 'u' and e == 0:
                        a_type = 'unanchoured'
                        break
                else:
                    a_type = 'extenderA'
#
            # update position of end of anchour based on anchour type and position of prior edits
            if a_type == 'extenderA':
                # Normal anchour using minimum WC base-pairing
                # ie anchour is created by prior insertions
                # pos = np.where(mRNA['edited'][max_a_slice] > 0)[0][0]
                # find 5'-most position of last insertion of prior gRNA edits 
                for p, (m, e) in enumerate(zip(max_seq[::-1], max_prior_edits[::-1])):
                    if e == 0 and m == 'u':
                        pos = len(max_seq)-p
                        a_end -= pos
                        break
                a_value = 1
            elif a_type == 'extenderB':
                # Normal anchour using minimum WC base-pairing
                # ie anchour is created by prior insertions
                # pos = np.where(mRNA['edited'][max_a_slice] > 0)[0][0]
                # find 5'-most position of last insertion of prior gRNA edits 
                for p, (m, e) in enumerate(zip(max_seq[::-1], max_prior_edits[::-1])):
                    if e == 0 and m == 'u':
                        pos = len(max_seq)-p
                        a_end -= pos
                        break
                a_value = 4
            elif a_type == 'initiator':
                # Initiator anchour of at least miniumum length
                # trim maximum anchour to first 'u' if it exists
                match = re.search('\d', mRNA['deletions'][max_a_slice][::-1])
                # search for deletions in anchour (this only seems to happen in A6_v1)
                if match:
                    dpos = match.start(0)
                else:
                    dpos = len(max_seq)
                # find first insertion in anchour
                epos = max_seq[::-1].find('u')
                if epos == -1:
                    epos = len(max_seq)
                # take the minimum position of these two
                pos = min(dpos, epos)
                a_end -= len(max_seq)-pos
                a_value = 2
            else:
                # otherwise gRNA has no edits to anchour it, unattached anchour
                # in analysis.py this is considered an extender
                # position end of anchour at minimum anchour length away form start of anchour
                a_end = a_start+6
                a_value = 3
#
            # add anchour region to mRNA and gRNA
            a_slice = slice(gRNA['mRNA_end']-a_end, gRNA['mRNA_end']-a_start)
            # a_slice = slice(gRNA['mRNA_end']-a_start-4, gRNA['mRNA_end']-a_start)
            mRNA['anchour_count'][a_slice] += 1
            x = mRNA['anchour'][a_slice]
            np.place(x, (x==0) | (x==3), a_value)
            gRNA['anchour'] = a_start, a_end, a_char[a_value]
#
            # add edited region to mRNA
            match = guiding_regex.match(gRNA['pairing'][::-1][a_end:])
            if match:
                e_start = match.start(0)+a_end
                e_end   = match.end(0)+a_end
                mRNA['edited'][gRNA['mRNA_end']-e_end:gRNA['mRNA_end']-e_start] += 1
        mRNA['anchour_count'] = np.clip(mRNA['anchour_count'], 0, 9) #value <0 becomes 0, value > 9 becomes 9

def identify_editing_groups(mRNAs, gRNAs): #causing some gRNAs missing
    def assign_group_no(gRNAs, group_no):
        for gRNA in gRNAs:
            gRNA['group_no'] = group_no
        group_no += 1
        return group_no
#
    gRNAs_by_mRNA = format_gRNAs_by_mRNA(gRNAs, mRNAs)
    editing_groups = dict([(i, []) for i in mRNAs])
    #print(gRNAs_by_mRNA['CYB'][0])
    #print(gRNAs_by_mRNA['ND7'])
#
    for mRNA_name, mRNA in sorted(mRNAs.items()):
#        if mRNA_name=='ND7':
#            print('ND7 is found')
#            print(''.join([str(int(i)) for i in mRNA['anchour_count']])) #all 0 for ND7, not a single anchor
        group_no = 0
        for m in re.finditer('[1-9]+', ''.join([str(int(i)) for i in mRNA['anchour_count']])):
            s, e = m.start(0), m.end(0)
            #gRNAs = [i for i in gRNAs_by_mRNA[mRNA_name] if s <= i['mRNA_end'] and i['mRNA_end'] <= e]
            gRNAs = [i for i in gRNAs_by_mRNA[mRNA_name]]
            if mRNA_name=='ND7':
                print(gRNAs)
#
            # for each cassette position of these gRNAs create a dictionary of cassette position and editing position
            cas_pos = {}
            for gRNA in gRNAs:
                pos = 'I'
                # pos = gRNA['cassette_pos']
                if pos not in cas_pos:
                    cas_pos[pos] = gRNA['mRNA_end']
                cas_pos[pos] = max(gRNA['mRNA_end'], cas_pos[pos])
            # group gRNAs with the same cassette position ordered by editing position
            for pos, v in sorted(cas_pos.items(), key=itemgetter(1)):
                # g = [i for i in gRNAs if i['cassette_pos'] == pos]
                g = gRNAs
                group_no = assign_group_no(g, group_no)
                editing_groups[mRNA_name].append({'start':v, 'end':e, 'gRNAs':g})
    new={k:[editing_groups[k][0]] for k in editing_groups if len(editing_groups[k]) > 0}
    for k in editing_groups:
        if len(editing_groups[k]) == 0:
            new[k]=editing_groups[k]
    return new
#new code not finished yet

def collapse(gRNAs,filter):
  #print(gRNAs)
  pass_gRNA=filter['pass_gRNA']
  start=filter['start_position']
  end=filter['end_position']
  collapsed={}
  for circle in gRNAs:
    if 'MAXI' in circle.upper() or circle in filter['pass_all']:
      collapsed[circle]=gRNAs[circle]
    else:
      collect=[]
      #tmp=[]
      #for gRNA in gRNAs[circle]:
      #  if int(gRNA['gRNA_end'])< start or int(gRNA['gRNA_start'])> end or int(gRNA['length'])>= 35:
      #    collect.append(gRNA)
      #  else:
      #    tmp.append(gRNA)
      #tmp=[gRNA for gRNA in gRNAs[circle]]
      tmp=gRNAs[circle]
      if len(tmp)!=0:
        maxlen=max([gRNA['length'] for gRNA in tmp])
        out=[gRNA for gRNA in tmp if gRNA['length'] == maxlen]
        collect+=out
      #select specific gRNAs
      if circle in pass_gRNA.keys():
        tmp_dict={f"{gRNA['mO_name']}_{gRNA['mRNA_name']}_{gRNA['mRNA_start']}":gRNA for gRNA in tmp}
        out=[tmp_dict[k] for k in tmp_dict if k in pass_gRNA[circle]]
        #print(out)
        collect+=out
      collapsed[circle]=collect
  return(collapsed)
#another way of collapsing: remove redundant gRNAs
def print_number_of_gRNAs(aligned_gRNAs, mRNAs):
    gRNAs_by_mRNA = format_gRNAs_by_mRNA(aligned_gRNAs, mRNAs) #aligned_gRNAs is the output of get_gRNAs
    for mRNA_name, mRNA in mRNAs.items():
        g = sorted(gRNAs_by_mRNA[mRNA_name], key=itemgetter('gapped_length'))
        g = sorted(g, key=itemgetter('mRNA_end'), reverse=True)
        print(mRNA_name,'number of gRNAs',len(g)) #number of gRNAs
##########################################################################################

#####put to include the config file
def main(config_file='config.yaml'):
    ############################################### FILES #########################################
    config = load_config(config_file)
    work_dir = get_directories(config)[1]

    minicircle_clean_file = f"{work_dir}/{config['minicircle clean fasta file']}"
    maxicircle_clean_file = f"{work_dir}/{config['maxicircle clean fasta file']}"
    hq_gRNAs_text_file   = f"{work_dir}/{config['high quality gRNAs text file']}"
    mini_align_file      = f"{work_dir}/{config['minicircle alignments file']}"
    maxi_align_file      = f"{work_dir}/{config['maxicircle alignments file']}"
    edited_mRNA_t_file   = f"{work_dir}/{config['edited mRNA with t fasta file']}"
    deletion_mRNA_file   = f"{work_dir}/{config['deletions mRNA text file']}"
    filter = config['collapse filter']
    ###new lines
    CSB_regexes = config['CSB regexes']
    print(CSB_regexes)
    edited_mRNA_small_u  = f"{work_dir}/{config['edited mRNA with u fasta file']}"
    edited_mRNA_small_t  = f"{work_dir}/{config['edited mRNA with t fasta file']}"
    feature_pickled      = f"{work_dir}/{config['features pickle file']}"
    maxicircle_gRNAs_xlsx= f"{work_dir}/{config['maxicircle_gRNAs_xlsx']}"
    print('done loading variables')
    ########################################## GET MRNAS #########################################
    mRNAs = get_mRNAs(edited_mRNA_small_t, deletion_mRNA_file)
    #print(mRNAs)
    maxicircle = get_maxicircle(maxicircle_clean_file)
    minicircles  = get_minicircles(minicircle_clean_file)
    CSB1, CSB2, CSB3 = identify_CSBs(minicircles,CSB_regexes)
    gRNAs = get_gRNAs(hq_gRNAs_text_file, minicircles, maxicircle, mRNAs)
    gRNAs=collapse(gRNAs,filter)
    identify_anchours(gRNAs, mRNAs)
    #print_number_of_gRNAs(gRNAs, mRNAs)
    editing_groups = identify_editing_groups(mRNAs, gRNAs)
    #for k in editing_groups:
    #  print(k,len(editing_groups[k]))
    #print(mRNAs)
    with gzip.open(feature_pickled, 'wb') as f:
            pickle.dump(minicircles, f)
            pickle.dump(mRNAs, f)
            pickle.dump(CSB1, f)
            pickle.dump(CSB2, f)
            pickle.dump(CSB3, f)
            pickle.dump(gRNAs, f)
            pickle.dump(editing_groups, f)
    ########################################## MAXICIRCLE #########################################
    #if 'Maxicircle' in gRNAs:
    #    maxi_gRNAs = gRNAs['Maxicircle'] #what's this?
    #    g = {}
     #   columns = ['strain', 'start on maxi', 'end on maxi', 'strand', 'gene', 'start on gene', 'end on gene']
     #   a = OrderedDict([(c, []) for c in columns])
    #    for gRNA in sorted(maxi_gRNAs, key=itemgetter('gRNA_start')):
    #        a['strain'].append(code)
    #        a['start on maxi'].append(gRNA['gRNA_start']) 
    #        a['end on maxi'].append(gRNA['gRNA_end']) 
    #        a['strand'].append(gRNA['strand']) 
     #       a['gene'].append(gRNA['mRNA_name']) 
     #       a['start on gene'].append(gRNA['mRNA_start']) 
    #        a['end on gene'].append(gRNA['mRNA_end']) 
    #    a_df = pd.DataFrame(a)
    #    a_df.to_excel('maxicircle_gRNAs.xlsx', index=False)
     #   pprint(g)