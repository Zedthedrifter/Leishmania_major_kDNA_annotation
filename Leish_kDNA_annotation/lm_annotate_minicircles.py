import re
import json
import gzip
import pickle
import pathlib
import datetime
import warnings
from copy import deepcopy
from collections import Counter, OrderedDict
from operator import itemgetter
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
from Bio.SeqFeature import SeqFeature, FeatureLocation, Reference
# warnings.filterwarnings('ignore')
plt.style.use('classic')
from .common import *

def get_features(features, no_maxi=True, mRNA=None):
    for mO_name, feature_list in features.items():
        if no_maxi and mO_name == 'Maxicircle':
            continue
        for feature in feature_list:
            if mRNA is None or feature['mRNA_name'] == mRNA:
                yield feature, mO_name

def annotate(minicircles, CSB1, CSB2, CSB3, aligned_gRNAs,species):
    cstrand = {'coding':1, 'template':-1}
    for mO_name, minicircle_record in minicircles.items():
        minicircle_record.description = species
        minicircle_record.annotations["molecule_type"] = "DNA"
        # minicircle_record.annotations['accession'] = ''
        # minicircle_record.annotations['version'] = ''
        # minicircle_record.annotations['keywords'] = ['kinetoplast DNA', 'kDNA', 'guide RNA', 'gRNA']
        # minicircle_record.annotations['source'] = 'kinetoplast Trypanosoma'
        # minicircle_record.annotations['organism'] = 'Leishmania peruviana'
        # minicircle_record.annotations['topology'] = 'circular'
        # minicircle_record.annotations['date'] = str(datetime.datetime.now().strftime("%d-%b-%Y")).upper()
        # minicircle_record.annotations['taxonomy'] = ['Trypanosomatid', 'etc']
        # ref = Reference()
        # ref.authors = ''
        # ref.title = 'kDNA genome Leishmania peruviana'
        # ref.journal = 'Unpublished'
        # minicircle_record.annotations['references'] = [ref]
        features  = [SeqFeature(FeatureLocation(CSB3[mO_name]['start'], CSB3[mO_name]['end']), type='CSB3')]
        if mO_name in CSB1:
            features += [SeqFeature(FeatureLocation(CSB1[mO_name]['start'], CSB1[mO_name]['end']), type='CSB1')]
        if mO_name in CSB2:
            features += [SeqFeature(FeatureLocation(CSB2[mO_name]['start'], CSB2[mO_name]['end']), type='CSB2')]
#
        if mO_name in aligned_gRNAs:
            features += [SeqFeature(FeatureLocation(int(i['gRNA_start']), int(i['gRNA_end']), strand=cstrand[i['strand']]),
                type='aligned_gRNA',
                qualifiers=OrderedDict([
                    ('name',              i['name']),
                    ('mRNA_',             "5'-{}-3'".format(i['mRNA_seq'])),
                    ('align',             "   {}   ".format(i['pairing'])),
                    ('seq__',             "3'-{}-5'".format(i['gRNA_seq'])),
                    ('length',            i['length']),
                    ('mRNA_start',        i['mRNA_start']+1),
                    ('mRNA_end',          i['mRNA_end']),
                    ('product',           i['mRNA_name']),
                    ('strand',            i['strand']),
                ])) for i in aligned_gRNAs[mO_name]]
#
        minicircle_record.features = sorted(features, key=lambda x: x.location.start)

def format_gRNAs_by_mRNA(gRNAs_by_mO, mRNAs):
    # return a copy of all gRNAs sorted by mRNA that have meet all the filter criteria
#
    gRNAs_by_mRNA = dict([(i, []) for i in mRNAs])
    for gRNA, mO_name in get_features(gRNAs_by_mO, no_maxi=False):
        # if mO_name == 'Maxicircle':
            try:
                gRNAs_by_mRNA[gRNA['mRNA_name']].append(gRNA)
            except KeyError:
                try:
                    gRNAs_by_mRNA[gRNA['mRNA_name']+'_v1'].append(gRNA)
                except:
                    print('here')
                    print(gRNA)
                    exit()
                gRNAs_by_mRNA[gRNA['mRNA_name']+'_v2'].append(gRNA)
    # sort by gapped length and mRNA_end
    for mRNA_name in gRNAs_by_mRNA:
        g = sorted(gRNAs_by_mRNA[mRNA_name], key=itemgetter('gapped_length'))
        gRNAs_by_mRNA[mRNA_name] = sorted(g, key=itemgetter('mRNA_end'), reverse=True)
    return gRNAs_by_mRNA

def output_edits(directory, suffix, mRNAs, aligned_gRNAs, editing_groups):
    gRNAs = format_gRNAs_by_mRNA(aligned_gRNAs, mRNAs)
    pathlib.Path(directory).mkdir(parents=True, exist_ok=True) 
    total_countM = 0
    for mRNA_name, mRNA_record in mRNAs.items():
        if mRNA_name == 'COX2':
            continue
        mRNA_seq = mRNA_record['seq']
        # pack reads into rows
        row = 0
        alignments = [[]]
        # allow gRNAs to extend past the 3' end of a mRNA (to make small RNAs work)
        full_length = mRNA_record['length']+10
        rightmost_limit = full_length
        nrows = 0
        rightmost = {}
        for group in editing_groups[mRNA_name][::-1]:
            #print(group)
            # try putting gRNAs in topmost rows
            row = 0
            while True:
                for gRNA in group['gRNAs']:
                    #print(gRNA)
                    #row=0
                    # check if row exists, if not create it
                    if nrows == row:
                        nrows += 1
                        if nrows > 100:
                            print('Too many rows in edit alignment: {}'.format(mRNA_name))
                            exit()
                        rightmost[row] = full_length
                    row += 1
                    # a gRNA overlaps a previous one, increment row and start again
                    if gRNA['mRNA_end'] >= rightmost[row-1]:
                        #
                        break
                else:
                    # all gRNAs fit so update alignments
                    break
                

            #row -= len(group['gRNAs'])
            for gRNA in group['gRNAs']:
                row=0 #used to be 0
                while rightmost[row]<=gRNA['mRNA_end']: #to avoid overlapping, rightmost is rightmost 5' end position
                    row+=1
                if len(alignments) <= row: #to make a new row (list) to collect gRNAs
                    alignments.append([])
                rightmost[row] = gRNA['mRNA_end']-gRNA['gapped_length'] #update the rightmost 5' of that row
                alignments[row].insert(0, gRNA) #populate the row (list) with gRNA
#
        out = []
        for j in [1000, 100, 10, 1]:
            out.append(''.join([str((i//j)%10) if i%j == 0 else ' ' for i in range(1, len(mRNA_seq)+1)]))
        out.append(''.join(['A' if i == 1 else 'I' if i == 2 else 'U' if i == 3 else '-' for i in mRNA_record['anchour']]))
        # out.append(''.join([str(int(i)) for i in mRNA_record['anchour_count']]))
        out.append(''.join(['-' if i or j != 'u' else 'M' for i, j in zip(mRNA_record['edited'], mRNA_seq)]))
        #print(mRNA_record['edited'], mRNA_seq)#
        countM =  out[-1].count('M') #count missing editing sites
        print(mRNA_name, countM)
        total_countM+=countM
        out.append(''.join(['E' if i else '-' for i, j in zip(mRNA_record['edited'], mRNA_seq)]))
        out.append(mRNA_record['deletions'])
        out.append(mRNA_seq)
        out.append(' '*mRNA_record['orf']+''.join(['{}  '.format(i) for i in mRNA_record['translate']]))
        for row in alignments:
            gRNA_name_align = [' ' for i in range(full_length)]
            pairing_align   = [' ' for i in range(full_length)]
            sequence_align  = [' ' for i in range(full_length)]
            expression      = [' ' for i in range(full_length)]
            for i in range(0, full_length, 10):
                gRNA_name_align[i] = '.'
            for gRNA in row:
                start = gRNA['mRNA_end']-gRNA['gapped_length']
                end   = gRNA['mRNA_end']
#
                if gRNA['mO_name'] == 'Maxicircle':
                    a = gRNA['anchour']
                    a_end = end
                else:
                    a = gRNA['anchour']
                    a_end = gRNA['mRNA_end'] # use end of aligned gRNA to position anchor correctly
                # if gRNA['cassette_pos'] == 'Maxi':
                #     gRNA['cassette_pos'] = ''
#
                if a is None:
                    info = [gRNA['mO_name']]
                    # info = [gRNA['mO_name'], gRNA['cassette_pos'], str(gRNA['group_no'])]
                else:
                    info = [gRNA['mO_name'], '{}{}'.format(a[2]*(a[1]-a[0]), ' '*a[0])]
                    # info = [gRNA['mO_name'], gRNA['cassette_pos'], str(gRNA['group_no']), '{}{}'.format(a[2]*(a[1]-a[0]), ' '*a[0])]
                gRNA_header = ' '.join(info)
                gRNA_name_align[a_end-len(gRNA_header):a_end] = list(gRNA_header)
                pairing_align[start:end]  = list(gRNA['pairing' ][:])
                sequence_align[start:end] = list(gRNA['gRNA_seq'][:])
#
            out.append(''.join(gRNA_name_align))
            out.append(''.join(pairing_align))
            out.append(''.join(sequence_align))
        #print(out)
#
        with open('{}/{}_{}.txt'.format(directory, mRNA_name, suffix), 'w') as f:
            outs = '\n'.join(out)
            f.write(outs)
            # plt.figure(figsize=(20, 9.5))
            # axes = plt.subplot(111)
            # axes.text(0, 0, outs, family='monospace', size=6)
            # plt.show()
    print(f'total number of missing editing sites {total_countM}')

def output_minicircles(directory, suffix, minicircle_list):
    pathlib.Path(directory).mkdir(parents=True, exist_ok=True) 
    SeqIO.write(minicircle_list, f'{directory}/annotated_mOs_{suffix}.gbk', 'genbank')

##########################################################################################
#####put to include the config file
def main(config_file='config.yaml'):
    date = datetime.datetime.now().strftime("%Y-%m-%d")
    ############################################### FILES #########################################
    config = load_config(config_file)
    work_dir, annotation_dir = get_directories(config)[1:3]
    feature_pickled      = f"{work_dir}/{config['features pickle file']}"
    species              = config['species']
    surfix               = config['surfix']
    gRNAs_text_file       = f"{annotation_dir}/gRNAs_{date}.txt"
    
    print('done loading variables')
    #load pickled file
    with gzip.open(feature_pickled, 'rb') as f:
        minicircles = pickle.load(f)
        mRNAs = pickle.load(f)
        CSB1 = pickle.load(f)
        CSB2 = pickle.load(f)
        CSB3 = pickle.load(f)
        aligned_gRNAs = pickle.load(f)
        editing_groups = pickle.load(f)
    #print(mRNAs)#
    annotate(minicircles, CSB1, CSB2, CSB3, aligned_gRNAs,species)
    minicircle_list = [v for k, v in sorted(minicircles.items())]
    output_minicircles(f'{annotation_dir}/Genbank', f"{surfix}_{date}", minicircle_list)
    output_edits(f"{annotation_dir}/Alignments/{surfix}_{date}", f"", mRNAs, aligned_gRNAs, editing_groups)
    #output_edits_2(aligned_gRNAs, mRNAs, config, f"{annotation_dir}/Alignments")
    #make a dataframe of gRNAs from aligned_gRNA dictionary
    idx,tmp_dict=0,{}
    for k1 in aligned_gRNAs:
    #  print(aligned_gRNAs[k1])
      for d in aligned_gRNAs[k1]:
        tmp_dict[idx]=d
        idx+=1
    #print(tmp_dict)
    dataframe_out(pd.DataFrame.from_dict(tmp_dict, orient='index'), gRNAs_text_file) #make an output for gRNAs
#################################################