import re
import numpy as np
import seaborn as sns
import pandas as pd
import warnings
import matplotlib.pyplot as plt
from pprint import pprint
from sys import exit
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Data import CodonTable
warnings.filterwarnings('ignore')

from .common import *

wcgu_mm = re.compile('([-\.]+)')
wcgu = re.compile('[-]+')
standard_re = re.compile('\|{6}[\|:]{19,}')
# wcgu = re.compile('[-\.]+')
# standard_re = re.compile('[\|:]{19,}\|{6}')
anchourless = re.compile('[\|:]{25,}')

min_gRNA_length = 25
max_gRNA_length = 70

def rand_jitter(arr):
    stdev = .01*(max(arr)-min(arr))
    return arr + np.random.randn(len(arr)) * stdev

def complement(seq):
    conv = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N', '-':'-'}
    return ''.join([conv[i] for i in seq])

def match_length(s):
    return len(s)-(s.count('-')+s.count('.'))

def scoring(s):
    # return len(s)-(s.count('.')+2*s.count('-')), -s.count('-')
    contigs = wcgu.split(s)
    max_contig = max([len(i) for i in contigs])
    # error_rate = (s.count('!')+s.count('I')+s.count('.'))/len(s)
    error_rate = (s.count('!')+s.count('.'))/len(s)
    matches = match_length(s)
    try:
        max_bi_contigs = max([len(i)+len(j) for i, j in zip(contigs[:-1], contigs[1:])])
    except:
        max_bi_contigs = 0
    contiguous = len(wcgu.split(s)[-1])
    return max_contig, matches, error_rate, max_bi_contigs, contiguous

def g_len(new_pairing, gRNA_seq):
    return len(new_pairing)-gRNA_seq[-len(new_pairing):].count('-')
def mm_pairing(mRNA, gRNA):
  pairings = {'GC':'|', 'CG':'|', 'AU':'|', 'UA':'|', 'GU':':', 'UG':':'}
  pairs = []
  for mb, gb in zip(mRNA, gRNA):
    bp = mb.upper()+gb.upper()
    if bp in pairings:
      pairs.append(pairings[bp])
    elif (gb == 'U' or gb == 'C') and mb == 'u':
      pairs.append('!')
    else:
      pairs.append('.')
  return ''.join(pairs)
def get_alignments(infiles, outfile, mRNAs, maxicircle, minicircles):
    """ 
    get raw ouput of minicircle and maxicircle alignments
    add insertions into mRNA and do T->U
    remove alignments that do not include insertions
    remove non-unique alignments
    """
    dd = {'b':'coding', 'n':'template', 'r':'NA', 'c':'NA'}
    unique = {}

    alignments = {}

    for infile in infiles:
        with open(infile) as f:
            for circle_name in f:
                circle    = circle_name.rstrip()[1:]
                mRNA_name, mRNA_seq, gDNA_seq, strand = next(f).rstrip().split()
                dd = {'c':'coding', 't':'template'}
                strand = dd[strand]
                # remove maxicircle position from alignments
                if 'MAXI' in circle.upper():
                    circle = 'Maxicircle'


                # find start position on circle
                # ignore coding strand gRNAs on minicircles
                # gDNA_seq = gRNA_seq[::-1]
                if 'Maxi' in circle:
                    if strand == 'template':
                        gRNA_start = str(maxicircle.seq).index(gDNA_seq)
                    else:
                        gRNA_start = str(maxicircle.seq).index(reverse_complement(gDNA_seq))
                else:
                    if strand == 'template':
                        gRNA_start = str(minicircles[circle].seq).index(gDNA_seq)
                    else:
                        continue
                gRNA_seq = complement(gDNA_seq).replace('T', 'U')
                gRNA_end = gRNA_start + len(gRNA_seq)

                # find position of alignment on edited mRNA
                mRNA_start = mRNAs[mRNA_name]['DNA_seq'].index(mRNA_seq)
                mRNA_end = mRNA_start + len(mRNA_seq)
                mRNA_seq = mRNAs[mRNA_name]['seq'][mRNA_start:mRNA_end]
                pairing   = mm_pairing(mRNA_seq, gRNA_seq)
                # combine overlapping alignments into longer alignments
                # need to do this because in align_minicircle.c and align_maxicircle.c
                # the max gRNA length was set to 50
                key = circle, strand, mRNA_name
                if key not in alignments:
                    # if 
                    a = {}
                    a['circle']     = circle
                    a['mRNA_name']  = mRNA_name
                    a['strand']     = strand
                    a['mRNA_seq']   = mRNA_seq
                    a['gRNA_seq']   = gRNA_seq
                    a['pairing']    = pairing
                    a['gRNA_start'] = gRNA_start
                    a['gRNA_end']   = gRNA_end
                    a['mRNA_start'] = mRNA_start
                    a['mRNA_end']   = mRNA_end
                    alignments[key] = [a]
                else:
                    # check this alignment against all other on this circle and mRNA
                    for a in alignments[key]:
                        # if gRNAs overlap on circle
                        if gRNA_start <= a['gRNA_end'] and a['gRNA_start'] <= gRNA_end:
                            start_overhang = a['gRNA_start']-gRNA_start
                            end_overhang   = a['gRNA_end']-gRNA_end
                            if strand == 'template':
                                # if length of overhang is the same on circle and mRNA then combine
                                if start_overhang == a['mRNA_start']-mRNA_start:
                                    if start_overhang > 0:
                                        a['gRNA_start'] = gRNA_start
                                        a['mRNA_start'] = mRNA_start
                                        a['gRNA_seq']   = gRNA_seq[:start_overhang]+a['gRNA_seq']
                                        a['mRNA_seq']   = mRNA_seq[:start_overhang]+a['mRNA_seq']

                                    if end_overhang < 0:
                                        a['gRNA_end']   = gRNA_end
                                        a['mRNA_end']   = mRNA_end
                                        a['gRNA_seq']   = a['gRNA_seq']+gRNA_seq[end_overhang:]
                                        a['mRNA_seq']   = a['mRNA_seq']+mRNA_seq[end_overhang:]
                                    break
                            elif strand == 'coding':
                                # if length of overhang is the same on circle and mRNA then combine
                                if start_overhang == mRNA_end-a['mRNA_end']:
                                    if start_overhang > 0:
                                        a['gRNA_start'] = gRNA_start
                                        a['mRNA_end']   = mRNA_end
                                        a['gRNA_seq']   = a['gRNA_seq']+gRNA_seq[-start_overhang:]
                                        a['mRNA_seq']   = a['mRNA_seq']+mRNA_seq[-start_overhang:]

                                    if end_overhang < 0:
                                        a['gRNA_end']   = gRNA_end
                                        a['mRNA_start'] = mRNA_start
                                        a['gRNA_seq']   = gRNA_seq[:-end_overhang]+a['gRNA_seq']
                                        a['mRNA_seq']   = mRNA_seq[:-end_overhang]+a['mRNA_seq']
                                    break
                    else:
                        a = {}
                        a['circle']     = circle
                        a['mRNA_name']  = mRNA_name
                        a['strand']     = strand
                        a['mRNA_seq']   = mRNA_seq
                        a['gRNA_seq']   = gRNA_seq
                        a['pairing']    = pairing
                        a['gRNA_start'] = gRNA_start
                        a['gRNA_end']   = gRNA_end
                        a['mRNA_start'] = mRNA_start
                        a['mRNA_end']   = mRNA_end
                        alignments[key].append(a)

    with open(outfile, 'w') as f:
        for key, c in sorted(alignments.items(), key=lambda k_v: k_v[0]):
            for a in c:
                a['pairing'] = mm_pairing(a['mRNA_seq'], a['gRNA_seq'])
                f.write('{} {} {}\n'.format(a['circle'], a['mRNA_name'], a['strand']))
                f.write('{}\n{}\n{}\n'.format(a['mRNA_seq'], a['pairing'], a['gRNA_seq']))

def complement(seq):
    conv = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N', '-':'-'}
    return ''.join([conv[i] for i in seq])

def find_likely_gRNAs(infile, canonical_outfile, noncanonical_outfile, maxi_outfile, 
    maxicircle, minicircles, mRNAs):
    anchour6 = re.compile('\|{6,}')
    anchour5 = re.compile('\|{5}:\|{1,}')
    anchour4 = re.compile('\|{4}:\|{2,}')
    # anchour3 = re.compile('\|{3}:\|{3,}')
    # anchour2 = re.compile('\|{2}:\|{4,}')
    # anchour1 = re.compile('\|{1}:\|{5,}')

    g_start = 450
    g_end = 525
    unique = {}
    for infile in [infile]:
        with open(infile) as f:
            for alignment in f:
                r = alignment.rstrip()
                circle, mRNA_name, strand = r.split()
                a = {}
                a['circle']    = circle
                a['mRNA_name'] = mRNA_name
                a['strand']    = strand
                a['mRNA_seq']  = next(f).rstrip()
                next(f)
                a['gRNA_seq']  = next(f).rstrip()
                # new pairing
                a['pairing']   = mm_pairing(a['mRNA_seq'], a['gRNA_seq'])

                gDNA_seq = a['gRNA_seq'].replace('U', 'T').replace('u', 'T')[::-1]
                try:
                    if 'Maxi' in circle:
                        if strand == 'template':
                            a['gRNA_start'] = str(maxicircle.seq).index(reverse_complement(gDNA_seq))
                        else:
                            a['gRNA_start'] = str(maxicircle.seq).index(gDNA_seq)
                    else:
                        if strand == 'template':
                            a['gRNA_start'] = str(minicircles[circle].seq).index(reverse_complement(gDNA_seq))
                        else:
                            # a['gRNA_start'] = str(minicircles[circle].seq).index(gDNA_seq)
                            continue
                except ValueError:
                    # alignment is done on untransformed minicircles sequences (CSB3 not at position 0)
                    # this can cause the gRNA_seq to overlap the end-start of the transformed sequences
                    # so need to check this. Not an issue because gRNAs don't occur at this position
                    continue

                # remove: 
                # unfiltered maxicircle sequences which aren't gRNAs
                # For minicircles:
                # coding strand
                # without a proper anchour
                # any too short or too long
                # duplicates
                # on minicircle not in correct region
                # does not edit

                # remove unfiltered maxicircle sequences which aren't gRNAs
                match = re.search('\|+', a['pairing'])
                if match and match.group(0) == a['pairing']:
                    continue
                if '|'*23 in a['pairing'] and mRNA_name != 'MURF2':
                    continue

                a['gRNA_end'] = a['gRNA_start']+len(a['pairing'])

                # find anchour
                match_pos = {}
                match6 = anchour6.search(a['pairing'][::-1])
                match5 = anchour5.search(a['pairing'][::-1])
                match4 = anchour4.search(a['pairing'][::-1])
                # match3 = anchour3.search(a['pairing'][::-1])
                # match2 = anchour2.search(a['pairing'][::-1])
                # match1 = anchour1.search(a['pairing'][::-1])
                if match6:
                    match_pos[6] = match6.start(0)
                if match5:
                    match_pos[5] = match5.start(0)
                if match4:
                    match_pos[4] = match4.start(0)
                # if match3:
                #     match_pos[3] = match3.start(0)
                # if match2:
                #     match_pos[2] = match2.start(0)
                # if match1:
                #     match_pos[1] = match1.start(0)
                if len(match_pos):
                    x = [(k, v) for k, v in sorted(match_pos.items(), key=lambda x: x[1])][0]
                    a['anchour_len'] = x[0]
                    a['anchour_pos'] = x[1]
                else:
                    # ignore if no anchour of 4 WC
                    continue

                a['length'] = len(a['pairing'])-a['anchour_pos']
                a['gRNA_end'] = a['gRNA_start']+a['length']
                # trim 5' end to start of anchour
                a['mRNA_seq'] = a['mRNA_seq'][:a['length']]
                # remove alignments that do not edit anything
                pos = str(mRNAs[mRNA_name]['seq']).find(a['mRNA_seq'])
                if ('u' not in mRNAs[mRNA_name]['seq'][pos:pos+a['length']] and
                    re.search(mRNAs[mRNA_name]['deletions'][pos:pos+a['length']], '\d') is None):
                    continue
                a['pairing']   = a['pairing'][:a['length']]
                a['gRNA_seq']  = a['gRNA_seq'][:a['length']]
                a['score']     = scoring(a['pairing'])[2]
                a['chosen']    = False
                a['redundant'] = False

                # Find region where gRNA is on minicircle
                if 'mO' in circle and (a['gRNA_end'] < g_start or a['gRNA_end'] > g_end):
                    a['region'] = 'outside'
                elif 'mO' in circle :
                    a['region'] = 'inside'
                else:
                    a['region'] = 'maxi'

                # too short alignments
                if ('mO' in circle and a['length'] < min_gRNA_length or a['length'] > max_gRNA_length or
                    'Maxi' in circle and a['length'] < 40 or a['length'] > max_gRNA_length):
                    a['size'] = 'too short'
                else:
                    a['size'] = 'correct'

                # find unique gRNAs
                if circle not in unique:
                    unique[circle] = {}
                key = a['mRNA_seq'], a['gRNA_seq'], a['mRNA_name']
                # test for uniqueness
                for k in list(unique[circle].keys()):
                    if key[2] in k[2]:
                        # if this gRNA sequence is in an existing gRNA sequence
                        if key[0] in k[0] and key[1] in k[1]:
                            # if this mRNA sequence and name are in an existing mRNA sequence ignore this alignment
                            break
                        elif k[0] in key[0] and k[1] in key[1]:
                            # if this mRNA sequence and name are are superset of an existing mRNA sequence delete existing alignment and replace with this one
                            del unique[circle][k]
                else: 
                    # new alignment or replace existing
                    unique[circle][key] = a

    # mark redundant alignments
    # put all alignments into a flat list
    all_alignments = []
    for mO_name, alignments in unique.items():
        all_alignments += list(alignments.values())

    # sort inside and maxi alignments by the length
    all_alignments = [a for a in sorted(all_alignments, key=itemgetter('length')) if a['region'] != 'outside' and a['size'] == 'correct']

    # remove redundant gRNAs
    for i in range(len(all_alignments)-1):
        a1 = all_alignments[i]
        for j in range(i+1, len(all_alignments)):
            a2 = all_alignments[j]
            if a1['mRNA_name'] == a2['mRNA_name'] and a1['mRNA_seq'] in a2['mRNA_seq']:
                key = a1['mRNA_seq'], a1['gRNA_seq'], a1['mRNA_name']
                a1['redundant'] = True
                break

    mg = []
    mO_has_gRNA = {}
    with open(canonical_outfile, 'w') as o, open(noncanonical_outfile, 'w') as n:
        for mO_name, alignments in sorted(unique.items(), key=lambda k_v: k_v[0]):
            if 'Maxi' in mO_name:
                for a in alignments.values():
                    if a['size'] == 'correct':
                        a['chosen'] = True
                        x = ' '*(70-len(a['pairing']))
                        o.write('{} {} {} {} {} {} {:.3f}\n'.format(a['circle'], a['mRNA_name'], a['strand'], a['region'], a['gRNA_start'], a['length'], a['score']))
                        o.write('{}{}\n{}{}\n{}{}\n'.format(x, a['mRNA_seq'], x, a['pairing'], x, a['gRNA_seq']))
                        o.write('\n')
                        mg.append(a)
                    else:
                        x = ' '*(70-len(a['pairing']))
                        n.write('{} {} {} {} {} {} {:.3f}\n'.format(a['circle'], a['mRNA_name'], a['strand'], a['region'], a['gRNA_start'], a['length'], a['score']))
                        n.write('{}{}\n{}{}\n{}{}\n'.format(x, a['mRNA_seq'], x, a['pairing'], x, a['gRNA_seq']))
                        n.write('\n')
            else:
                for a in alignments.values():
                    if a['region'] == 'inside' and a['size'] == 'correct' and a['redundant'] == False:
                        a['chosen'] = True
                        x = ' '*(70-len(a['pairing']))
                        o.write('{} {} {} {} {} {} {:.3f}\n'.format(a['circle'], a['mRNA_name'], a['strand'], a['region'], a['gRNA_start'], a['length'], a['score']))
                        o.write('{}{}\n{}{}\n{}{}\n'.format(x, a['mRNA_seq'], x, a['pairing'], x, a['gRNA_seq']))
                        o.write('\n')
                        mO_has_gRNA[mO_name] = True
                    else:
                        x = ' '*(70-len(a['pairing']))
                        n.write('{} {} {} {} {} {} {:.3f}\n'.format(a['circle'], a['mRNA_name'], a['strand'], a['region'], a['gRNA_end'], a['length'], a['score']))
                        n.write('{}{}\n{}{}\n{}{}\n'.format(x, a['mRNA_seq'], x, a['pairing'], x, a['gRNA_seq']))
                        n.write('\n')

    with open(maxi_outfile, 'w') as o:
        for a in sorted(mg, key=itemgetter('gRNA_start')):
            o.write(f"{a['gRNA_start']} {a['gRNA_end']} {a['strand']} {a['mRNA_name']}\n")

    # print(len(mO_has_gRNA))
    # data = {'score':[], 'length':[], 'region':[], 'position':[], 'size':[], 'chosen':[], 'redundant':[]}
    # alignments = [a for key in unique.values() for a in key.values()]
    # for a in sorted(alignments, key=itemgetter('length'), reverse=False):
    #     data['score'].append(100*a['score'])
    #     data['length'].append(a['length'])
    #     data['region'].append(a['region'])
    #     data['size'].append(a['size'])
    #     data['chosen'].append(a['chosen'])
    #     data['redundant'].append(a['redundant'])
    #     data['position'].append(a['gRNA_end'])
    # df = pd.DataFrame(data)
    # df['length_jitter'] = rand_jitter(df.length)
    # df['score_jitter'] = rand_jitter(df.score)
    # df['size_cat'] = df.length.apply(lambda x: '>=40nt' if x >= 40 else '30-39nt' if x >= 30 else '25-29')
    # print(len(df))

    # plt.figure()
    # sns.stripplot('redundant', 'length', data=df.query('region != "outside" and size == "correct"'))

    # plt.figure()
    # ax = sns.scatterplot('length_jitter', 'score_jitter', 'chosen', data=df.query('region == "inside" and size == "correct"'), linewidth=0)
    # ax.set_xlabel('gRNA length (jittered')
    # ax.set_ylabel('% of mismatches between gRNA and mRNA (jittered')
    # ax.set_title('Decision on whether gRNAs on minicircles are chosen or not')
    # plt.savefig('Figures/gRNA_chosen.pdf')

    # plt.figure()
    # ax = sns.scatterplot('length_jitter', 'score_jitter', 'chosen', data=df.query('region == "maxi"'), linewidth=0)
    # ax = sns.scatterplot('length', 'score', 'chosen', data=df.query('region == "maxi"'), linewidth=0)

    # figures for ms
    # plt.figure()
    # ax = sns.scatterplot('position', 'score_jitter', 'size_cat', palette=['C2', 'C1', 'C0'], data=df.query('region != "maxi" and size == "correct"'), linewidth=0)
    # ax.set_xlabel("Distance from 5' end of CSB3 to 5' end of gRNA on template strand")
    # ax.set_ylabel('% of mismatches between gRNA and mRNA')
    # ax.set_title('{}\nShown are all potential gRNAs with a correct anchour.\nHighly probable gRNAs are located within a narrow region on minicircles'.format(code))
    # ax.axvline(g_start)
    # ax.axvline(g_end)
    # plt.savefig(f'Figures/gRNA_postion_on_minicircle_{code}.pdf')


#codes = ('LCA04', 'HR78', 'HR434', 'LC1412')
#strains = ('Lperuviana_LCA04', 'Lperuviana_HR78', 'Lhybrid_HR434', 'Lbraziliensis_LC1412')

#for code, strain in zip(codes, strains):
#    print(code)
#   maxicircle = get_maxicircle('Maxicircle/{}_maxicircle.fasta'.format(strain))
#    minicircles = get_minicircles('Minicircles/{}_minicircles.fasta'.format(strain))
#
#    mRNAs = get_mRNAs('mRNA/edited_mRNA_small_t_{}.fasta'.format(code), 
#        'mRNA/deletion_mRNA_{}.txt'.format(code))
#    get_alignments(('Minicircles/alignments_mini_{}.txt'.format(code), 
#        'Maxicircle/alignments_maxi_{}.txt'.format(code)), 
#        'gRNAs/gRNAs_{}.txt'.format(code), mRNAs, maxicircle, minicircles)
#    find_likely_gRNAs((f'gRNAs/gRNAs_{code}.txt',), f'gRNAs/canonical_{code}.txt',
#        f'gRNAs/noncanonical_{code}.txt', f'gRNAs/maxi_gRNAs_{code}.txt', maxicircle, minicircles, mRNAs)

plt.show()



##
##
##
def main(config_file='config.yaml'):
    ############################################### FILES #########################################
    config = load_config(config_file)
    work_dir = get_directories(config)[1]

    minicircle_file      = f"{work_dir}/{config['minicircle clean fasta file']}"
    maxicircle_clean_file = f"{work_dir}/{config['maxicircle clean fasta file']}"
    mini_align_file      = f"{work_dir}/{config['minicircle alignments file']}"
    maxi_align_file      = f"{work_dir}/{config['maxicircle alignments file']}"
    edited_mRNA_t_file   = f"{work_dir}/{config['edited mRNA with t fasta file']}"
    deletion_mRNA_file   = f"{work_dir}/{config['deletions mRNA text file']}"
    hq_gRNAs_pickle_file = f"{work_dir}/{config['high quality gRNAs pickle file']}"
    hq_gRNAs_text_file   = f"{work_dir}/{config['high quality gRNAs text file']}"
    hq_gRNAs_fasta_file  = f"{work_dir}/{config['high quality gRNAs fasta file']}"
    filtered_gRNAs_text_file  = f"{work_dir}/{config['filtered gRNAs text file']}"
    canonical_outfile    = f"{work_dir}/{config['canonical_outfile']}"
    noncanonical_outfile = f"{work_dir}/{config['noncanonical_outfile']}"
    maxi_outfile         = f"{work_dir}/{config['maxi_outfile']}"

    ########################################## PARAMETERS #########################################
    # parameters defining high quality gRNAs
    filter = config['high quality gRNAs filter']

    # number of nt upstream and downstream of start of gRNA to output for motif calling and nt bias scoring 
    up = config['upstream']  
    down = config['downstream']  


    ####################################### LOAD SEQUENCES #########################################
    minicircles = get_minicircles(minicircle_file)
    maxicircle = get_maxicircle(maxicircle_clean_file)
    mRNAs = get_mRNAs(edited_mRNA_t_file, deletion_mRNA_file)
    get_alignments((mini_align_file,maxi_align_file), filtered_gRNAs_text_file, mRNAs, maxicircle, minicircles)
    #
    find_likely_gRNAs(filtered_gRNAs_text_file, canonical_outfile, noncanonical_outfile, maxi_outfile, 
    maxicircle, minicircles, mRNAs)