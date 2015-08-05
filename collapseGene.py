# abundanceToolkit.py 
"""
This toolkit has several functions that help compare RNA seq runs.
"""


import csv
import os
import json
import numpy as np



####

def abundance_read(infile):
  readline, lineNum = [], -1
  with open(infile, 'r') as fIn:
    for line in fIn:
      if line:
        lineNum = lineNum + 1
        if lineNum > 0 and float(line.split(None)[-1]) != 0:
          readline.append(line.split(None))
  return readline



##################

import os
#
def merge_abundances(filelist, output='testOutput.txt', returnoutput=False):
  # Given a directory, this script loads all *.txt in a given folder and
  # treats them as abundances -- no other .txt files should be in that folder
  if filelist is None:
    fils = ['output_ERR852089.fastq/abundance.txt', 
                'output_ERR852099.fastq/abundance.txt']
  else:
    fils = os.listdir(filelist)
    fils = [f for f in fils if f.split('.')[-1] == 'txt']
  o_data = []
  for afile in fils:
    file_h = open(afile)
    a_list = []
    a_list.append(afile)
    csv_reader = csv.reader(file_h, delimiter='\t')
    for row in csv_reader:
      a_list.append(row[4])
    o_data.append((n for n in a_list))
    file_h.close()
  #
  with open('testOutput', 'w') as op_file:
    csv_writer = csv.writer(op_file, delimiter='\t')
    for row in list(zip(*o_data)):
      csv_writer.writerow(row)
  op_file.close()
  if returnoutput:
    return o_data

###############

def parse_gname(gi):
  # From abundance file from Kallisto
  # Example: gi|559098400|ref|NM_001287053.1|
  s = gi.split('|')
  return int(s[1]), s[3]



def get_genedict(refgfile='gene2refseq'):
  """
  Return the default gene dictionary.
  """
  genedict = {}
  # Load the reference sequence to genome file
  with open(refgfile, 'r') as fIn:
    for line in fIn:
      if line:
        splitline = line.split(None)
        g_num = int(splitline[1])
        nm_num = splitline[3]
        tx_num = int(splitline[4])
        # Add each gene to the gene dict, then add its transcripts
        if g_num not in genedict.keys():
          genedict[g_num] = {'nm_nums': [], 'gi_nums': []}
        if nm_num not in genedict[g_num]['nm_nums']:
          genedict[g_num]['nm_nums'].append(nm_num)
        if tx_num not in genedict[g_num]['gi_nums']:
          genedict[g_num]['gi_nums'].append(tx_num)
  return genedict



def load_abfile(abfile, countonly=True):
  """
  Load the abundance file -- Kallisto-like input.
  """
  # Should now have all of the transcripts for all the genes
  #
  # Get the matrix from the abfile and assign it to the genes
  lol, lineNum = [], -1
  with open(abfile, 'r') as fIn:
    for line in fIn:
      if line:
        lineNum = lineNum + 1
        splitline = line.split(None)
        if lineNum > 1:
          # But only get the 
          if float(splitline[3]) > 0:
            gname = parse_gname(splitline[0])
            if countonly is False:
              thing = [gname[0], gname[1], int(splitline[1]),
                                              int(splitline[2])]
              for k in [float(i) for i in splitline[3:]]:
                thing.append(k)
            else:
              thing =[gname[0], gname[1]]
              for k in [float(i) for i in splitline[1:]]:
                thing.append(k)
            lol.append(thing)
  return lol



def collapse_genes(abfile, genedict='gene2refseq', countonly=True):
  """
  Given an abundance file (similar to Kallisto format, except colums are:
  target id (gi|...|ref|NM_...|) length eff_length counts0 counts1 ... countsN
  """
  def add_entry(ab_dict, lin, g, countonly=True):
    ab_dict[g]['gi_nums'].append(lin[0])
    ab_dict[g]['nm_nums'].append(lin[1])
    if countonly is False:
      ab_dict[g]['length'].append(lin[2])
      ab_dict[g]['eff_length'].append(lin[3])
      ab_dict[g]['count'].append(lin[4])
    else:
      ab_dict[g]['count'].append(lin[2])
    return ab_dict
    #
  # Get the list of lists (abfile matrix)
  if type(abfile) == str:
    lol = load_abfile(abfile, countonly)
  elif type(abfile) == list:
    lol = abfile
  # lol[i] = [gi(int), nm(str), length(int), eff_len(int), count(float)]
  if type(genedict) == str:
    genedict = get_genedict(genedict)
  elif type(genedict) == dict:
    genedict = genedict
  # Should have all transcript info
  samples = [{} for i in range(len(lol[0])-2)] # Log all samples
  lost = [] # In case nothing fits
  # Now, for each transcript in lol log its counts and length in a new dict
  # called ab_gene
  for l in lol: # For each line (each transcript!)
    for i in range(2,len(l)): # For each sample (each column except first 2)
      for g in genedict.keys(): # Scan through genes
        # If this particular transcript is in the gene dict, log it and counts
        ##print(l)
        #print(l[i])
        if l[0] in genedict[g]['gi_nums'] or l[1] in genedict[g]['nm_nums']:
          # Log this gene
          if g not in samples[i-2].keys():
            if countonly is False:
              samples[i-2][g] = {'gi_nums': [], 'nm_nums': [], 'eff_length': [],
                            'length': [], 'count': []}
            else:
              samples[i-2][g] = {'gi_nums': [], 'nm_nums': [], 'count': []}
          samples[i-2] = add_entry(samples[i-2], l, g, countonly)
        else:
          lost.append(l)
  # Now have all genes assigned
  return np.unique(samples)


def simple_collapse(samples, genedict='gene2refseq', show=None):
  """
  Return simple counts by gene instead of transcipt.
  """
  # Get the most common transcript
  def most_common_tx(dict_elem):
    max_gi, gi_name, max_nm, nm_name = 0, None, 0, None
    for gi in dict_elem['gi_nums']:
      if dict_elem['gi_nums'].count(gi) > max_gi:
        max_gi = dict_elem['gi_nums'].count(gi)
        gi_name = gi
    for nm in dict_elem['nm_nums']:
      if dict_elem['nm_nums'].count(nm) > max_nm:
        max_nm = dict_elem['nm_nums'].count(nm)
        nm_name = nm
    if man_nm > max_gi:
      return nm_name
    else:
      return gi_name
  #
  # Collapse gene-wise stats into single element
  sample_trunc = []
  for s in samples:
    for g in s.keys():
      gene_trunc[g] = {'counts': sum(s[g]['counts']),
                       'eff_length_max': max(s[g]['eff_length']),
                       'most_common': most_common_tx(s[g])}
  if show is None:
    print(sample_trunc)
  else:
    json.dump(sample_trunc, open(show, 'w'))
  return



def compare_variability():
  
  
  return


  




  























