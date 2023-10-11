#!/usr/bin/env python3.8

import Leish_kDNA_annotation as lka #Leish_kDNA_annotation a directory where the .py files are stored
from os import system
import subprocess as sbp
import argparse

def prepare(config_file):
  lka.clean_mini_and_maxicircles(config_file)
  lka.lm_mRNA_process(config_file)
  system(f'align_maxi {config_file}') 
  system(f'align_mini {config_file}')
  

def leish(config_file):
  lka.lm_hq_gRNAs(config_file)
  #lka.lm_filter_alignments(config_file) #quite useless actually but good practice
  lka.lm_feature_id(config_file)
  lka.lm_annotate_minicircles(config_file)
  
  
def main(config_file,work_dir):
  prepare(config_file) #both
  leish(config_file)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='run kDNA annotation')
  parser.add_argument('config_file', help='Name of the config file', metavar='config_file')
  parser.add_argument('work_dir', help='Name of the working directiory', metavar='work_dir')
  options = parser.parse_args()
  
  main(options.config_file, options.work_dir)
  

