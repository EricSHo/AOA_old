# Import packages
import pandas as pd
import re
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

# Assume at the root directory AOA
path = "./data"

for bacterium in ['hInfluenzae','sAuresus','pGingivalis']:
'
  # STEP 1: Set up and read files
  # Set up
  sequence_RNAfold = [] # sequences from RNAfold
  deg_ids_RNAfold = []  # deg_ids from RNAfold
  sequence_mxfold = []  # sequences from mxfold
  deg_ids_mxfold = []   # deg_ids from mxfold
  rna_file = f'{path}/{bacterium}_eg_RNAfold.txt'
  mxfold_file = f'{path}/{bacterium}_eg_mxfold2.txt'
  
  # Read fasta file and store sequences in list 
  for seq_record in SeqIO.parse(open(rna_file),'fasta'): # RNAfold
    # Add record to list
    sequence_RNAfold.append(str(seq_record.seq))
    deg_ids_RNAfold.append(seq_record.id)
  
  for seq_record in SeqIO.parse(open(mxfold_file),'fasta'): # mxfold
    sequence_mxfold.append(str(seq_record.seq))  # sequences from mxfold
    deg_ids_mxfold.append(seq_record.id)
  
  # STEP 2: Separate lines by delimiter
  # Set up
  input_list_RNAfold = [] # Store input sequences RNAfold
  prediction_list_RNAfold = [] # Store secondary predicition RNAfold
  input_list_mxfold = [] # Store input sequences mxfold
  prediction_list_mxfold = [] # Store input sequences mxfold
  
  for i in range(len(sequence_RNAfold)): 
    count = 0
    s = sequence_RNAfold[i]
    if count == 0:
      input_seq = re.split('[.()]', s) # get input sequence
      input_list_RNAfold.append(input_seq[0])
  
      prediction = s[len(input_seq[0]):2*len(input_seq[0])] # get secondary prediction
      prediction_list_RNAfold.append(prediction)
  
      count += 1
  
  for i in range(len(sequence_mxfold)): 
    count = 0
    s = sequence_mxfold[i]
    if count == 0:
      input_seq = re.split('[.()]', s) # get input sequence
      input_list_mxfold.append(input_seq[0])
  
      prediction = s[len(input_seq[0]):2*len(input_seq[0])] # get secondary prediction
      prediction_list_mxfold.append(prediction)
  
      count += 1
  
  # STEP 3: Masking     FIX
  # Set up 
  masked_seq_RNAfold = [] # Store masked sequences
  masked_seq_mxfold = []
  
  for i in range(len(input_list_RNAfold)): 
    s_RNAfold = input_list_RNAfold[i].lower() # Everything starts off masked (aka lower case)
    p_RNAfold = prediction_list_RNAfold[i]
  
    # String to list
    list_sequence = list(s_RNAfold)
  
    # Inner for loop (to iterate through prediction list)
    for k in range(len(s_RNAfold)):
      if p_RNAfold[k] == '.':
        unmask = list_sequence[k].upper()
        list_sequence[k] = unmask
  
    # Append to masked sequence list, convert list back to string
    masked_seq_RNAfold.append(''.join(list_sequence))
  
  for i in range(len(input_list_mxfold)): 
    # Set up
    s_mxfold = input_list_mxfold[i].lower()
    p_mxfold = prediction_list_mxfold[i]
  
    # String to list
    list_sequence = list(s_mxfold)
  
    # Inner for loop (to iterate through prediction list)
    for k in range(len(s_mxfold)):
      if p_mxfold[k] == '.':
        unmask = list_sequence[k].upper()
        list_sequence[k] = unmask
  
    # Append to masked sequence list, convert list back to string
    masked_seq_mxfold.append(''.join(list_sequence))
  
  print(len(masked_seq_RNAfold))
  print(len(masked_seq_mxfold))
  
  # Write to file
  # File for RNAfold
  ofname = f"{path}/{bacterium}_eg_masked_RNAfold.txt"
  ofile = open(ofname, "w")
  
  for i in range(len(masked_seq_RNAfold)):
    ofile.write(">" + deg_ids_RNAfold[i] + "\n" +masked_seq_RNAfold[i] + "\n")
  
  ofile.close()
  
  # File for mxfold
  ofname = f"{path}/{bacterium}_eg_masked_mxfold.txt"
  ofile = open(ofname, "w")
  
  for i in range(len(masked_seq_mxfold)):
    ofile.write(">" + deg_ids_mxfold[i] + "\n" +masked_seq_mxfold[i] + "\n")
  
  ofile.close()
  
