"""Assembly algorithm made for a course exercice
"""

__authors__ = "Jules Garreau"
__contact__ = "jules.garreau00@gmail.com"
__version__ = "1.0.0"
__copyright__ = "copyleft"
__date__  = "13/04/21"

import assembly_fct as fct

#make a list of all the sequence from the txt file
path = __file__.split("assembly.py")
path = path[0]
path = path + "sequence3.txt"

with open(path,"r") as sequence:
    sequence_string = sequence.read()
seq_list = sequence_string.split()

#initialize the direction to concatenate r = right, l = left
concatenate = "r"

#first contig is the first sequence of the list
contig = seq_list[0]
del seq_list[0]

#algorithm will turn untill all the sequence are concatenated
while seq_list:

    best_overlap_index = 0
    best_seq = 0
    
    #compute the best alignement to the right or the the left
    for i in range(len(seq_list)):
        right_overlap_index, contig = fct.seq_comp_right(contig,seq_list[i])
        left_overlap_index, contig = fct.seq_comp_left(contig,seq_list[i])
        
        if right_overlap_index > best_overlap_index and right_overlap_index > left_overlap_index:
            
            best_overlap_index = right_overlap_index
            best_seq = i
            concatenate = "r"
            
        if left_overlap_index > best_overlap_index and right_overlap_index < left_overlap_index:
            
            best_overlap_index = left_overlap_index
            best_seq = i   
            concatenate = "l"
    
    if concatenate == "r":
        contig = contig + seq_list[best_seq][best_overlap_index:]
    
    if concatenate == "l":
        contig = seq_list[best_seq] + contig[best_overlap_index:]

    del seq_list[best_seq]


#print the sequence assembled
print(contig)




