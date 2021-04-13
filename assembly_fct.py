def seq_comp_right(contig,seq_test):
    """return the best alignement between two sequence to the right side
    
    if the two sequences show a similarity of 90% the algorithm will
    considere that they are "the same" and the mismatch will be corrected
    based on the tested sequence.
    this function call the functions _error_comp() and _mismatch_repair()   
    """
    for i in range(len(seq_test),0,-1):

        if seq_test[0:i] == contig[len(contig)-i:len(contig)]:
            return i, contig

        elif  _error_comp(seq_test[len(seq_test)-i:len(seq_test)],
        list(contig[0:i]),i,"left") >= 90:
            
            contig = _mismatch_repair(contig,
                    seq_test[len(seq_test)-i:len(seq_test)],
                    list(contig[0:i]), i, "left")

            return i, contig
    #return 0 if there is no alignement
    return 0, contig

def seq_comp_left(contig,seq_test):
    """return the best alignement between two sequence to the left side
    
    if the two sequences show a similarity of 90% the algorithm will
    considere that they are "the same" and the mismatch will be corrected
    based on the tested sequence.
    this function call the functions _error_comp() and _mismatch_repair()   
    """
    for i in range(len(seq_test),0,-1):

        if seq_test[len(seq_test)-i:len(seq_test)] == contig[0:i]:
            return i, contig

        elif  _error_comp(seq_test[len(seq_test)-i:len(seq_test)],
        list(contig[0:i]),i,"left") >= 90:
            
            contig = _mismatch_repair(contig,
                    seq_test[len(seq_test)-i:len(seq_test)],
                    list(contig[0:i]), i, "left")

            return i, contig
        
    #return 0 if there is no alignement
    return 0, contig

def _error_comp(seq_test,seq_ref_list,i,orientation):
    """function returning the percent of similarity between two sequence
    
    the similarity is considered as the percentage of similare nucleotide
    between the two sequence
    """
    correct = 0
    for j in range(0,len(seq_ref_list)):
        if seq_ref_list[j] == seq_test[j]:
            correct += 1

    pourcent = (correct * 100)/len(seq_ref_list)

    return pourcent
    
def _mismatch_repair(contig, seq_test,seq_ref_list,i,orientation):
    """correct the mismatch in the contig

    if the similarity between two sequence is more than 90%
    then this function will change the nucleotide in the contig
    with the nucleotide at the same position from the aligned sequence.
    it return the corrected contig
    """
    seq_ref_correct = ""
    for j in range(0,len(seq_ref_list)):
        if seq_ref_list[j] != seq_test[j]:
            seq_ref_list[j] = seq_test[j]

    seq_ref_correct = "".join(seq_ref_list)
    if orientation == "right":
        contig = contig[0:len(contig)-i] + seq_ref_correct
    elif orientation == "left":
        contig = seq_ref_correct + contig[i:]
    
    return contig


