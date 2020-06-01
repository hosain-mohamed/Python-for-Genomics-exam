""" 
    - This was a final exam in (Python for genomics at Coursera which is a course 
    part of Genome Data Science Specializrion : for asnwering multiple questions about 
    FASTA file containg number of sequences like : 
        - How many recoreds ? 
        - what is the longest recored ? 
        - what is the shortest recored ? 
        - what is the longest ORF and it's starting position ? 
        - what is most frequent repeats in the FASTA file
    - there is tow files attached to work with ("dna.example.fasta","dna2.fasta")    
 """

filename = "dna2.fasta"

def arrange_seq():
    with open(filename,"r") as f:
        seq_dict ={}
        seq_lines = f.readlines()
        for line in seq_lines:
            if line.strip().startswith(">"):
                key = line.split()[0][1:]
                seq_dict[key]=""
            else:
                seq_dict[key]+=line.strip()
        return seq_dict

seq_dict = arrange_seq() 

print("- Total records are " , len(seq_dict.keys()))

def seq_lengthes():
    longest_length = -1
    longest_indentifiers = []
    shortest_length = float("inf")
    shortest_identifiers = []
    for value in seq_dict.values():
        longest_length = max(longest_length,len(value))
        shortest_length = min(shortest_length,len(value))
    for key,value in seq_dict.items():
        if len(value) == longest_length:
            longest_indentifiers.append(key)
        if len(value) == shortest_length:
            shortest_identifiers.append(key)
    print("- the longest identifiers  are \n % s \nthere's length is % d "
     % (longest_indentifiers , longest_length))
    
    print("- the shorts identifiers  are \n % s \nthere length is % d " 
    % (shortest_identifiers , shortest_length))

seq_lengthes()
def longest_ORF_with_start(t,seq):
    longest_ORF = -1
    longest_ORF_start_pos = -1
    value = seq[t-1:]
    seq_list =[]
    seq = seq[t-1:]
    for i in range(0,len(seq) -(len(seq)%3),3):
        base = seq[i:i+3]
        seq_list.append(base)
    i = t
    while i < len(value) - 2:
        if value[i:i+3] =="ATG":
            start = i
            i = i+3
            while i < len(value) -2:
                if(value[i:i+3]) in ("TAA","TAG","TGA"):
                    if(longest_ORF<i+3 - start):
                        longest_ORF = i+3 - start
                        longest_ORF_start_pos = start
                    break
                i+=1
        i+=1    
    return [longest_ORF,longest_ORF_start_pos]
def longest_ORF_began_with(index,seq):
    longest_ORF = -1
    longest_ORF_start_pos = -1
    seq = seq[index-1:]
    seq_list =[]
    for i in range(0,len(seq) -(len(seq)%3),3):
        base = seq[i:i+3]
        seq_list.append(base)
    i = 0 
    start_codon = False 
    while i < len(seq_list):
        if seq_list[i] in ("ATG","atg") and start_codon == False:
            start = i *3
            start_codon = True
            i+=1
        if seq_list[i] in ('TAA','TAG','TGA','taa','tag','tga') and start_codon==True:
            if (i*3 +3 -start) > longest_ORF :
                longest_ORF = i*3 +3 -start
                longest_ORF_start_pos = start + index
            start_codon = False
            i+=1
        else:
            i+=1 
    return[longest_ORF,longest_ORF_start_pos]

def longest_ORF_with_identifier(index,key):
    value = seq_dict[key]
    lst = longest_ORF_began_with(index,value)
    print("- in this key longest ORF is %s and it's starting positon is %s" 
    %(lst[0],lst[1]))

longest_ORF_with_identifier(3,"gi|142022655|gb|EQ086233.1|16")

def get_seq_by_identifier(key):
    value = seq_dict[key]
    return value

def longest_ORF_in_file(index): 
    longest_ORF = -1 
    longest_ORF_start_pos = -1
    longest_identifier = ''
    for key , value in seq_dict.items():
        lst = longest_ORF_began_with(index,value)
        if lst[0] > longest_ORF:
            longest_ORF = lst[0]
            longest_ORF_start_pos = lst[1]
            longest_identifier = key
    
    print("- longest ORF is %s \n it's start position is %s \n it's  identifier is %s"
    %(longest_ORF,longest_ORF_start_pos,longest_identifier))
     
longest_ORF_in_file(3)

def arrange_repeats_with_length_n(seq,n):
    most_frequent_count = -1
    repeats = {}
    for i in range(len(seq)-2):
        base = seq[i:i+n]
        repeats[base] = repeats.get(base,0)
        repeats[base] +=1
        if repeats[base] > most_frequent_count:
            most_frequent_count = repeats[base]
    return [most_frequent_count,repeats]
st="GTCGATCGACACGACGCTCGCGCAGCGCGACGCGAAGGCCGCGTGAGCGCACGACGCGCGTCACACACCACGAGCACAACGAACACGACCCCACTCTCACGGAGCCGACCATGGCCGACCTTCGCTGCACCATCGCGGGCATCACTTCGCCGAACCCTTTCTGGCTGGCGTCCGCGCCGCCGACCGACAAGGCCTACAACGTGAACCGCGCGTTCGAGGCGGGCTGGGGCGGGGTCGTCTGGAAGACGCTCGGGCTCGATCCGCATGTCGTCAACGTCAGTTCGCGCTATGGCGCGGTGCAGTGGAACGGCCAGCGCATCGCGGGGCTGAACAACATCGAGCTGATCACCGACCGTCCGCTCGACGTGAACCTGAGAGAGATCGCGCAGGTGAAGCGCGACTGGCCGGACCGCGCGCTGATCGTGTCGCTGATGGTGCCGTGCAACGAGCGCGACTGGAAATGGATCCTGCCGCTCGTCGAGGATACGGGCGCCGACGCGGTCGAGCTGAACTTCGGTTGTCCGCACGGGATGAGCGAGCGCGGGATGGGCGCGGCGGTCGGGCAGGTGCCCGAATATGTGGAGATGGTCACGCGCTGGGTGAAGGAAGGCACGAAGCTGCCGTGCCTCGTGAAGCTCACGCCGAACATCAGCGACATCCGGATGGGGTCGCGCGCCGCGTACAAGGGCGGCGCGGACGGCGTGTCGCTGATCAACACGATCAACTCGATCGTCGCGGTCGATCTCGACCATATGGCGCCGATGCCGACGGTCGACGGCAAGGGCACGCACGGCGGCTATTGCGGCCCGGCGGTCAAGCCGATCGCATTGAACATGGTCGCGGAGATCGCACGTGACCCGGAAACGCCGAACCTGCCGATCTCGGGCATCGGCGGCATCTCGTCATGGCGCGACGCGGCGGAGTTCATGGTGCTCGGCGCCGGCAGCGTGCAGGTGTGCACCGCCGCGATGCATTACGGATTCCGGATCGTGTCGGACCTGGCCGACGGATTGTCGAACTGGATGGACGAGAAGGGCTACGCGACGCTCGACGACATTCGCGGCCGCGCGGTGCCGAACGTGACCGACTGGAAATACCTGAACCTGAAATACGACATCAAGGCGCGTATCGACCAGGACCGCTGCATCCAGTGCGGGTTGTGCCATATCGCGTGCGAGGACACGTCGCACCAGGCGATCACCGCGACGAAGGACGGCGTGCGGCATTTCGAAGTGGTCGATTCGGCGTGCGTCGGGTGCAATCTTTGCATGCATGTGTGTCCGGTCGAGCAATGCATCACGATGGAGCGTGTCGATTCGGGCGACTACGCGAACTGGACCACGCATCCGAACAATCCGGCGAGCGCGGAGGCGGGGGCGAGTGCAGGCGCGGCGGCACCCGAGAAGCACGCGAAGAAGGCTGCTTGACGGCGTCCGGCGATGCGGGCCATCCTGCATCGCCGCCTTTCGTTCCACCCGGGCCGGCATCGAGTGATGCCGGCGTTGACGTTTTCGTGGAGTGAGTCAGATGAATCACGCAGCGAATCCCGCCGATCCCGATCGCGCCGCGGCGCAGGGCGGCAGCCTGTACAACGACGATCTCGCGCCGACGACGCCGGCGCAGCGCACGTGGAAGTGGTATCACTTCGCGGCGCTGTGGGTCGGGATGGTGATGAACATCGCGTCGTACATGCTCGCGGCCGGGCTGATCCAGGAAGGCATGTCGCCGTGGCAGGCGGTGACGACGGTGCTGCTCGGCAACCTGATCGTGCTCGTGCCGATGCTGCTGATCGGCCATGCGGGCGCGAAGCACGGGATTCCGTACGCGGTGCTCGTGCGCGCGTCGTTCGGCACGCAGGGGGCGAAGCTGCCGGCGCTGCTGCGCGCGATCGTCGCGTGCGGCTGGTACGGGATCCAGACCTGGCTCGGCGGCAGCGCGATCTATACGCTGCTGAACATCCTGACCGGCAACGCGCTGCATGGCGCCGCGCTGCCGGTCATCGGCATCGGGTTCGGGCAGCTCGCATGCTTCCTCGTGTTCTGGGCGCTGCAGCTCTACTTCATCTGGCATGGCACCGATTCGATCCGCTGGCTCGAAAGCTGGTCGGCGCCGATCAAGGTCGTGATGTGCGTGGCGCTGGTGTGGTGGGCAACGTCGAAGGCGGGCGGCTTCGGCACGATGCTGTCGGCGCCGTCGCAGTTTGCCGCAGGCGGCAAGAAAGCCGGGCTGTTCTGGGCGACCTTCTGGCCGGGGCTGACCGCGATGGTCGGCTTCTGGGCGACGCTCGCGCTGAACATCCCCGACTTCACGCGCTTCGCGCATTCGCAGCGCGACCAGGTGATCGGCCAGTCGATCGGGCTGCCGTTGCCGATGGCGCTGCTGTCGGTGGTGTCGGTCGTCGTGACGTCGGCGACCGTCGTGATCTACGGCAACGCGATCTGGGATCCGATCGACCTGACGAGCCGGATGACGGGCATCGGCGTGGGCATCGCGCTCGTGATCCTCACGCTCGACACGATGTGCTGCAACCTCGCCGCGAATCTCGTCGGCCCGGCGTACGACTTCTCGAGCCTGTGGCCGAAGGCGATCTCGTACCGCACCGGCGGGATGATCACCGCGACGCTCGCGATCGTGATGATGCCGTGGAAGATCCTCGCGACGACGGACGGCTACATCTTCACCTGGCTCGTCGGCTACTCGGCGCTGCTCGGGCCCGTGGCGGGGATCCTGATGGTCGACTACTTCCTGATTCGCGGCACGCGGCTCGACACGCGCGCGCTGTTCGACGAGCGCGGCGGCTTCAGCTACGCGCGCGGCTGGAACCCGGCCGCGCTGGCCGCGCTCGCGGTCGGCGTGCTGCCGAACCTGCCCGGCTTCCTGCACACGGCGTTTCCGGCGTCGTTTCCGAACGTGCCGGCGTTCTTCAACACGCTTTACACGTACGCGTGGTTCGTCGGCCTCGTGCTGGCGTCATGCGTGTACGGCACCTGGATGAAGTGGCGCGCCGGACAGCACGCGCAGATCGCGAGCGCCTGATTCGGCACCCGACAGTCAACGAGGAGGCAACCCCATGGCAATCCTGATTCGTGGCGGCACCGTGGTCGATGCGGACCGTTCCTACCGCGCGGACGTGCTCTGCGCAGCCCCGGAGGACGGCGGCACGATCCTGCAGATCGCCGGGCAGATCGATGCGCCGGCCGGCGCGACCGTCGTCGATGCGCACGACCAGTACGTGATGCCGGGCGGCATCGATCCGCATACGCACATGGAACTGCCGTTCATGGGCACGACCGCGAGCGACGATTTCTACTCGGGTACGGCCGCCGGGCTCGCGGGCGGCACGACGAGCATCATCGACTTCGTGATCCCGAGCCCGAAGCAGCCGCTGATGGACGCGTTCCATGCCTGGCGCGGCTGGGCCGAGAAGGCGGCGGCCGACTACGGCTTCCACGTGGCCGTGACGTGGTGGGACGAGAGTGTGCACCGCGACATGGGCACGCTCGTGCGCGAACACGGCGTGTCGAGCTTCAAGCACTTCATGGCGTACAAGAACGCGATCATGGCCGACGACGAGGTGCTCGTGAACAGCTTCTCGCGTTCGCTCGAACTCGGCGCGTTGCCGACCGTGCATGCGGAGAACGGCGAGCTCGTGTTCCAGTTGCAGAAGGCGCTGCTCGCGCGCGGGATGACGGGGCCGGAGGCGCATCCGCTGTCGCGGCCGCCGGAGGTCGAGGGTGAGGCGGCGAATCGTGCGATCCGCATTGCGCAGGTGCTCGGCGTGCCGGTGTATATCGTGCATGTGTCCGCGAAGGACGCGGTCGATGCGATCACGAAGGCGCGCAGCGAAGGGCTGCGCGTGTTCGGCGAGGTGCTGCCGGGCCATCTGGTGATCGACGAGGCCGTCTATCGCGATCCGGACTGGACACGTGCGGCCGCGCACGTGATGAGCCCGCCGTTCCGCTCGGCCGAGCACCGCGAGGCGCTGTGGCGCGGGCTGCAGGCAGGGCAGCTGCATACGACGGCAACCGACCACTGCGTGTTCTGCGCGTCGCAGAAGGCGATGGGCCGCGAGGATTTCACGAAGATCCCGAACGGCTGCGGCGGTGTCGAGGATCGCATGTCGGTGCTGTGGCATCACGGCGTGAATCATGGCCGCATCACGCCGAACGAGTTCGTGCGGATCACGTCGACGAACGCCGCGCAGATCTTCAACCTGTATCCGCGCAAGGGCGCCGTGCAGGTGGGCGCCGATGCCGACCTCGTCGTGTGGGACCCGGCCGCGACCAGGACGATCTCGGTGAAGACGCATCACCAGCAGGTCGATTTCAACGTGTTCGAGGGGATGACCGTACAAGGCGTCGCAACCCACACGCTCACGCGCGGCGCGCTCGCGTGGGCCGACGGCGATCTGCGTGCCGTGCGCGGCGCGGGCCGCTATCTGAAGCGCCCGCCGGCAGCCAGCTACTACGAGGCCGCGCGGATCGCGAACCGGCTGCGCGAACCGCATCCGGTCGAGCGCGCCGGTTGAGCGTTGCGTATCGCGCGGGGCGTGTCGGTTCGAACGACACGCCCCGCGCATGTTTGAGCGTGCGTTTACGTTCGTGCCGGCACCGCTCGTGCCGCCTCTCCGCATCACGCCCATCCTCTCAATATTTGGGATGAATTGAGCGCGATCGCGCCTTGCCGATCTCCGGATACATAGAACAACTGAGCAAGTCGATGAAACACGCGATGTCGCGCAAATGCGACCATTTTGTTTGCGTTGTCGACGTGCATGCCGGCGAGTAATATCCACCGACGGCGT"
arrange_repeats_with_length_n(st,6)

def most_freq_repeats_in_file_with_length_n(n):
    total_seq = ''
    for value in seq_dict.values():
        total_seq = total_seq + value
    lst = arrange_repeats_with_length_n(total_seq,n)
    most_frequent_count = lst[0]
    repeats = lst[1]
    most_frequent_sequences = []
    for key,value in repeats.items():
        if value == most_frequent_count:
            most_frequent_sequences.append(key)
    return [lst[0],most_frequent_sequences]
 

print("- the length frequent repeat is %s of length 6 " 
    %( most_freq_repeats_in_file_with_length_n(6)[0])) 

print("- the number of frequent repeats of length 12 is %s " 
    %( most_freq_repeats_in_file_with_length_n(7)[1])) 

