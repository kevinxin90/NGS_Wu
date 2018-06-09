from Bio.Blast.Applications import NcbiblastnCommandline
from io import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import defaultdict, Counter
import csv


def parse_fastq(file_name):
    """
    This function aims at parsing multiple fastq sequences contained in one file
    Steps:
    1) Read in file
    2) Each 4 lines is a group, with line1 been FASTA title line, 
                                    line 2 been raw sequence letters
                                    line 3 been "+" character
                                    line 4 been quality values
    3) make sure line1, line2,line3,line4 is correct format
    4) extract line 2 only, 
    
    Params
    ======
    file_name: 
        the file path that stores fastq info
    
    Return
    ======
    list of FASTA sequences as str
    """
    with open(file_name) as f:
        data = f.readlines()
        # count line number, initialize with 0
        line_num = 0
        # final_output
        results = []
        for i in range(0, len(data), 4):
            if data[i].startswith("@") and data[i+2].startswith("+"):
                results.append(data[i+1].strip('\n'))
            else:
                print('The {} line does not comply with the fasta format'.format(data[i+1]))
    return results

def get_seq_pos(fasta_list):
    seq_pos_list = []
    for i, _seq in enumerate(fasta_list):
        if i%500 ==0:
            print('{} has been processed!'.format(i))
        seq1 = SeqRecord(Seq(_seq))
        SeqIO.write(seq1, "seq1.fasta", "fasta")
        # Run BLAST and parse the output as XML
        try:
            output = NcbiblastnCommandline(query="seq1.fasta", subject="parent.fasta", outfmt=5)()[0]
            blast_result_record = NCBIXML.read(StringIO(output))
            # Print some information on the result
            if blast_result_record.alignments != []:
                hsps = blast_result_record.alignments[0].hsps
                if len(hsps) == 2:
                    results = []
                    for hsp in hsps:
                        results.append(hsp.sbjct_start)
                        results.append(hsp.sbjct_end)
                    seq_pos_list.append(sorted(results))  
        except:
            print('failed to blast!')
            continue  
    return seq_pos_list

def count_hsp_num(fasta_list):
    """
    count number of alignments
    and group them based on the alignment number
    """
    k = defaultdict(list)
    for i, _seq in enumerate(fasta_list):
        if i%500 ==0:
            print('{} has been processed!'.format(i))
        seq1 = SeqRecord(Seq(_seq))
        SeqIO.write(seq1, "seq1.fasta", "fasta")
        # Run BLAST and parse the output as XML
        output = NcbiblastnCommandline(query="seq1.fasta", subject="parent.fasta", outfmt=5)()[0]
        blast_result_record = NCBIXML.read(StringIO(output))
        # Print some information on the result
        if blast_result_record.alignments != []:
            hsps = blast_result_record.alignments[0].hsps

            k[len(hsps)].append(_seq)
        else:
            k[-1].append(_seq)
    return k

def count_align_length(input_sequence):
    seq1 = SeqRecord(Seq(input_sequence))
    SeqIO.write(seq1, "seq1.fasta", "fasta")
    # Run BLAST and parse the output as XML
    output = NcbiblastnCommandline(query="seq1.fasta", subject="parent.fasta", outfmt=5)()[0]
    blast_result_record = NCBIXML.read(StringIO(output))
    results = []
    for hsp in blast_result_record.alignments[0].hsps:
        results.append(hsp.align_length)
        results.append(hsp.align_length)
    return results

def seq_pos(input_sequence):
    seq1 = SeqRecord(Seq(input_sequence))
    SeqIO.write(seq1, "seq1.fasta", "fasta")
    # Run BLAST and parse the output as XML
    output = NcbiblastnCommandline(query="seq1.fasta", subject="parent.fasta", outfmt=5)()[0]
    blast_result_record = NCBIXML.read(StringIO(output))
    # Print some information on the result
    results = []
    for hsp in blast_result_record.alignments[0].hsps:
        results.append(hsp.sbjct_start)
        results.append(hsp.sbjct_end)
    return sorted(results)

def extractMiddle(seq_pos_list):
    results = []
    for _item in seq_pos_list:
        results.append(str(_item[1]) + '-' + str(_item[2]))
    return results

def main_func(input_file_path, output_file_path):
    # read in sequence for comparison
    PARENTAL_SEQUENCE = SeqRecord(Seq("AGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGTAAGTCTCGACGAAACAAGGATGCTGTTAGAGTTTCACAAACCTAAATCCGGTTGCTTTAGTTGTCTTAAAAGCTTTGAGCAAAAACCTTGATTGTTCTCCAGTGGAGGTGCAGTCACTGCCCTCTATCCGTTGGTCATTTCATTTGTGTCCTTGCCTGTTGGCAAGACTCCACTGAAACCTCTCTGGGAGATTGGTAGGTGGAGGGGGCAGGAGGCCCTACTTAGAAAGTGTCATTGAAGCCAATCCTTCTAACTGACCACCTCTGCCCTCCTAATAATTCTGGTGTGAAGGCGTAATGATGTGGGCTTCAGGGTCTTTGTTCTTCCTCCCCTAAGTCTTCAGAATGGGTAGTTGGGAGTAAGGGTGGTAGAAGGGGAACTGGATGAAGTGGACATGGTGGGGGTCTTCCCATAGAGGGTCCCTCATTGACTAGAGCTCAGATCTATTACCCTGTTATCCCTAGAGCTCCTTTTTCCTTTTCAGATTGCACCACCCGAAACACCTGACTCCAAAGTTCGTATGGTTCTCGAGAAGAAGAAGGGCCCAGCATCGCTGGACCGCCTGAGGCCCAATTCAAGGTGAGGGTCTTTATTGTTTTCCAAGTCGAGTAGGGATAACAGGGTAATCTCGACTGAACATGGAGGAATTGAGGTTGGGTATTTCCCCTGAGGTAGGAAAAAGGCTGGGTCAGTTTCCCGTTAGCCGTCAAGTCCTCATCACATCTTTAAGCCTTCCATGCAGGATAAAGGGCTGCAGAGCTATTTTCAAATTGACATCAAACTGGATTTCTGTTGACTTCGTCTTCCCTTTTTAAGGTCCACAGAAGAAGATGGGAAGGAAAGAAGTCTGAGGGCATCTTATTTGCACTCCGCTGTCATTTCTAAGGAAGGCCTTTAATGCCAAATTCTCATCTTTTATGTCCCCACTAAATCCTAAGGTTCTTGAACTTCTGATCAGACAGCCAAAAAATGAACCATCAACTAGCTTAACCTAACATATGTGAGGATAGAGGACTGGGACAGCTCTCTGGGCCACTGGAGAGTCAGACAGGCCTGCCCTCTGTGTGATTGACCGCGGTCTCTTTCTTCCAGGAGCGCACCATCTTCTTCAAGGACGACGG"), id="seq2")
    SeqIO.write(PARENTAL_SEQUENCE, "parent.fasta", "fasta")

    # read file
    fasta_list = parse_fastq(input_file_path)


    seq_pos_list = get_seq_pos(fasta_list)

    middle_pos = extractMiddle(seq_pos_list)

    count_middle_pos = Counter()

    for _pos_type in middle_pos:
        count_middle_pos[_pos_type] += 1

    with open(output_file_path,'w') as csvfile:
        fieldnames=['interval', 'count']
        writer=csv.writer(csvfile)
        writer.writerow(fieldnames)
        for key, value in count_middle_pos.items():
            writer.writerow([key, str(value)]) 




