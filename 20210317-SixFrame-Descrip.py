from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys

if(len(sys.argv))<2:
    sys.exit("You must provide the name of a FASTA file as a command-line argument for this script")
TargetName = sys.argv[1]
MinLength = 100
for Transcript in SeqIO.parse(TargetName, "fasta"):
    for strand, nuc in [(+1, Transcript.seq), (-1, Transcript.seq.reverse_complement())]:
        for frame in range(3):
            length = 3*((len(Transcript)-frame) // 3)
            counter = 65
            for pro in nuc[frame:frame+length].translate().split("*"):
                if len(pro) >= MinLength:
                    accession = "%s:str%i:frm%i" % (Transcript.id, strand, frame)
                    descrip = Transcript.description.replace(Transcript.id,'')
                    sequence = pro
                    print(">%s%s%s" % (accession, chr(counter), descrip))
                    print(sequence)
                    counter = counter+1
