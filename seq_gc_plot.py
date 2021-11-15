from Bio import SeqIO
import pylab
from Bio.SeqUtils import GC

sizes = [len(rec) for rec in SeqIO.parse("sequence.gb", "genbank")]
#-----------------------------------SEQUENCE LENGTH PLOT------------------------------------------------

pylab.hist(sizes, bins=20,color='lightblue' ,histtype='barstacked', label='str', edgecolor='black', linewidth=1.0)
pylab.title("%i Alphacoronavirus Complete Genome\nLengths %i to %i" \
            % (len(sizes),min(sizes),max(sizes)))
pylab.xlabel("Sequence length (bp)")
pylab.ylabel("Count")
pylab.show()

#-----------------------------------GC CONTENT PLOT------------------------------------------------

gc_values = sorted(GC(rec.seq) for rec in SeqIO.parse("sequence.gb", "genbank"))
pylab.plot(gc_values)
pylab.title("%i Alphacoronavirus Complete Genome\nGC%% %0.1f to %0.1f" \
            % (len(gc_values),min(gc_values),max(gc_values)))
pylab.xlabel("Genes")
pylab.ylabel("GC%")
pylab.show()

#-----------------------------------NUCLEOTIDE PLOT------------------------------------------------
#UNABLE TO SHOWS NUCLEOTIDE PLOT, WILL FIX THIS CODE LATER ON

# with open("human_bat_alphacoronavirus.gbk") as in_handle:
#     record_iterator = SeqIO.parse(in_handle, "genbank")
#     rec_one = next(record_iterator)
#     rec_two = next(record_iterator)
#
# window = 7
# seq_one = str(rec_one.seq).upper()
# seq_two = str(rec_two.seq).upper()
# data = [[(seq_one[i:i + window] != seq_two[j:j + window])
#          for j in range(len(seq_one) - window)]
#          for i in range(len(seq_two) - window)]
#
# pylab.gray()
# pylab.imshow(data)
# pylab.xlabel("%s (length %i bp)" % (rec_one.id, len(rec_one)))
# pylab.ylabel("%s (length %i bp)" % (rec_two.id, len(rec_two)))
# pylab.title("Dot plot using window size %i\n(allowing no mis-matches)" % window)
# pylab.show()