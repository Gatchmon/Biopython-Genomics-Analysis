from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from matplotlib import pyplot as plt
from reportlab.lib.colors import red, grey, orange, green, brown, blue, lightblue, purple

from Bio import SeqIO
from Bio.SeqUtils import GC
import pylab

record = SeqIO.read("NC_003436.gbk","gb")
table = 11
min_pro_len = 100
def find_orfs_with_trans(seq, trans_table, min_protein_length):
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end-aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame+aa_start*3
                        end = min(seq_len,frame+aa_end*3+3)
                    else:
                        start = seq_len-frame-aa_end*3-3
                        end = seq_len-frame-aa_start*3
                    answer.append((start, end, strand,
                                   trans[aa_start:aa_end]))
                aa_start = aa_end+1
    answer.sort()
    return answer
orf_list = find_orfs_with_trans(record.seq, table, min_pro_len)
for start, end, strand, pro in orf_list:
    print("%s...%s - length %i, strand %i, %i:%i" \
          % (pro[:30], pro[-3:], len(pro), strand, start, end))


A_rec = SeqIO.read("NC_022103.gbk", "gb")
B_rec = SeqIO.read("NC_002645.gbk", "gb")
C_rec = SeqIO.read("NC_030292.gbk", "gb")
D_rec = SeqIO.read("NC_003436.gbk", "gb")

A_colors = [red]*7
B_colors = [green]*7
C_colors = [blue]*9
D_colors = [orange]*7

name = "Genome Diagram of GB Records"
gd_diagram = GenomeDiagram.Diagram(name)
max_len = 0
for record, gene_colors in zip([A_rec, B_rec, C_rec, D_rec], [A_colors, B_colors, C_colors, D_colors]):
    max_len = max(max_len, len(record))
    gd_track_for_features = gd_diagram.new_track(1,
                            name=record.name,
                            greytrack=True,
                            start=0, end=len(record))
    gd_feature_set = gd_track_for_features.new_set()

    i = 0
    for feature in record.features:
        if feature.type != "gene":
            #Exclude this feature
            continue
        gd_feature_set.add_feature(feature, sigil="ARROW",
                                   color=gene_colors[i], label=True,
                                   name= str(i+1),
                                   label_position="start",
                                   label_size = 6, label_angle=0)
        i+=1

gd_diagram.draw(format="linear", pagesize='A4', fragments=1,
                start=0, end=max_len)
gd_diagram.write(name + ".pdf", "PDF")

#------------------------ CREATE GENOME DIAGRAM FROM GENBANK FILE -------------------------------------

record1 = SeqIO.read('NC_002645.gbk', 'genbank')

gd_diagram = GenomeDiagram.Diagram('Human Coronavirus 229E')
gd_track_for_features = gd_diagram.new_track(1, name='Annotated Features')
gd_feature_set = gd_track_for_features.new_set()

for feature in record.features:
    if feature.type != 'gene':
        continue
    if len(gd_feature_set)%2 == 0:
        color = colors.blue
    else:
        color = colors.lightblue
    gd_feature_set.add_feature(feature, color=color, label=True)

gd_diagram.draw(format='linear', orientation='landscape', pagesize='A4'
               , start=0, end=len(record))
gd_diagram.write('virus_linear.png', 'png')
gd_diagram.write('virus_linear.pdf', 'pdf')
gd_diagram.draw(format='circular', circular=True, pagesize=(20*cm,20*cm),
                start=0, end=len(record), circle_core=0.7)
gd_diagram.write('virus_circular.pdf', 'pdf')

gc_content = [40.82, 38.51, 38.26, 40.24, 39.00, 41.78, 40.99, 41.81, 42.01, 39.28, 34.46]
acc_num = ['NC_022103','NC_018871.1','NC_002645','NC_032730','NC_030292','NC_010438','NC_028811.1','NC_028833','NC_003436','NC_009988.1','NC_005831']
label = ['Colacovirus', 'Decacovirus', 'Duvinacovirus', 'Luchacovirus', 'Minacovirus', 'Minunacovirus', 'Myotacovirus',
            'Nyctacovirus', 'Pedacovirus', 'Rhinacovirus', 'Setracovirus']

plt.scatter(acc_num, gc_content, color='g')
plt.title('GC Content - Subgenus')
plt.xlabel('Subgenus')
plt.ylabel('GC Content')
plt.show()







