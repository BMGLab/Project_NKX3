#!/usr/bin/env bash
# 16_rootcause_09_construct_reads.sh
# Which reads did STAR leave unmapped, and are they plasmid (construct)-derived?
#  (1) unmapped read pairs = trimmed FASTQ pairs absent from STAR's genome BAM
#  (2) classify every unmapped pair by exact 25-mer matches to construct references:
#        FP     : GFP-family fluorescent-protein ORF (EGFP and EYFP-type variants;
#                 reads in this dataset carry the EYFP-type codon Tyr203)
#        bsd    : blasticidin-S deaminase ORF (blasticidin resistance)
#        neo    : neomycin phosphotransferase II ORF (G418 resistance of pcDNA3.1 vectors)
#        junctions observed in this dataset: CMV/T7 5'UTR-attB1-NKX3.1 ATG and
#                 NKX3.1 3'UTR(20 nt)-TOPO-attB2-FP
#      Reference ORFs are standard published sequences; their presence in the
#      reads is validated here (fraction of reference 25-mers seen in the reads).
#  (3) per gene: STAR-unmapped reads that reached the toTranscriptome BAM, by class
# usage: 16_rootcause_09_construct_reads.sh <sample_id> <run_tag>
set -euo pipefail
SAMPLE=$1; TAG=$2
WORK=/mnt/ssd1/Projects/NKX3_wd/work/rootcause
OUT=/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/rootcause
PY=/home/biolab/miniconda3/envs/nkx3_py/bin/python
SAMTOOLS=/usr/local/bin/samtools
SD=$WORK/star_${TAG}/$SAMPLE
CD=$WORK/construct_${TAG}/$SAMPLE
RD=$OUT/construct/${SAMPLE}_${TAG}
FQ1=$WORK/fastq_trim/$SAMPLE/${SAMPLE}_1_val_1.fq.gz
FQ2=$WORK/fastq_trim/$SAMPLE/${SAMPLE}_2_val_2.fq.gz
mkdir -p "$CD" "$RD"
export TMPDIR=$WORK/tmp; mkdir -p "$TMPDIR"

if [ ! -s "$CD/unmapped_R2.fq" ]; then
  echo "[$(date)] (1) unmapped read pairs"
  zcat "$FQ1" | awk 'NR%4==1{split($1,a," "); print substr(a[1],2)}' | sort -u -S 8G --parallel=8 -T "$TMPDIR" > "$CD/all_names.txt"
  comm -23 "$CD/all_names.txt" "$SD/genome_names.txt" > "$CD/unmapped_names.txt"
  rm -f "$CD/all_names.txt"
  for m in 1 2; do
    F=$FQ1; [ $m = 2 ] && F=$FQ2
    zcat "$F" | awk 'NR==FNR{n["@"$1]; next} FNR%4==1{split($1,a," "); p=(a[1] in n)} p' "$CD/unmapped_names.txt" - > "$CD/unmapped_R$m.fq"
  done
fi
wc -l < "$CD/unmapped_names.txt" | awk '{print "unmapped pairs:",$1}'
# read -> gene for STAR-unmapped reads present in the toTranscriptome BAM
if [ ! -s "$CD/txonly_read_gene.tsv" ]; then
  $SAMTOOLS view -@ 4 "$SD/tx_only_records.bam" | cut -f1,3 | awk -F'\t' -v OFS='\t' 'NR==FNR{g[$1]=$3; next}{print $1, g[$2]}' "$WORK/ref/tx2gene_name.tsv" - | sort -u -S 4G -T "$TMPDIR" > "$CD/txonly_read_gene.tsv"
fi

echo "[$(date)] (2-3) classify"
$PY -I - "$CD" "$RD" "$SAMPLE" <<'EOF'
import sys, collections
cd, rd, sample = sys.argv[1:4]
def rc(s): return s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]
EGFP = ("ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAG"
        "CTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGAC"
        "TTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATC"
        "GAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAAC"
        "TTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCC"
        "GCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA")
def mutate(s, codon, new):  # codon: 1-based codon index in the ORF
    i = 3 * (codon - 1); return s[:i] + new + s[i + 3:]
EYFP = EGFP
for c, n in [(66, 'GGC'), (69, 'CTG'), (73, 'GCC'), (204, 'TAC')]:   # T65G V68L S72A T203Y (GFP numbering)
    EYFP = mutate(EYFP, c, n)
BSD = ("ATGGCCAAGCCTTTGTCTCAAGAAGAATCCACCCTCATTGAAAGAGCAACGGCTACAATCAACAGCATCCCCATCTCTGAAGACTACAGCGTCGCCAGCGCAGCTCTCTCTAGCGACGGCCGC"
       "ATCTTCACTGGTGTCAATGTATATCATTTTACTGGGGGACCTTGTGCAGAACTCGTGGTGCTGGGCACTGCTGCTGCTGCGGCAGCTGGCAACCTGACTTGTATCGTCGCGATCGGAAATGAG"
       "AACAGGGGCATCTTGAGCCCCTGCGGACGGTGCCGACAGGTGCTTCTCGATCTGCATCCTGGGATCAAAGCCATAGTGAAGGACAGTGATGGACAGCCGACGGCAGTTGGGATTCGTGAATTG"
       "CTGCCCTCTGGTTATGTGTGGGAGGGCTAA")
NEO = ("ATGATTGAACAAGATGGATTGCACGCAGGTTCTCCGGCCGCTTGGGTGGAGAGGCTATTCGGCTATGACTGGGCACAACAGACAATCGGCTGCTCTGATGCCGCCGTGTTCCGGCTGTCAGCG"
       "CAGGGGCGCCCGGTTCTTTTTGTCAAGACCGACCTGTCCGGTGCCCTGAATGAACTGCAGGACGAGGCAGCGCGGCTATCGTGGCTGGCCACGACGGGCGTTCCTTGCGCAGCTGTGCTCGAC"
       "GTTGTCACTGAAGCGGGAAGGGACTGGCTGCTATTGGGCGAAGTGCCGGGGCAGGATCTCCTGTCATCTCACCTTGCTCCTGCCGAGAAAGTATCCATCATGGCTGATGCAATGCGGCGGCTG"
       "CATACGCTTGATCCGGCTACCTGCCCATTCGACCACCAAGCGAAACATCGCATCGAGCGAGCACGTACTCGGATGGAAGCCGGTCTTGTCGATCAGGATGATCTGGACGAAGAGCATCAGGGG"
       "CTCGCGCCAGCCGAACTGTTCGCCAGGCTCAAGGCGCGCATGCCCGACGGCGAGGATCTCGTCGTGACCCATGGCGATGCCTGCTTGCCGAATATCATGGTGGAAAATGGCCGCTTTTCTGGA"
       "TTCATCGACTGTGGCCGGCTGGGTGTGGCGGACCGCTATCAGGACATAGCGTTGGCTACCCGTGATATTGCTGAAGAGCTTGGCGGCGAATGGGCTGACCGCTTCCTCGTGCTTTACGGTATC"
       "GCCGCTCCCGATTCGCAGCGCATCGCCTTCTATCGCCTTCTTGACGAGTTCTTCTGA")
JUNC = ("GACTCACTATAGGGAGACCCAAGCTGGCTAGTTAAGCTGAGCATCAACAAGTTTGTACAAAAAAGCAGGCTCCGAATTCGCCCTTGGATCCATGCTCAGGGTTCCGGAGCCGC"  # 5'UTR-attB1-NKX3.1
        "NNNNN"
        "GGAGCCCAGCTTTTTGGTAATGCCAGCTCAGGTGACAACCGAATTCAAGGGCGAATTCGACCCAGCTTTCTTGTACAAAGTGGTTGATGCTGTTAACATGGTGAGCAAGGGCGAGGAG")  # NKX3.1-attB2-FP
refs = {'FP': [EGFP, EYFP], 'bsd': [BSD], 'neo': [NEO], 'NKX3.1_construct_junctions': [JUNC]}
with open(f"{rd}/../construct_references.fa", 'w') as o:   # written for downstream scripts / inspection
    for name, seq in [('FP_EGFP_ORF', EGFP), ('FP_EYFPtype_ORF(T65G,V68L,S72A,T203Y)', EYFP), ('bsd_ORF', BSD),
                      ('neo_nptII_ORF', NEO), ('NKX3.1_construct_junctions(5prime|NNNNN|3prime)', JUNC)]:
        o.write(f">{name}\n{seq}\n")
K = 25
kidx = {}
refk = {}
for name, seqs in refs.items():
    ks = set()
    for s in seqs:
        for i in range(len(s) - K + 1):
            k = s[i:i + K]
            if 'N' in k: continue
            ks.add(k); ks.add(rc(k))
    refk[name] = ks
    for k in ks: kidx.setdefault(k, name)
seen = collections.Counter()
cls = {}
n = 0
with open(f"{cd}/unmapped_R1.fq") as f1, open(f"{cd}/unmapped_R2.fq") as f2:
    while True:
        h = f1.readline()
        if not h: break
        s1 = f1.readline().strip(); f1.readline(); f1.readline()
        f2.readline(); s2 = f2.readline().strip(); f2.readline(); f2.readline()
        n += 1
        hits = collections.Counter()
        for s in (s1, s2):
            for i in range(len(s) - K + 1):
                k = s[i:i + K]; r = kidx.get(k)
                if r: hits[r] += 1; seen[k] += 1
        if hits:
            r, c = hits.most_common(1)[0]
            if c >= 3: cls[h[1:].split()[0]] = r
# reference validation: fraction of reference 25-mers (sense) observed >= 2 times
val = []
for name, seqs in refs.items():
    for j, s in enumerate(seqs):
        ks = [s[i:i + K] for i in range(len(s) - K + 1) if 'N' not in s[i:i + K]]
        obs = sum(1 for k in ks if seen[k] + seen[rc(k)] >= 2)
        val.append(dict(sample=sample, reference=name + (['', '_EYFPvariant'][j] if name == 'FP' else ''), n_25mers=len(ks),
                        frac_25mers_observed_ge2=round(obs / len(ks), 3) if ks else 'NA'))
c = collections.Counter(cls.values())
with open(f"{rd}/09_construct_reference_validation.tsv", 'w') as o:
    o.write('sample\treference\tn_25mers\tfrac_25mers_observed_ge2\n')
    for v in val: o.write('\t'.join(str(v[x]) for x in ['sample', 'reference', 'n_25mers', 'frac_25mers_observed_ge2']) + '\n')
with open(f"{rd}/09_unmapped_pairs_by_construct_class.tsv", 'w') as o:
    o.write('sample\tclass\tread_pairs\tfrac_of_STAR_unmapped_pairs\n')
    for k in list(refs) + ['none']:
        v = c.get(k, 0) if k != 'none' else n - len(cls)
        o.write(f"{sample}\t{k}\t{v}\t{v / n:.4f}\n")
    o.write(f"{sample}\tTOTAL_STAR_unmapped_pairs\t{n}\t1\n")
# per gene: STAR-unmapped reads that reached the toTranscriptome BAM
pg = collections.defaultdict(collections.Counter)
for l in open(f"{cd}/txonly_read_gene.tsv"):
    r, g = l.rstrip('\n').split('\t')
    pg[g]['txonly_reads'] += 1
    pg[g][cls.get(r, 'none')] += 1
rows = sorted(pg.items(), key=lambda x: -x[1]['txonly_reads'])
with open(f"{rd}/09_per_gene_STARunmapped_reads_by_class.tsv", 'w') as o:
    cols = ['txonly_reads'] + list(refs) + ['none']
    o.write('gene_name\t' + '\t'.join(cols) + '\tfrac_construct\n')
    for g, cc in rows:
        fc = (cc['txonly_reads'] - cc['none']) / cc['txonly_reads']
        o.write(g + '\t' + '\t'.join(str(cc[x]) for x in cols) + f"\t{fc:.3f}\n")
print(open(f"{rd}/09_construct_reference_validation.tsv").read())
print(open(f"{rd}/09_unmapped_pairs_by_construct_class.tsv").read())
print(''.join(open(f"{rd}/09_per_gene_STARunmapped_reads_by_class.tsv").readlines()[:30]))
EOF
echo "[$(date)] done"
