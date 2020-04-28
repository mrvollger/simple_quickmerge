import os
import json
import random
import tempfile
import pysam

SMKDIR = os.path.dirname(workflow.snakefile)
shell.executable("/bin/bash")
shell.prefix("source %s/env.cfg; set -eo pipefail; " % SMKDIR)

SSD_TMP_DIR = "/data/scratch/ssd"
if "TMPDIR" in os.environ:
    TMPDIR = os.environ['TMPDIR']
elif "TMPDIR" in config:
    TMPDIR = config['TMPDIR']
elif os.path.exists(SSD_TMP_DIR):
    TMPDIR = SSD_TMP_DIR
else:
    TMPDIR = tempfile.gettempdir()


configfile: "merge.yaml"
SMS = list(config.keys())
Rs={}; Qs={}
for SM in SMS:
	Rs[SM] = os.path.abspath(config[SM]["r"]) 
	Qs[SM] = os.path.abspath(config[SM]["q"]) 

WORKDIR = "quickmerge_out"
workdir: WORKDIR

Ls = [200000]

#
# MOST OF THE MERGED SEQEUCNE COMES FORM THE QUERY
#

wildcard_constraints:
	SM="|".join(SMS),
	L="\d+"

rule all:
	input:
		merged = expand("{SM}_{L}_merged.fasta", SM=SMS, L=Ls),	
		fai = expand("{SM}_{L}_merged.fasta.fai", SM=SMS, L=Ls),	

def get_r(wc):
	SM = str(wc.SM); return( Rs[SM] )
def get_q(wc):
	SM = str(wc.SM); return( Qs[SM] )
def clean_fasta(fastaf, outf, key):
	with pysam.FastxFile(fastaf) as fasta, open(outf, mode='w+') as out:
		counter = 1
		for rec in fasta:
			out.write( ">{}_{}_{:08}\n{}\n".format( rec.name.strip(), key, counter, rec.sequence.upper()) )
			counter += 1

rule fix_fasta:
	input:
		r=get_r,
		q=get_q,
	output:
		r = "temp/{SM}.r.fasta",
		q = "temp/{SM}.q.fasta",
	run:
		clean_fasta(input["r"], output["r"], "r")
		clean_fasta(input["q"], output["q"], "q")
		

rule nucmer:
	input:
		r = rules.fix_fasta.output.r,
		q = rules.fix_fasta.output.q,
	output:
		delta = "nucmer_out/{SM}.delta",
	threads: 64
	shell:"""
mkdir -p nucmer_out
nucmer -t {threads} -l 100 --prefix=nucmer_out/{wildcards.SM} {input.r} {input.q}
"""

rule nuc_filter:
	input:
		delta = rules.nucmer.output.delta,
	output:
		delta = protected("nucmer_out/{SM}.{L}.delta"),
	threads:1
	shell:"""
delta-filter -r -q -l {wildcards.L} {input.delta} > {output.delta}
"""


	
rule quick_merge:
	input:
		r = rules.fix_fasta.output.r,
		q = rules.fix_fasta.output.q,
		delta = rules.nuc_filter.output.delta,
	output:
		aln= "qm_out/aln_summary_{SM}_{L}.tsv",
		anc= "qm_out/anchor_summary_{SM}_{L}.txt",
		param = "qm_out/param_summary_{SM}_{L}.txt",
		merged = "qm_out/merged_{SM}_{L}.fasta",
	threads:1
	shell:"""
D=$(readlink -f {input.delta})
Q=$(readlink -f {input.q})
R=$(readlink -f {input.r})
mkdir -p qm_out && cd qm_out

quickmerge -d $D -q $Q -r $R \
	-hco 5.0 -c 1.5 \
	-l 20000000 \
	-ml 1000000 \
	-p {wildcards.SM}_{wildcards.L}

"""

rule final_fasta:
	input:
		merged = rules.quick_merge.output.merged , 
	output:
		merged = "{SM}_{L}_merged.fasta",
		fai  = "{SM}_{L}_merged.fasta.fai",
	threads:1	
	shell:"""
ln {input.merged} {output.merged}
samtools faidx {output.merged}
"""

# make fake config
if(True):
	import yaml
	qc_cfg = {}
	L=200000
	for SM in SMS:
		qc_cfg[SM] = {"asm": os.path.abspath(f"{SM}_{L}_merged.fasta") }
	open("qc.yaml", "w+").write(yaml.dump(qc_cfg))
	
subworkflow qc:
	workdir:
		WORKDIR 
	snakefile:
		"/net/eichler/vol26/home/mvollger/projects/hifi_asm/qc_asm.smk"
		

rule ideograms:
	input:
		qc( expand("results/{SM}.ideogram.pdf", SM=SMS) )



