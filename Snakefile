import os

configfile: 'config/config.yml'
cpus = config['ncores']
host_ref = config['host_ref']

INPUT_FILE = config['input_file']
OUTPUT = config['output_dir']
SAMPLES = []
STUDIES = []

samp2study = {}

os.system("chmod -R +x scripts")

with open(INPUT_FILE) as f:
    for line in f:
        cols = line.rstrip().split("\t")
        samp2study[cols[0]] = cols[1]
        SAMPLES.append(cols[0])
        STUDIES.append(cols[1])

for sample in samp2study:
    dirname = OUTPUT+"/"+samp2study[sample]+"/"+sample+"/logs"
    if not os.path.exists(dirname):
        os.makedirs(dirname)

rule targets:
    input:
        expand([OUTPUT+"/{study}/{sample}/done.txt"], zip, study=STUDIES, sample=SAMPLES)
               
rule ena_download:
    output:
        fwd = "{output}/{study}/{sample}/{sample}_1.fastq.gz",
        rev = "{output}/{study}/{sample}/{sample}_2.fastq.gz"
    params:
        outdir = "{output}/{sample}"
    conda:
        "config/envs/ena_download.yml"
    resources:
        ncores = cpus
    shell:
        """
        fastq-dl --cpus {resources.ncores} -a {input} --provider ena --only-provider --outdir {params.outdir}
        """

rule metagen_qc:
    input:
        fwd = "{output}/{study}/{sample}/{sample}_1.fastq.gz",
        rev = "{output}/{study}/{sample}/{sample}_2.fastq.gz"
    output:
        "{output}/{study}/{sample}/done.txt"
    params:
        bwa_ref = host_ref,
        fwd = "{output}/{study}/{sample}/{sample}_clean_1.fastq.gz",
        rev = "{output}/{study}/{sample}/{sample}_clean_2.fastq.gz"
    conda:
        "config/envs/metagen_qc.yml"
    resources:
        ncores = cpus
    shell:
        """
        ./scripts/metagen-fastqc.sh -t {resources.ncores} -c {params.bwa_ref} -f {input.fwd} -r {input.rev}
        mv {params.fwd} {input.fwd}; mv {params.rev} {input.rev}
        touch done.txt
        """
