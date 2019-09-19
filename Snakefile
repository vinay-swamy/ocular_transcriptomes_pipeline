'''
config file:
    sampleFile: tab seperated file  containing sample info
    refFasta_url: link to refernece fastq_path
    salmon_version:
    sratoolkit_version:
notes:
-the 5.2 version requires specifying directorys in output section of rule iwth directory(). Biowulf currently using 5.1
-need to make a rule to download all Gencode refs

06/10/19 changes
- I added way to much downstream stuff to the pipeline, assuming that the splicing/trnascriptome stuff was more accurate
  than it actually was, so I removed a lot stuff,  hmmer, converting to b37, getting a gff3, to focus more on improving
  the accuracy of the exon detections
- using the comprehensive gencode annotation for everything
- no longer realigning to stringtie gtf, this caused a lot of problems, an it makes a little more sense over all now,
  as both tools are blind to each other
- strategy for imporving - use a model trained on exons that get longer/shorter present in reference annotation
    - a model is trained for each sample, and then potentially discovered exons are run through it for each smaple
    - the model uses the bp level count info at the boundary of an exon getting longer/shorter.
    - more often an exon is found th emore likely its actually there
- after I get this working, then ill add the downstream stuff back later

02/12/19 Changes
-all software versioned in config file
-added hmmscan protein domain search

01/15/19 Changes
- converted STAR rules into shell from python
- made outputs more organized
- created a synthetic set by sampling 5 samples from all tissues >10 samples, incliding eye tissues, might wanna change that;
    using synthetic set as body set for gtf and rmats
- making mulitple tissue comparisons against the synth for rMATs instead of pair wise comparisons
- restructured outputs so things are a little more organized
- only using paired samples for rmats
- finally got rid of the shitty salmon command
**************REMEMBER TO CHANGE THE WORKING DIR IN THE CONFIG FILE IF YOU RERUN*********
'''
import subprocess as sp

def readSampleFile(samplefile):
    # returns a dictionary of dictionaries where first dict key is sample id and second dict key are sample  properties
    res={}
    with open(samplefile) as file:
        file.readline()#skip the first line because it has a header
        for line in file:
            info=line.strip('\n').split('\t')
            res[info[0]]={'paired':True if info[1]=='y' else False, 'tissue':info[2],'subtissue':info[3]}
    return(res)

def subtissue_to_gtf(subtissue, sample_dict):
    res=[]
    for sample in sample_dict.keys():
        if sample_dict[sample]['subtissue']==subtissue :
            res.append('st_out/{}.gtf'.format(sample))
    return (res)
def subtissue_to_sample_agg(subtissue, sample_dict):
    res=[]
    for sample in sample_dict.keys():
        if sample_dict[sample]['subtissue'] == subtissue:
            res.append('data/rawST_tx_quant_files/{}/{}/quant.sf'.format(subtissue, sample))
    return(res)

def subtissue_to_sample_all( sample_dict):
    res=[]
    for sample in sample_dict.keys():
        subtissue=sample_dict[sample]['subtissue']
        res.append('data/filter_tx_quant_files/{}/{}/quant.sf'.format(subtissue, sample))
    return(res)

def salmon_input(id,sample_dict,fql):
    paired=sample_dict[id]['paired']
    id= fql + 'fastq_files/' + id
    if paired:
        return('-1 {s}_1.fastq.gz -2 {s}_2.fastq.gz'.format(s=id))
    else:
        return('-r {}.fastq.gz'.format(id))
def build_to_fasta_file(build):
    if build == 'gencode':
        return('ref/gencode_tx_ref.fa')# hard coded for now, will need to add more later
def build_subtissue_lookup(subtissue, build=''):
    if subtissue == 'all_tissues.combined':
        return 'data/gtfs/all_tissues.combined.gtf'
    else:
        if build == 'rawST_tx_quant_files':
            return 'data/gtfs/raw_tissue_gtfs/{}_st.gtf'.format(subtissue)
        else :
            return 'data/gtfs/filtered_tissue/{}.gtf'.format(subtissue)



#sample information
sample_file=config['sampleFile']
sample_dict=readSampleFile(config['sampleFile'])# sampleID:dict{path,paired,metadata}
tissue_file=config['tissueFile']
subtissue_file=config['subtissueFile']
with open(tissue_file) as tf, open(subtissue_file) as sf:
    tissues= [line.strip('\n') for line in tf]
    subtissues= [line.strip('\n') for line in sf]
sample_names=sample_dict.keys()
rmats_events=['SE','RI','MXE','A5SS','A3SS']
eye_tissues=['Retina_Adult.Tissue', 'Retina_Fetal.Tissue', 'RPE_Adult.Tissue', 'RPE_Fetal.Tissue', 'Cornea_Adult.Tissue','Cornea_Fetal.Tissue']
#software version info
salmon_version=config['salmon_version']
stringtie_version=config['stringtie_version']
STAR_version=config['STAR_version']
rmats_version=config['rmats_verson']
R_version=config['R_version']
TransDecoder_version=config['TransDecoder_version']
samtools_version=config['samtools_version']
gffcompare_version=config['gffcompare_version']
mosdepth_version=config['mosdepth_version']
bedtools_version=config['bedtools_version']
python_env=config['python_env']
#commonly used files/paths
working_dir=config['working_dir']
STARindex='ref/STARindex'
ref_tx_fasta='ref/gencode_tx_ref.fa'
ref_GTF='ref/gencode_comp_ano.gtf'
ref_genome='ref/gencode_genome.fa'
fql=config['fastq_path']
bam_path=config['bam_path']
stringtie_full_gtf='data/gtfs/all_tissues.combined.gtf'
win_size=config['window_size']

rule all:
    input:'data/rmats/all_tissues_psi.tsv', 'data/rmats/all_tissues_incCounts.tsv', stringtie_full_gtf,\
    'data/exp_files/all_tissues_complete_quant.rdata','data/seqs/transdecoder_results/best_orfs.transdecoder.pep',\
    'rdata/novel_exon_classification.rdata'

'''
****PART 1**** download files and align to genome
-still need to add missing fastq files
-gffread needs indexed fasta
-need to add versioning of tools to yaml{DONE}
04/08/2019 - added the mysql command to get a comprehensive refseq gtf from ucsc, which is somehow different than the one
from ncbi. see https://bioinformatics.stackexchange.com/questions/2548/hg38-gtf-file-with-refseq-annotations
'''
rule downloadAnnotation:
    output: reftxfa=ref_tx_fasta, refgtf=ref_GTF, ref_gen=ref_genome,prot_seq='ref/gencodeProtSeq.fa',\
     gencode_gff='ref/gencodeGFF3.gff', ensbl_gtf='ref/ensembl_ano.gtf', refseq='ref/refseq_ncbi.gff3', ucsc='ref/ucsc.gtf'
    shell:
        '''
        wget -O - {config[ref_tx_fasta_url]} | gunzip -c - > {output.reftxfa}
        wget -O - {config[refGTF_url]} | gunzip -c - > {output.refgtf}
        wget -O - {config[ref_genome_url]} | gunzip -c - > /tmp/gencodePA_tmp.fa
        wget -O - {config[refProtSeq_url]} | gunzip -c - > /tmp/gencodeProtSeq.fa
        wget -O - {config[refGFF3_url]} | gunzip -c > {output.gencode_gff}
        wget -O - {config[ensembl_gtf_url]} | gunzip -c - > {output.ensbl_gtf}
        wget -O - {config[refseq_ncbi_url]} | gunzip -c - > {output.refseq}
        module load mysql
        module load ucsc
        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from refGene" hg38 |\
         cut -f2- | genePredToGtf -source=hg38.refGene.ucsc file stdin {output.ucsc}
        module load python/3.6
        python3 scripts/filterFasta.py /tmp/gencodePA_tmp.fa ref/chroms_to_remove {output.ref_gen}
        python3 scripts/clean_fasta.py /tmp/gencodeProtSeq.fa {output.prot_seq}
        module load {samtools_version}
        samtools faidx ref/gencode_genome.fa
        '''

# This is manily for rerunning on biowulf,
rule build_gffread:
    output:'gffread/gffread'
    shell:
        '''
        git clone https://github.com/gpertea/gclib
        git clone https://github.com/gpertea/gffread
        cd gffread
        make release
        '''



rule build_STARindex:
    input: ref_genome, ref_GTF
    output:directory(STARindex)
    shell:
        '''
        module load {STAR_version}
        mkdir -p ref/STARindex
        STAR --runThreadN 16 --runMode genomeGenerate --genomeDir {output[0]} --genomeFastaFiles {input[0]} --sjdbGTFfile {input[1]} --sjdbOverhang 100

        '''

rule run_STAR_alignment:
    input: fastqs=lambda wildcards: [fql+'fastq_files/{}_1.fastq.gz'.format(wildcards.id),fql+'fastq_files/{}_2.fastq.gz'.format(wildcards.id)] if sample_dict[wildcards.id]['paired'] else fql+'fastq_files/{}.fastq.gz'.format(wildcards.id),
        index=STARindex
    output:temp(bam_path+'STARbams/{id}/raw.Aligned.out.bam'), bam_path+'STARbams/{id}/raw.Log.final.out'
    shell:
        '''
        id={wildcards.id}
        mkdir -p STARbams/$id
        module load {STAR_version}
        STAR --runThreadN 8 --genomeDir {input.index} --outSAMstrandField intronMotif  --readFilesIn {input.fastqs} \
        --readFilesCommand gunzip -c --outFileNamePrefix STARbams/$id/raw. --outSAMtype BAM Unsorted
        '''

rule sort_bams:
    input:bam_path + 'STARbams/{id}/raw.Aligned.out.bam'
    output:bam_path + 'STARbams/{id}/Sorted.out.bam'
    shell:
        '''
        module load {samtools_version}
        samtools sort -o {output[0]} --threads 7 {input[0]}
        samtools index -b {output}
        '''
rule calculate_cov:
    input:bam_path+'STARbams/{id}/Sorted.out.bam'
    output: 'coverage_files/{id}/cov.per-base.bed.gz'
    shell:
        '''
        module load mosdepth
        sample={wildcards.id}
        mosdepth coverage_files/$sample/cov {input[0]}
        '''

'''
****PART 2**** build Transcriptome, and process
-Reminder that STAR makes the bam even if the alignment fails
-following CHESS paper - run each bam individually through stringtie, then merge to a tissue level, then merge into 1
 use gffcompare at each meerge step;
-12/13/18
    - tried GFFcompare on all samples first gave 300k tx's but salmon couldn't map to them, so will now use
      stringtie merge on a per tissue level, which will cut out a lot of transcripts, then merge with gffcompare.
    - moved gffread > tx into its own rule
-12/14/18
    -st-merge at at it default tpm cut off did nothing, so now going to do what chess ppl did and filter it at 1tpm per tissue
-01/22/19
    - use gffcompare gtf not stringtie-merge gtf because gffcompare is significantly better than stringtie at mapping
    back to genes. GFFcompare found 20K novel tx vs 18K on st-merge, with the same number of transcript.
    - at the initial merge step with stringtie, filteing out transcripts with at least 1 tpm in a third of the samples

'''


rule run_stringtie:
    input:bam_path + 'STARbams/{sample}/Sorted.out.bam'
    output:'st_out/{sample}.gtf'
    shell:
        '''
        module load {stringtie_version}
        stringtie {input[0]} -o {output[0]} -p 8 -G {ref_GTF}
        '''
#gffread v0.9.12.Linux_x86_64/
rule merge_gtfs_by_tissue:
    input: lambda wildcards: subtissue_to_gtf(wildcards.subtissue, sample_dict)
    output: 'data/gtfs/raw_tissue_gtfs/{subtissue}_st.gtf'
    shell:
        '''

	    module load {stringtie_version}
        stringtie --merge  -l {wildcards.subtissue}_MSTRG  -T 1  -F 0 {input} |\
         tr '*' '+' > {output}
        '''

'''
This chunk runs twice, where we first quantify the raw transcriptome, then filter, rebuild and requantify

'''
#************************************************************************************************************************
rule make_tx_fasta:
    input: tool='gffread/gffread',gtf= lambda wildcards: build_subtissue_lookup(wildcards.subtissue, wildcards.build)
    output: 'data/seqs/{build}/{subtissue}_tx.fa'
    shell:
        '''
        ./gffread/gffread -w {output} -g {ref_genome}  {input.gtf}
        '''


rule build_salmon_index:
    input: 'data/seqs/{build}/{subtissue}_tx.fa'
    output: directory('data/salmon_indices/{build}/{subtissue}')
    shell:
        '''
        module load {salmon_version}
        salmon index -t {input} -i {output} --type quasi --perfectHash -k 31
        '''


rule run_salmon:
    input: fastqs=lambda wildcards: [fql+'fastq_files/{}_1.fastq.gz'.format(wildcards.sampleID),fql+'fastq_files/{}_2.fastq.gz'.format(wildcards.sampleID)] if sample_dict[wildcards.sampleID]['paired'] else fql+'fastq_files/{}.fastq.gz'.format(wildcards.sampleID),
        index='data/salmon_indices/{build}/{subtissue}'
    params: cmd=lambda wildcards: salmon_input(wildcards.sampleID,sample_dict,fql),\
     outdir=lambda wildcards: 'data/{}/{}/{}'.format(wildcards.build, wildcards.subtissue, wildcards.sampleID)
    output: 'data/{build}/{subtissue}/{sampleID}/quant.sf','data/{build}/{subtissue}/{sampleID}/quant_bootstraps.tsv.gz'
    shell:
        '''
        id={wildcards.sampleID}
        module load {salmon_version}
        salmon quant -p 8 -i {input.index} -l A --gcBias --seqBias --numBootstraps 100  {params.cmd} -o {params.outdir}
        python3 scripts/convertBootstrapsToTsv.py {params.outdir} {params.outdir}
        '''
#***************************************************************************************************************************************
rule aggregate_salmon_counts_and_filter_gtf:
    input: qfiles= lambda wildcards: subtissue_to_sample_agg(wildcards.subtissue, sample_dict),\
     gtf='data/gtfs/raw_tissue_gtfs/{subtissue}_st.gtf'
    params: qfolder= lambda wildcards: 'data/rawST_tx_quant_files/'+ wildcards.subtissue
    output: 'data/exp_files/{subtissue}_tx_quant.tsv.gz', 'data/gtfs/filtered_tissue/{subtissue}.gtf'
    shell:
        '''
        module load {R_version}
        Rscript scripts/aggCounts_filterGtf.R {working_dir} {input.gtf} {params.qfolder} {output}
        '''



#***************************************************************************************************************************************
'''
rMATs part.
'''
rule preprMats_running:
    input: expand(bam_path + 'STARbams/{id}/Sorted.out.bam',id=sample_names),
    params: bam_dir=bam_path
    output:expand('ref/rmats_locs/{subtissue}.rmats.txt',subtissue=subtissues)
    shell:
        #include trailing / for bam_dir
        '''
        mkdir -p ref/rmats_locs
        cat ref/subtissues.txt | while read t
        do
            prefix={params.bam_dir}/STARbams/
            suffix=/Sorted.out.bam
            grep $t {sample_file} |\
              awk ' $3 == "y" {{print $1}}'  |\
              sed -e "s|^|$prefix|g" - |\
              sed -e "s|$|$suffix|g" - |\
              tr '\\n' ',' > {params.bam_dir}/ref/rmats_locs/$t.rmats.txt
        done
        '''

rule runrMATS:
    input: loc='ref/rmats_locs/{subtissue}.rmats.txt',idx='ref/STARindex', gtf='data/gtfs/filtered_tissue/{subtissue}.gtf'
    output:expand('rmats_out/{{subtissue}}/{event}.MATS.JC.txt', event=rmats_events)
    # might have to change read length to some sort of function
    shell:
        '''
        subtissue={wildcards.subtissue}
        module load {rmats_version}
        rmats --b1 {input.loc} --b2 ref/rmats_locs/synth.rmats.txt  -t paired  \
        --nthread 8  --readLength 130 --gtf {input.gtf} --bi {input.idx} --od rmats_out/$subtissue
        '''

rule process_rmats_output:
    input: expand('rmats_out/{sub_tissue}/{event}.MATS.JC.txt', sub_tissue=subtissues, event= rmats_events)
    params: rmats_od='rmats_out/', rm_locdir='ref/rmats_locs/'
    output: 'data/rmats/all_tissues_psi.tsv', 'data/rmats/all_tissues_incCounts.tsv'
    shell:
        '''
        module load {R_version}
        Rscript scripts/process_rmats_output.R {working_dir} {sample_file} {params.rmats_od} {params.rm_locdir} {output}
        '''
#***************************************************************************************************************************************



#gffread v0.9.12.Linux_x86_64/

rule merge_tissue_gtfs:
    input: expand('data/gtfs/filtered_tissue/{subtissue}.gtf',subtissue=subtissues)
    output: stringtie_full_gtf, 'data/gtfs/all_tissues.stringtie_merge.gtf', 'data/gffcomp_dir/all_tissues.tracking'
    shell:
        '''
        module load {stringtie_version}
        stringtie --merge -G {ref_GTF}  -o {output[1]} {input}
        mkdir -p data/gffcomp_dir
        module load {gffcompare_version}
        gffcompare -r {ref_GTF} -o data/gffcomp_dir/all_tissues {input}
        module load {R_version}
        Rscript scripts/fix_gene_id.R {working_dir} data/gffcomp_dir/all_tissues.combined.gtf {ref_GTF} {output[0]}
        '''


rule merge_filtered_salmon_quant:
    input: qfiles=subtissue_to_sample_all(sample_dict), track_file='data/gffcomp_dir/all_tissues.tracking'
    params: qdir='data/filter_tx_quant_files'
    output: 'data/exp_files/all_tissues_complete_quant.rdata', 'data/misc/gfc_TCONS_to_st_MSTRG.tsv'
    shell:
        '''
        module load {R_version}
        Rscript scripts/TCONS_to_tissueMSTRG.R {working_dir} {params.qdir} {sample_file} {input.track_file} {output}
        '''

rule run_trans_decoder:
    input:'data/seqs/all/all_tissues.combined_tx.fa'
    output:'data/seqs/transdecoder_results/all_tissues.combined_transdecoderCDS.gff3', \
    'data/seqs/transdecoder_results/transcripts.fasta.transdecoder.pep'
    shell:
        '''
        rm -rf TransDecoder
        git clone https://github.com/TransDecoder/TransDecoder.git
        cd TransDecoder
        module load {TransDecoder_version}
        mkdir -p ../data/seqs/transdecoder_results/
        ./util/gtf_genome_to_cdna_fasta.pl ../data/gtfs/all_tissues.combined.gtf ../ref/gencode_genome.fa > transcripts.fasta
        ./util/gtf_to_alignment_gff3.pl ../data/gtfs/all_tissues.combined.gtf > transcripts.gff3
        TransDecoder.LongOrfs -t transcripts.fasta
        TransDecoder.Predict --single_best_only -t transcripts.fasta
        ./util/cdna_alignment_orf_to_genome_orf.pl \
            transcripts.fasta.transdecoder.gff3 \
            transcripts.gff3 \
            transcripts.fasta > ../data/seqs/transdecoder_results/all_tissues.combined_transdecoderCDS.gff3
        mv transcripts.fasta.transdecoder.*  ../data/seqs/transdecoder_results/
        '''
#
rule clean_pep:
    input:'data/seqs/transdecoder_results/transcripts.fasta.transdecoder.pep'
    output:pep='data/seqs/transdecoder_results/best_orfs.transdecoder.pep', meta_info='data/seqs/transdecoder_results/pep_fasta_meta_info.tsv'#, len_cor_tab='data/seqs/len_cor_tab.tsv'
    shell:
        '''
        python3 scripts/clean_pep.py {input} {output.pep} {output.meta_info}
        '''
                #python3 scripts/fix_prot_seqs.py /tmp/tmpvs.fasta  {output.pep} {output.len_cor_tab}
rule catagorize_novel_exons:
    input: stringtie_full_gtf
    output: 'rdata/novel_exon_classification.rdata'
    shell:
        '''
        module load {R_version}
        Rscript scripts/classify_novel_exons.R {working_dir} {stringtie_full_gtf} {sample_file} {output}
        '''


#q
# '''
# part? prep for ML step
#df
# '''
#
# rule makeExonBeds:
#     input: 'rdata/{build}_tx_quant.Rdata'
#     output:'data/bed_files/{build}_alternative_exons.bed'
#     shell:
#         '''
#         module load {bedtools_version}
#         bash scripts/merge_beds_distinct.sh {input} {output}
#         '''
# rule intersect_and_spread:
#     input:'data/bed_files/gencode_alternative_exons.bed', 'coverage_files/{sample}/cov.per-base.bed.gz'
#     output:'data/cleaned_cov/{sample}_bp_features.tsv.gz'
#     shell:
#         '''
#         module load bedtools
#         cut -f1,2,3,4 {input[0]} |\
#          bedtools intersect -loj -a {input[1]} -b stdin |\
#          awk ' $6 != "-1"' - |\
#          python3 scripts/makePerBaseFeatureTable.py {working_dir} {output}
#         '''
# rule train_base_model:
#     input: 'data/cleaned_cov/{sample}_bp_features.tsv.gz'
#     output: 'models/{sample}_xgb_trd.pck'
#     shell:
#         '''
#         python3 scripts/train_model.py {working_dir} {input} {output}
#         '''
