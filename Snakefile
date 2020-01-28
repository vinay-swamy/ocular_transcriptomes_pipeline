'''
config file:
    sampleFile: tab seperated file  containing sample info
    refFasta_url: link to refernece fastq_path
    salmon_version:
    sratoolkit_version:
notes:
-the 5.2 version requires specifying directorys in output section of rule iwth directory(). Biowulf currently using 5.1
-need to make a rule to download all Gencode refs

11/05/19
    REevesied 
-09/20/19
    Revised strategy: merge w/o genome guide with 1 TPM as c/o. quantify with salmon, remove tx with avg count< 1, and high
    bootstrap variance. filter removed transcripts from gtf, and then requantify. and then merge gtfs together.
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

-01/22/19
    - use gffcompare gtf not stringtie-merge gtf because gffcompare is significantly better than stringtie at mapping
    back to genes. GFFcompare found 20K novel tx vs 18K on st-merge, with the same number of transcript.
    - at the initial merge step with stringtie, filteing out transcripts with at least 1 tpm in a third of the samples

01/15/19 Changes
- converted STAR rules into shell from python
- made outputs more organized
- created a synthetic set by sampling 5 samples from all tissues >10 samples, incliding eye tissues, might wanna change that;
    using synthetic set as body set for gtf and rmats
- making mulitple tissue comparisons against the synth for rMATs instead of pair wise comparisons
- restructured outputs so things are a little more organized
- only using paired samples for rmats
- finally got rid of the shitty salmon command
-Reminder that STAR makes the bam even if the alignment fails
-following CHESS paper - run each bam individually through stringtie, then merge to a tissue level, then merge into 1
 use gffcompare at each meerge step;

-12/14/18
    -st-merge at at it default tpm cut off did nothing, so now going to do what chess ppl did and filter it at 1tpm per tissue

-12/13/18
    - tried GFFcompare on all samples first gave 300k tx's but salmon couldn't map to them, so will now use
      stringtie merge on a per tissue level, which will cut out a lot of transcripts, then merge with gffcompare.
    - moved gffread > tx into its own rule


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
            res.append('data/st_filt/{}.gtf'.format(sample))
    return (res)


def salmon_input(id,sample_dict,fql):
    paired=sample_dict[id]['paired']
    id= fql + 'fastq_files/' + id
    if paired:
        return('-1 {s}_1.fastq.gz -2 {s}_2.fastq.gz'.format(s=id))
    else:
        return('-r {}.fastq.gz'.format(id))


def subtissue_to_sample(subtissue, sample_dict):
    res=[]
    [res.append(f'data/salmon_quant/{subtissue}/{sample}/quant.sf') for sample in sample_dict.keys() if sample_dict[sample]['subtissue'] == subtissue ]
    return(res)




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
#python_env=config['python_env']
#commonly used files/paths
working_dir=config['working_dir']
STARindex='ref/STARindex'
ref_tx_fasta='ref/gencode_tx_ref.fa'
ref_GTF='ref/gencode_comp_ano.gtf'
ref_genome='ref/gencode_genome.fa'
fql=config['fastq_path']
bam_path=config['bam_path']
stringtie_full_gtf='data/gtfs/all_tissues.combined.gtf'
hmmer_version=config['hmmer_version']
#win_size=config['window_size']

rule all:
    input:  'data/rdata/novel_exon_classification.Rdata', 'data/rmats/all_tissues_psi.tsv', 'data/all_tissue_quant.Rdata',  'data/novel_loci/hmmer/seq_hits.tsv', 'data/novel_loci/novel_loci_blast_results.tsv', 'data/shiny_data/DNTX_db.sql' 
    #expand('data/gtfs/raw_tissue_gtfs/{subt}.combined.gtf', subt=subtissues)
    # input:stringtie_full_gtf,'data/exp_files/all_tissue_quant.tsv.gz','data/rmats/all_tissues_psi.tsv', 'data/rmats/all_tissues_incCounts.tsv', 'data/seqs/transdecoder_results/best_orfs.transdecoder.pep', 'data/rdata/novel_exon_classification.Rdata'

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
        awk '$3 == "transcript"' ref/gencode_comp_ano.gtf | cut -f1,7,4,5 > ref/gencode_comp_ano_trim.tsv
        '''

rule clean_phylop_and_snps:
    input: snps=expand('ref/snps/bed_chr_{chrom}.bed.gz', chrom=list(range(23))[1:]), pp='ref/phylop_20/hg38.phyloP20way.bw'
    output:snps='ref/snps/hg38.snps.all.sorted.bed.gz', pp='ref/phylop_20/hg38.phyloP20way.sorted.bed.gz'
    shell:
        '''
        rm -rf snp_tmp
        mkdir snp_tmp
        module load ucsc
        module load bedops
        moulde load bedtools
        bigWigToWig {input.pp} | wig2bed | bedtools sort -i - | gzip -c - > {output.pp}
        for i in ref/snps/bed_chr_*.bed.gz ; do zcat $i | tail -n+2   > snp_tmp/snps.bed 
        awk '$2 != $3{OFS="\t"; print $0}' snp_tmp/snps.bed  > snp_tmp/good_snps.bed 
        awk '$2 == $3 {OFS="\t";print $1, $2, $3+1, $4, $5, $6}' snp_tmp/snps.bed | cat - snp_tmp/good_snps.bed | sort-bed --max-mem 24G - | gzip -c - > ref/snps/hg38.snps.all.sorted.bed.gz
        rm -rf snp_tmp
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


#hand added these synth samples.
rule make_rmats_synth:
    input: expand(bam_path+'STARbams/{id}/Sorted.out.bam', id=['SRS648866', 'SRS648919','SRS649535','SRS649567','SRS649619','SRS649622','SRS649657'])
    output:'ref/rmats_locs/synth.rmats.txt'
    shell:
        '''
        echo {input} | tr ' ' ',' > {output}
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

rule run_stringtie:
    input:bam_path + 'STARbams/{sample}/Sorted.out.bam'
    output:'st_out/{sample}.gtf'
    shell:
        '''
        module load {stringtie_version}
        stringtie {input[0]} -o {output[0]} -p 8 -G {ref_GTF}
        '''


'''
****PART 2**** build Transcriptome, and process
'''


#gffread v0.9.12.Linux_x86_64/
rule filter_sample_gtf:
    input:'st_out/{sample}.gtf'
    output:'data/st_filt/{sample}.gtf'
    shell:
        '''
        module load {stringtie_version}
        stringtie --merge -l {wildcards.sample} -T 1 -F0 {input} > {output}
        '''




rule merge_gtfs_by_tissue:
    input:lambda wildcards: subtissue_to_gtf(wildcards.subtissue, sample_dict)
    params: outdir= lambda wildcards: 'data/gtfs/raw_tissue_gtfs/{}'.format(wildcards.subtissue)
    output: 'data/gtfs/raw_tissue_gtfs/{subtissue}.combined.gtf', 'data/gtfs/raw_tissue_gtfs/{subtissue}.tracking'
    shell:
        '''
	    module load {gffcompare_version}
        gffcompare -r {ref_GTF} -p {wildcards.subtissue} -o {params.outdir} {input}
        '''

'''
keep only transcripts that are present in samples form at least three different studies, and correct the transcript id so refernce trasncripts have ENST ids, and novel transcirpts have a novel ID
'''
rule filter_tissue_gtfs_gffcompare:
    input:gtf='data/gtfs/raw_tissue_gtfs/{subtissue}.combined.gtf',track='data/gtfs/raw_tissue_gtfs/{subtissue}.tracking'
    params: repeats = 'ref/UCSC_grch38_repats.bed'
    output: gtf='data/gtfs/raw_tissue_gtfs/{subtissue}.gfcfilt.gtf', det_df = 'data/misc/raw_dd/{subtissue}.dd.tsv.gz'
    shell:
        '''
        module load {R_version}
        Rscript scripts/filter_gtf_gffcompare.R \
            --workingDir {working_dir} \
            --gtfFile {input.gtf} \
            --refGtfFile {ref_GTF} \
            --trackingFile {input.track} \
            --targetTissue {wildcards.subtissue} \
            --sampleTableFile {sample_file} \
            --repeatsFile {params.repeats} \
            --outGtfFile {output.gtf} \
            --outDetDf {output.det_df}
        '''

rule make_tx_fasta:
    input: tool = 'gffread/gffread', gtf = lambda wildcards: f'data/gtfs/raw_tissue_gtfs/{wildcards.subtissue}.gfcfilt.gtf' if wildcards.subtissue != 'all_tissues.combined' else 'data/gtfs/all_tissues.combined.gtf'
    output: 'data/seqs/{subtissue}_tx.fa'
    shell:
        '''
        ./gffread/gffread -w {output} -g {ref_genome}  {input.gtf}
        '''

rule build_salmon_index:
    input: 'data/seqs/{subtissue}_tx.fa'
    output: directory('data/salmon_indices/{subtissue}')
    shell:
        '''
        module load {salmon_version}
        salmon index -t {input} -i {output} --type quasi --perfectHash -k 31
        '''

rule run_salmon:
    input: fastqs=lambda wildcards: [fql+'fastq_files/{}_1.fastq.gz'.format(wildcards.sampleID),fql+'fastq_files/{}_2.fastq.gz'.format(wildcards.sampleID)] if sample_dict[wildcards.sampleID]['paired'] else fql+'fastq_files/{}.fastq.gz'.format(wildcards.sampleID),
        index='data/salmon_indices/{subtissue}'
    params: cmd=lambda wildcards: salmon_input(wildcards.sampleID,sample_dict,fql),\
     outdir=lambda wildcards: 'data/salmon_quant/{}/{}'.format(wildcards.subtissue, wildcards.sampleID)
    output: 'data/salmon_quant/{subtissue}/{sampleID}/quant.sf','data/salmon_quant/{subtissue}/{sampleID}/quant_bootstraps.tsv.gz'
    shell:
        '''
        id={wildcards.sampleID}
        module load {salmon_version}
        salmon quant -p 8 -i {input.index} -l A --gcBias --seqBias --numBootstraps 100 --validateMappings {params.cmd} -o {params.outdir}
        python3 scripts/convertBootstrapsToTsv.py {params.outdir} {params.outdir}
        '''



'''
Use the variance of salmon quantification to remove transcripts that have a high variability in their quantification
I noticed that reference transcripts have a much lower variability than novel ones, so we are going to remove anything 
abouve the 95th percentile of ref quant variance

'''

rule filter_gtf_salmonvar:
    input: salmon_quant= lambda wildcards:subtissue_to_sample(wildcards.subtissue, sample_dict), gtf='data/gtfs/raw_tissue_gtfs/{subtissue}.gfcfilt.gtf'
    params: quant_path=lambda wildcards: f'data/salmon_quant/{wildcards.subtissue}/'
    output:gtf='data/gtfs/raw_tissue_gtfs/{subtissue}.compfilt.gtf'
    shell:
        '''
        module load {R_version}
        Rscript scripts/filter_gtf_salmonvar.R \
            --workingDir {working_dir} \
            --pathToQuantFiles {params.quant_path} \
            --gtfFile {input.gtf} \
            --outGtfFile {output.gtf} 
        '''


'''
use gfcompare to combine all tissue specifc gtfs to make the all_tissues.cmobined final gtf, then use the script to parse the tracking table, and fix the transcript id's similar to before.

'''


rule merge_tissue_gtfs:
    input: gtfs=expand('data/gtfs/raw_tissue_gtfs/{subtissue}.compfilt.gtf',subtissue=subtissues)
    params: gffc_prefix='all_tissues'
    output: gtf='data/gtfs/all_tissues.combined.gtf', raw_track_file='data/gffcomp_dir/all_tissues.tracking', tx_converter_tab='data/misc/TCONS2MSTRG.tsv'
    shell:
        '''
        mkdir -p data/gffcomp_dir
        module load {gffcompare_version}
        gffcompare -r {ref_GTF} --strict-match  -p DNTX -o data/gffcomp_dir/{params.gffc_prefix} {input.gtfs}
        module load {R_version}
        Rscript scripts/TCONS_to_tissueMSTRG.R \
            --workingDir {working_dir} \
            --sampleTableFile {sample_file} \
            --trackFile {output.raw_track_file} \
            --rawGtfFile data/gffcomp_dir/all_tissues.combined.gtf \
            --refGtfFile {ref_GTF} \
            --tcons2mstrgFile {output.tx_converter_tab} \
            --outGtfFile {output.gtf} 
        '''


'''
fix the tissue specific gtfs so they have the all_tissue novel transcript ids(DNTX); make salmon exp Rdata with the new tissue specific gtf

'''

rule track_CHESS_transcripts:
    input: base=expand('data/gtfs/raw_tissue_gtfs/{subtissue}.combined.gtf', subtissue=subtissues), \
    gfcfilt=expand('data/gtfs/raw_tissue_gtfs/{subtissue}.gfcfilt.gtf', subtissue=subtissues), \
    compfilt=expand('data/gtfs/raw_tissue_gtfs/{subtissue}.compfilt.gtf', subtissue=subtissues), \
    final=expand('data/gtfs/final_gtfs/{subtissue}.gtf', subtissue=subtissues)
    output:'data/union_gtfs/final.combined.gtf'
    shell:
        '''
        module load {gffcompare_version}
        mkdir -p data/union_gtfs/
        gffcompare -r {ref_GTF} -o data/union_gtfs/base {input.base}
        gffcompare -r {ref_GTF} -o data/union_gtfs/gfcfilt {input.gfcfilt}
        gffcompare -r {ref_GTF} -o data/union_gtfs/compfilt {input.compfilt}
        gffcompare -r {ref_GTF} -o data/union_gtfs/final {input.final}
        '''



rule clean_tissue_gtfs_clean_salmon_quant:
    input:tissue_gtf='data/gtfs/raw_tissue_gtfs/{subtissue}.compfilt.gtf',  salmon_quant= lambda wildcards: subtissue_to_sample(wildcards.subtissue,  sample_dict), tx_converter_tab='data/misc/TCONS2MSTRG.tsv', det_df='data/misc/raw_dd/{subtissue}.dd.tsv.gz'
    params: quant_path=lambda wildcards: f'data/salmon_quant/{wildcards.subtissue}/'
    output:gtf='data/gtfs/final_gtfs/{subtissue}.gtf', exp_file='data/exp_files/{subtissue}.RDS', out_detdf='data/misc/final_dd/{subtissue}.dd.tsv'
    shell:
        '''
        module load {R_version}
        Rscript scripts/fix_tissue_gtf_txids_make_expfiles.R \
            --workingDir {working_dir} \
            --tissueGtfFile {input.tissue_gtf} \
            --targetTissue {wildcards.subtissue} \
            --quantPath {params.quant_path} \
            --trackFile {input.tx_converter_tab} \
            --detDfFile {input.det_df} \
            --expFile {output.exp_file} \
            --outGtfFile {output.gtf} \
            --outDetFile {output.out_detdf}
        '''


'''
Finally, merge the tissue specific quant into an all tissues quant

'''

rule merge_all_salmon_quant:
    input:expand('data/exp_files/{subtissue}.RDS', subtissue= subtissues),tx_converter_tab='data/misc/TCONS2MSTRG.tsv'
    params: quant_path='data/exp_files/'
    output: 'data/all_tissue_quant.Rdata'
    shell:
        '''
        module load {R_version}
        Rscript scripts/merge_salmon_quant.R  --workingDir {working_dir} --pathToQuant {params.quant_path} --outExpFile {output}
        '''

'''
rMATs part.Moo
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
              awk ' $2 == "y" {{print $1}}'  |\
              sed -e "s|^|$prefix|g" - |\
              sed -e "s|$|$suffix|g" - |\
              tr '\\n' ',' > {params.bam_dir}/ref/rmats_locs/$t.rmats.txt
        done
        '''

rule runrMATS:
    input: loc = 'ref/rmats_locs/{subtissue}.rmats.txt', idx = 'ref/STARindex', gtf = 'data/gtfs/final_gtfs/{subtissue}.gtf', synthfile='ref/rmats_locs/synth.rmats.txt' 
    params: outdir= lambda wildcards: f'data/rmats_out/{wildcards.subtissue}'
    output: expand('data/rmats_out/{{subtissue}}/{event}.MATS.JC.txt', event=rmats_events)
    # might have to change read length to some sort of function
    shell:
        '''
        subtissue={wildcards.subtissue}
        module load {rmats_version}
        rmats --b1 {input.loc} --b2 {input.synthfile}  -t paired  \
        --nthread 8  --readLength 130 --gtf {input.gtf} --bi {input.idx} --od {params.outdir}
        '''

rule process_rmats_output:
    input: expand('data/rmats_out/{sub_tissue}/{event}.MATS.JC.txt', sub_tissue=[x for x in subtissues if x != 'Cornea_Fetal.Tissue'], event= rmats_events)
    params: rmats_od='data/rmats_out/', rm_locdir='ref/rmats_locs/'
    output: psi='data/rmats/all_tissues_psi.tsv', incCounts='data/rmats/all_tissues_incCounts.tsv'
    shell:
        '''
        module load {R_version}
        Rscript scripts/process_rmats_output.R \
            --workingDir {working_dir} \
            --sampleTableFile {sample_file} \
            --rmatsOutDir {params.rmats_od} \
            --rmLocDir {params.rm_locdir} \
            --psiOutFile {output.psi} \
            --incCountsOutFile {output.incCounts}
        '''


'''
pulled this from the trinnotate pipeline, makes a gff of the of using protein translations
'''


rule run_trans_decoder:
    input:'data/seqs/all_tissues.combined_tx.fa'
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
        TransDecoder.LongOrfs -m 60 -t transcripts.fasta
        TransDecoder.Predict --single_best_only  -t transcripts.fasta
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

'''

have to manually download this from UCSC
table browser > group=repeats, track=RepeatMasker
'''

rule catagorize_novel_exons:
    input: gtf=stringtie_full_gtf, gff3='data/seqs/transdecoder_results/all_tissues.combined_transdecoderCDS.gff3', full_pep='data/seqs/transdecoder_results/best_orfs.transdecoder.pep'
    params: repeats='ref/UCSC_grch38_repats.bed', novel_loci_txids='data/novel_loci/novel_loci_txids.txt'
    output: classfile='data/rdata/novel_exon_classification.Rdata', gtfano='data/gtfs/all_tissues.combined_NovelAno.gtf',                   novel_loci_pep='data/novel_loci/novel_loci.pep'
    shell:
        '''
        module load {R_version}
        module load {bedtools_version}
        Rscript scripts/classify_novel_exons.R \
            --workingDir {working_dir} \
            --gfcGtfFile {stringtie_full_gtf} \
            --sampleTableFile {sample_file} \
            --repeatBedFile {params.repeats} \
            --gff3File {input.gff3} \
            --exonClassOutFile {output.classfile} \
            --gtfAnoOutFile {output.gtfano}\
            --novelLociTxIds {params.novel_loci_txids}
        python3 scripts/select_entry_from_fasta.py --infasta {input.full_pep} --txToKeep {params.novel_loci_txids} --outfasta {output.novel_loci_pep}
        '''



rule prep_shiny_data:
    input: ano_gtf = 'data/gtfs/all_tissues.combined_NovelAno.gtf', \
        det_dfs = expand('data/misc/final_dd/{subtissue}.dd.tsv', subtissue=subtissues), \
        tc2m = 'data/misc/TCONS2MSTRG.tsv', \
        all_exp_file='data/all_tissue_quant.Rdata', \
        snps='ref/snps/hg38.snps.all.sorted.bed.gz', \
        pp='ref/phylop_20/hg38.phyloP20way.sorted.bed.gz', \
        gff3='data/seqs/transdecoder_results/all_tissues.combined_transdecoderCDS.gff3'
    params: dd_stem='data/misc/final_dd/REPLACE.dd.tsv'
    output:db_file='data/shiny_data/DNTX_db.sql'  , rdata_file='data/shiny_data/shiny_data.Rdata'
    shell:
        '''
        rm -rf data/shiny_data/debug
        mkdir -p data/shiny_data/debug
        touch data/shiny_data/debug/no_cds_bu_marked_as_pc.txt
        touch data/shiny_data/debug/multi_tx_under_same_id.txt
        touch data/shiny_data/debug/failed_to_find_CDS_start_or_end.txt
        touch data/shiny_data/debug/special_is_funky.txt
        touch data/shiny_data/debug/bad_genes.txt
        module load {R_version}
        module load {bedtools_version}
        Rscript scripts/prep_data_for_shiny.R \
            --workingDir {working_dir} \
            --gtfFile {input.ano_gtf} \
            --tcons2mstrgFile {input.tc2m} \
            --quantFile {input.all_exp_file} \
            --sampleTableFile {sample_file} \
            --ddStem {params.dd_stem} \
            --snpsBed {input.snps} \
            --phylopBed {input.pp} \
            --gff3File {input.gff3}
            --outRdata {output.rdata_file} \
            --outDbFile {output.db_file}

        '''



rule blastp_novel_loci:
    input: pep='data/novel_loci/novel_loci.pep'
    output: results='data/novel_loci/novel_loci_blast_results.tsv'
    shell:
        '''
        module load blast
        blastp -query {input.pep} -db /fdb/blastdb/swissprot  -max_target_seqs 250 -max_hsps 3 -outfmt "6 qseqid sseqid qlen slen nident evalue bitscore"  -num_threads 8 > {output.results}
        '''

rule build_pfm_hmmDB:
     params: url=config['pfam_db']
     output:'ref/hmmer/Pfam-A.hmm'
     shell:
         '''
         wget -O - {params.url} | gunzip -c - > {output}
         module load {hmmer_version}
         hmmpress {output}
         '''


rule run_hmmscan:
     input: pfam='ref/hmmer/Pfam-A.hmm',  pep='data/novel_loci/novel_loci.pep'
     output:tab='data/novel_loci/hmmer/seq_hits.tsv',
         dom='data/novel_loci/hmmer/domain_hits.tsv',
         pfm='data/novel_loci/hmmer/pfam_hits.tsv'
     shell:
         '''
         module load {hmmer_version}
         hmmscan --cpu 24 --tblout {output.tab} --domtblout {output.dom} --pfamtblout {output.pfm} {input.pfam} {input.pep}
         '''



