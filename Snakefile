'''
config file:
    sampleFile: tab seperated file  containing sample info
    refFasta_url: link to refernece fastq_path
    salmon_version:
    sratoolkit_version:
notes:
-the 5.2 version requires specifying directorys in output section of rule iwth directory(). Biowulf currently using 5.1
-need to make a rule to download all Gencode refs

***IF YOU CHANGE A RULE NAME MAKE SURE TO CHECK cluster.json ****
Things to do
-rewrite sonneson_low_usage
-rewrite script for removing_tx

'''
import subprocess as sp
import itertools as it


def readSampleFile(samplefile):
    # returns a dictionary of dictionaries where first dict key is sample id and second dict key are sample  properties
    res={}
    with open(samplefile) as file:
        for line in file:
            info=line.strip('\n').split('\t')
            res[info[0]]={'files':info[1].split(','),'paired':True if info[2]=='y' else False, 'tissue':info[3],'subtissue':info[4]}
    return(res)


def genrMATsinput(subtissue,type):
    all_combs= list(it.combinations(subtissue,2))
    res=['rmats_out/{}_VS_{}'.format(x[0],x[1]) for x in all_combs]
    return(res)


def subtissue_to_bam(subtissue,type ,sample_dict):
    res=[]
    # if subtissue=='body':
    #     with open(config['synth_body']) as sb:
    #         res=['STARbams_realigned/'+line.strip()+'/Aligned.out.bam' for line in sb]
    #     return(res)
    for sample in sample_dict.keys():
        if sample_dict[sample]['subtissue']==subtissue and sample_dict[sample]['paired']:
            res.append('STARbams_realigned/{}/Aligned.out.bam'.format(sample))
    #if not res:

    return (res)
#does string tie look
def tissue_to_gtf(tissue, sample_dict):
    res=['st_out/' +sample+'.gtf' for sample in sample_dict.keys() if sample_dict[sample]['tissue']==tissue ]
    return(res)



configfile:'config.yaml'
sample_dict=readSampleFile(config['sampleFile'])# sampleID:dict{path,paired,metadata}
#subtissues_SE=["RPE_Stem.Cell.Line","RPE_Cell.Line","Retina_Adult.Tissue","RPE_Fetal.Tissue","Cornea_Adult.Tissue","Cornea_Fetal.Tissue","Cornea_Cell.Line","Retina_Stem.Cell.Line",'body']
subtissues_PE=["Retina_Adult.Tissue", "RPE_Cell.Line", "RPE_Adult.Tissue",'Cornea_Adult.Tissue','body' ]# add body back in  at some point
tissues=['Retina','RPE','Cornea','body']
sample_names=sample_dict.keys()
loadSRAtk="module load {} && ".format(config['sratoolkit_version'])
loadSalmon= "module load {} && ".format(config['salmon_version'])

STARindex='ref/STARindex'
ref_fasta='ref/gencodeRef.fa'
ref_GTF='ref/gencodeAno.gtf'
ref_GTF_basic='ref/gencodeAno_bsc.gtf'
ref_PA='ref/gencodePA.fa'
fql=config['fastq_path']
stringtie_full_gtf='ref/all_tissues.combined.gtf'

rule all:
    input: genrMATsinput(subtissues_PE,'_PE'), expand('quant_files/{sampleID}/quant.sf',sampleID=sample_names)
    #,'smoothed_filtered_tpms.csv'
'''
****PART 1**** download files and align to genome
-still need to add missing fastq files
-gffread needs indexed fasta
-need to add versioning of tools to yaml
'''
rule downloadGencode:
    output:ref_fasta,ref_GTF_basic,ref_PA
    shell:
        '''

        wget -O ref/gencodeRef.fa.gz {config[refFasta_url]}
        wget -O ref/gencodeAno_bsc.gtf.gz {config[refGTF_basic_url]}
        wget -O ref/gencodePA_tmp.fa.gz {config[refPA_url]}
        gunzip ref/gencodeRef.fa.gz
        gunzip ref/gencodeAno_bsc.gtf.gz
        gunzip ref/gencodePA_tmp.fa.gz
        module load python/3.6
        python3 scripts/filterFasta.py ref/gencodePA_tmp.fa ref/chroms_to_remove ref/gencodePA.fa
        module load samtools
        samtools faidx ref/gencodePA.fa

        '''


rule build_STARindex:
    input: ref_PA, ref_GTF_basic
    output:STARindex
    shell:
        '''
        module load STAR
        mkdir -p ref/STARindex
        STAR --runThreadN 16 --runMode genomeGenerate --genomeDir {output[0]} --genomeFastaFiles {input[0]} --sjdbGTFfile {input[1]} --sjdbOverhang 100

        '''



rule run_STAR_alignment:
    input: lambda wildcards: [fql+'fastq_files/{}_1.fastq.gz'.format(wildcards.id),fql+'fastq_files/{}_2.fastq.gz'.format(wildcards.id)] if sample_dict[wildcards.id]['paired'] else fql+'fastq_files/{}.fastq.gz'.format(wildcards.id),STARindex
    output:temp('STARbams/{id}/raw.Aligned.out.bam')
    run:
        id=wildcards.id
        sp.run('mkdir -p STARbams/{}'.format(id),shell=True)
        STARcmd_pref= 'module load STAR && STAR --runThreadN 8 --genomeDir ref/STARindex --outSAMstrandField intronMotif'
        STARcmd_suf=' --readFilesCommand gunzip -c --outFileNamePrefix STARbams/{}/raw. --outSAMtype BAM Unsorted '.format(id)
        # might have to add that back later
        if sample_dict[id]['paired']:
            paired=' --readFilesIn {} {} '.format(input[0],input[1])
        else:
            paired=' --readFilesIn {}'.format(input[0])
        sp.run(STARcmd_pref+paired+STARcmd_suf,shell=True)
        #samtools_cmd= 'module load samtools && samtools sort -o Aligned.out.bam --threads 7 raw.Aligned.out.bam

rule sort_bams:
    input:'STARbams/{id}/raw.Aligned.out.bam'
    output:'STARbams/{id}/Aligned.out.bam'
    shell:
        '''
        module load samtools
        samtools sort -o {output[0]} --threads 7 {input[0]}
        '''
'''
****PART 2**** Align with STAR and build Transcriptome
-Reminder that STAR makes the bam even if the alignment fails
-following CHESS paper - run each bam individually through stringtie, then merge to a tissue level, then merge into 1
 use gffcompare at each meerge step;
-12/13/18
    - tried GFFcompare on all samples first gave 300k tx's but salmon couldn't map to them, so will now use
      stringtie merge on a per tissue level, which will cut out a lot of transcripts, then merge with gffcompare.
    - moved gffread > tx into its own rule

'''


rule run_stringtie:
    input: 'STARbams/{sample}/Aligned.out.bam'
    output:'st_out/{sample}.gtf'
    shell:
        '''
        module load stringtie
        stringtie {input[0]} -o {output[0]} -p 8 -G ref/gencodeAno_bsc.gtf
        '''

#gffread v0.9.12.Linux_x86_64/

rule merge_gtfs_by_tissue:
    input: lambda wildcards: tissue_to_gtf(wildcards.tissue, sample_dict)
    output: 'ref/{tissue}_st.gtf'
    shell:
        '''
        pattern={wildcards.tissue}
        awk -v pattern="$pattern" '$4==pattern' sampleTableV2.5.tab | cut -f1 | awk '$0="st_out/"$0".gtf"' -  > ref/{wildcards.tissue}_gtf.loc
        module load stringtie
        stringtie --merge -G ref/gencodeAno_bsc.gtf -o {output[0]} {input}
        '''
rule merge_tissue_gtfs:
    input: expand('ref/{tissue}_st.gtf',tissue=tissues)
    output: stringtie_full_gtf
    shell:
        '''
        gffcompare -r ref/gencodeAno_bsc.gtf -o ref/all_tissues {input}
        '''

rule make_tx_fasta:
    input: stringtie_full_gtf
    output: 'ref/combined_stringtie_tx.fa'
    shell:
        '''
        ./gffread/gffread -w {output[0]} -g {ref_PA} {input[0]}
        '''

'''
****PART 4**** rMATS
-the rmats shell script double bracket string thing works even though it looks wrong
-updated STAR cmd to match rmats source
-right now we have little to no body samples that are single ended, so gotta deal wit that< only running paired samples rn
'''
rule rebuild_star_index:
    input: ref_PA, stringtie_full_gtf
    output:'ref/STARindex_stringtie'
    shell:
        '''
        module load STAR
        mkdir -p {output[0]}
        STAR --runThreadN 16 --runMode genomeGenerate --genomeDir {output[0]} --genomeFastaFiles {input[0]} --sjdbGTFfile {input[1]} --sjdbOverhang 100
        '''



rule realign_STAR:
    input: lambda wildcards: [fql+'fastq_files/{}_1.fastq.gz'.format(wildcards.id), fql+'fastq_files/{}_2.fastq.gz'.format(wildcards.id)] if sample_dict[wildcards.id]['paired'] else fql + 'fastq_files/{}.fastq.gz'.format(wildcards.id), 'ref/STARindex_stringtie'
    output:'STARbams_realigned/{id}/Aligned.out.bam'
    run:
        id=wildcards.id
        sp.run('mkdir -p STARbams_realigned/{}'.format(id),shell=True)
        STARcmd_pref= 'module load STAR && STAR  --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --alignSJDBoverhangMin 6 '
        STARcmd_pref+=' --alignIntronMax 299999 --runThreadN 8 --genomeDir {} --sjdbGTFfile {} '.format(input[-1],'ref/combined_final.gtf ')
        STARcmd_suf=' --readFilesCommand gunzip -c --outFileNamePrefix STARbams_realigned/{}/ '.format(id)
        if sample_dict[id]['paired']:
            paired=' --readFilesIn {} {} '.format(input[0],input[1])
        else:
            paired=' --readFilesIn {}'.format(input[0])
        sp.run(STARcmd_pref+paired+STARcmd_suf,shell=True)


rule preprMats_running:# this is going to run ultiple times, but should not be a problem
    input: lambda wildcards: subtissue_to_bam(wildcards.st,'PE',sample_dict)
    output:'ref/{st}.rmats.txt'
    shell:
        '''
        module load R
        Rscript scripts/preprMATS.R {config[sampleFile]}
        '''


rule runrMATS:
    input: 'ref/{tissue1}.rmats.txt','ref/{tissue2}.rmats.txt','ref/STARindex_stringtie',stringtie_full_gtf
             #,'ref/{tissue1}.rmats.txt','ref/{tissue2}.rmats.txt'
    output: 'rmats_out/{tissue1}_VS_{tissue2}'
    # might have to change read length to some sort of function
    shell:
        #**need to fix this for bam mode**
        '''
        wc={wildcards.tissue1}
        id=${{wc: -2}}
        paired=PE
        if [ $id = $paired ]
        then
           flag=paired
        else
            flag=single
        fi
        module load rmats
        rmats --b1 {input[0]} --b2 {input[1]}  -t $flag --readLength 130 --gtf {input[3]} --bi {input[2]} --od {output[0]}
        '''
'''
PART 5 - quantify new transcripts
'''

rule build_salmon_index:
    input:  'ref/combined_stringtie_tx.fa'
    output:'ref/salmonindex_st'
    run:
        salmonindexcommand=loadSalmon + 'salmon index -t {} --gencode -i {} --type quasi --perfectHash -k 31'.format(input[0],output[0])
        sp.run(salmonindexcommand, shell=True)

rule run_salmon:
    input: lambda wildcards: [fql+'fastq_files/{}_1.fastq.gz'.format(wildcards.sampleID),fql+'fastq_files/{}_2.fastq.gz'.format(wildcards.sampleID)] if sample_dict[wildcards.sampleID]['paired'] else fql+'fastq_files/{}.fastq.gz'.format(wildcards.sampleID),
        'ref/salmonindex_st'
    output: 'quant_files/{sampleID}/quant.sf'
    run:
        id=wildcards.sampleID
        #tissue=wildcards.tissue
        paired=sample_dict[id]['paired']
        if paired:
            salmon_command=loadSalmon + 'salmon quant -p 4 -i {} -l A --gcBias --seqBias  -1 {} -2 {} -o {}'.format(input[2],input[0],input[1],'quant_files/{}'.format(id))
        else:
            salmon_command=loadSalmon + 'salmon quant -p 4 -i {} -l A --gcBias --seqBias -r {} -o {}'.format(input[1],input[0],'quant_files/{}'.format(id))
        sp.run(salmon_command,shell=True)
