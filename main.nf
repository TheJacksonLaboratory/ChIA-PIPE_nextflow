#!/usr/bin/env nextflow
// Author : Sai Lek
// Email  : Sai.Lek@jax.org 
// Date   : 10/10/2021

import Helpers
import Logos

logo = new Logo()
println logo.show()


// ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Parameter Defaults ~ ~ ~ ~ ~ ~
def setParamDefaults() {
    params.help = false

    // Configurable variable parameters specific to individual runs:
    params.fastqInputs = null // paths to input fastqs. Comma delimited (e.g., /path/to/R1,/path/to/R2).
    params.outputDir   = null // Base directory where output will be stored.
    params.bwaIndex    = null // absolute path to bwa index.
    params.run         = null
    params.runType     = null
    params.ipFactor    = null
    params.cellType    = null
    params.inputCtrl   = null
    params.blackList   = null
    params.peakCaller  = null
    params.basicFolder = null
    params.tmpdir      = null

    // NOTE: See nextflow.config for ChIA-PIPESupportFiles actual default!
    params.CHIAPETSupportFiles = ''  
}


setParamDefaults()


def helpMessage() {
    log.info"""
    =========================================
      ${workflow.manifest.name} v${workflow.manifest.version}
    =========================================
    ${workflow.manifest.description}

    Usage:
      The typical command for running the pipeline is as follows:
        nextflow -c path/to/params.cfg run ${workflow.projectDir}/${workflow.manifest.name} -profile sumner
            (The params.cfg file needs to have the following mandatory parameters
             OR they need to specified on the command line.)

    Mandatory:
        --fastqInputs           Paths to input fastqs. Comma delimited (e.g., /path/to/R1,/path/to/R2)
        --outputDir             Base directory where output will be stored.
        --bwaIndex              Absolute path to Bowtie2 index file

    Optional:
        --CHIAPETSupportFiles   CHIAPET  tool container images
        --tmpdir                Path to hold temporary files for software execution.
        --email                 The email address to send the pipeline report.
        -name                   Name for the pipeline run. If not specified Nextflow will
                                    automatically generate a random mnemonic.

        -profile                Environment config to use.
                                    [ choices: standard (local), sumner ]
    """.stripIndent()
}


// Show help message
if (params.help){
    helpMessage()
    exit 0
}


// ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Param File and Format Checking ~ ~ ~ ~ ~ ~
// Required parameters
if ( ! params.fastqInputs ) {
    exit 1, "Parameter ERROR: fastqInputs ($params.fastqInputs) must be a comma separated list of two _R1 and _R2 fastq filenames."
}
if ( ! params.outputDir ) {
    exit 1, "Parameter ERROR: Output directory parameter must be specified."
}


// Check bwaIndex, append file extensions to check exists
if ( ! params.bwaIndex ) {
    exit 1, "Parameter ERROR: bwaIndex absolute path must be specified."
} else {
    bwaFile = params.bwaIndex + '.fa.bwt'
    if ( ! file(bwaFile).exists() ) {
        exit 1, "Parameter ERROR: bwa index file for param ($params.bwaIndex) does not exist."
    }
}


// Optional params 
if ( ! file(params.CHIAPETSupportFiles).exists() ) {
    exit 1, "Parameter ERROR: CHIAPETSupportFiles directory ($params.CHIAPETSupportFiles) does not exist."
}



// ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Summary Info ~ ~ ~ ~ ~ ~
// Header info
def summary = [:]
summary['Pipeline']         = workflow.manifest.name
summary['Description']      = workflow.manifest.description
if(workflow.revision) {
    summary['Pipeline Release'] = workflow.revision
}
summary['Run Name']         = workflow.runName
summary['User']             = workflow.userName
summary['Config Profile']   = workflow.profile
summary['Config Files']     = workflow.configFiles
summary['Command Line']     = workflow.commandLine
summary['Nextflow Info']    = "v${nextflow.version}, build: ${nextflow.build}, on ${nextflow.timestamp}"
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Workflow dir']     = workflow.projectDir
if(workflow.containerEngine) {
    summary['Container Engine'] = "$workflow.containerEngine"
}


// Pipeline Params:
summary['Parameters......']  = ''
if(params.email) {
    summary['.  Email']      = params.email
}
summary['.  Output dir']     = params.outputDir
summary['.  FASTQ inputs']   = params.fastqInputs
summary['.  Bwa Index']      = params.bwaIndex
summary['.  CHIAPETSupportFiles'] = params.CHIAPETSupportFiles
summary['.  TMP dir']        = params.tmpdir
summary['Run Start Time']    = workflow.start


// print summary header:
println Summary.show(summary)


// ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Opening Variables and Channels ~ ~ ~ ~ ~
def timestamp = new Date().format("yyyyMMdd'T'hhmmSSS") // `date +"%Y%M%dT%H%M%N"`
def tmpdir    = params.tmpdir


// ~~~~~ Determine sampleID and outputDir ~~~~
def fqPairRE = ~/_R1.*\..*f[ast]*q.*$/
fqin = params.fastqInputs.tokenize(",")
fqR1 = file(fqin[0])

fqR1p = fqR1.toAbsolutePath().toString()
fqR2p = file(fqin[1]).toAbsolutePath().toString()


def sampleID = ( fqR1.name - fqPairRE )


// ~~~~~~~~~ Initial Channel of SampleID and Fastq Data ~~~~
Channel.of( sampleID, fqR1p, fqR2p )
       .toList()
       .set { sample_fastqs_ch }


// ~~~~~~~~~~~~~~~~ Primary publishDir == sample_outdir ~~~~~
def sample_outdir = "${params.outputDir}/${timestamp}_${sampleID}/"


// defined names
pairlabel="singlelinker.paired"
pair_map_qual="30"
pair_suffix="UU"
singlabel="singlelinker.single"
single_map_qual="10"
single_suffix="UxxU"
nonelabel="none"
none_map_qual="30"
none_suffix="UU"
selfbp="8000"
extbp="500"


cis_file="${params.run}.e500.clusters.cis.gz"
be3_file="${params.run}.e500.clusters.cis.BE3"
be2_file="${params.run}.e500.clusters.cis.BE2"


hic_file="ChIA-PET_${params.genome}_${params.cellType}_${params.ipFactor}_${params.run}_${params.runType}_pairs.hic"
url="http://ctencode01.jax.org/chiapet/dev/${hic_file}"


out_file="${params.run}.final_stats.tsv"


// ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  Processes ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
/*
Detect Linker & pigziiping	(cpu stag & pigz)
Mapping paired tags		(cpu memaln, pair, span, dedup, dedup span, cluster)
Bam2pair & juicer		(bam2pairs & juicer)
Mapping single tag		(cpu memaln, pair, span, dedup, dedup span)
Mapping none tag		(cpu memaln, pair, span, dedup, dedup span)
Converting file formats		(samtools)
Make bedgraph			(cpu bedtools)
Make bigwig			(kentUtils)
Peak Calling			(spp)
Peak Calling			(macs2)
Final stats			(shell cmds)
*/


if (params.runType == 'miseq' || params.runType == 'novaseqsp' || params.runType == 'qcseq'){
   resource = 'level_1'
   thread   = 20
}
else if (params.runType == 'hiseq' || params.runType == 'nextseq' || params.runType == 'novaseq'){
   resource = 'level_2'
   thread   = 30
}
else {
   resource = 'level_3'
   thread   = 40
} 


// ~~~~~~~~~~~~~~~~~~~~~~~ Linker Detection ~~~~~~~~~~~~~~~~~~~~~~
process detect_linker_pigz {
    tag "$sampleID"
    label 'cpuprg'
    label "$resource"

    publishDir "${sample_outdir}", pattern: "*.cpu", mode: 'copy'
    publishDir "${sample_outdir}", pattern: "*.stat", mode: 'copy'

    input:
    tuple sampleID, fqR1, fqR2 from sample_fastqs_ch

    output:
    tuple sampleID, file("${sampleID}.${pairlabel}.fastq*") \
          into ( linker_pe_fq )
    tuple sampleID, file("${sampleID}.${singlabel}.fastq*") \
          into ( linker_se_fq )
    tuple sampleID, file("${sampleID}.none.fastq*") into none_fq 
   
    tuple sampleID, file("*.stat") into stat_cpuprg 
    tuple sampleID, file("*.1.log") into log_cpuprg

    script:
    log.info "--- ChIA-PIPE start ---"
    log.info "-----Linker Detection & pigz on ${sampleID} -----"
    """
    echo "--- ChIA-PIPE start ---" >> ${params.run}.1.log
    echo "--- linker detection ---" >> ${params.run}.1.log
    echo "`date`" >> ${params.run}.1.log

    cpu stag -W -T 18 \
        -t $thread \
        -O ${params.run} \
        ${fqR1} ${fqR2} \
        >> ${params.run}.1.log \
        2>&1
     
    echo "--- linker detection completed ---" >>${params.run}.1.log
    echo "`date`" >> ${params.run}.1.log

    cpu stat -s -p -T 18 \
        -t $thread \
        ${params.run}.cpu \
        2>> ${params.run}.1.log \
        1>${params.run}.stat
     
    echo echo "--- statistics done ---"  >>${params.run}.1.log
    echo "`date`" >> ${params.run}.1.log
    echo "--- pigziiping ---" >> ${params.run}.1.log

    pigz -p \
        $thread \
        ${params.run}.singlelinker.paired.fastq \
        >> ${params.run}.1.log \
        2>&1

    pigz -p \
        $thread \
        ${params.run}.singlelinker.single.fastq \
        >> ${params.run}.1.log \
        2>&1

    pigz -p \
        $thread \
        ${params.run}.none.fastq \
        >> ${params.run}.1.log \
        2>&1

    pigz -p \
        $thread \
        ${params.run}.conflict.fastq \
        >> ${params.run}.1.log \
        2>&1

    pigz -p \
        $thread \
        ${params.run}.tied.fastq \
        >> ${params.run}.1.log \
        2>&1

    echo echo "--- pigziiping done ---"  >>${params.run}.1.log
    echo "`date`" >> ${params.run}.1.log

    echo "--- Mapping Start ---" >> ${params.run}.1.log
    echo "`date`" >> ${params.run}.1.log

    echo "START  ${sampleID} cpu memaln .." >> ${params.run}.1.log
    echo  "Mapping paired tags .." >> ${params.run}.1.log

    """
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~ Mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
process map_pair {
    tag "$sampleID"
    label 'cpuprg'
    label "$resource"

    publishDir "${sample_outdir}", pattern: "*.span.xls", mode: 'copy' 
    publishDir "${sample_outdir}", pattern: "*.dedup.lc", mode: 'copy' 
    publishDir "${sample_outdir}", pattern: "*.juice.gz", mode: 'copy' 
    publishDir "${sample_outdir}", pattern: "*.stat.xls", mode: 'copy' 
    publishDir "${sample_outdir}", pattern: "*.trans.gz", mode: 'copy' 
    publishDir "${sample_outdir}", pattern: "*.BE*", mode: 'copy' 

    input:
    tuple sampleID, file(fqpe) from linker_pe_fq
    tuple sampleID, file(logtxt) from log_cpuprg
    
    output:
    tuple sampleID, file("${sampleID}.${pairlabel}.${pair_suffix}.nr.bam") \
          into ( map_pe_bam, map_pe_bam1 )
    tuple sampleID, file("${sampleID}.${pairlabel}.${pair_suffix}.span.xls") \
          into xls_pe_bam
    tuple sampleID, file("${sampleID}.${pairlabel}.${pair_suffix}.nr.span.xls") \
          into xls_nrpe_bam
    tuple sampleID, file("*.clusters.*.gz") \
          into cluster_gz
    tuple sampleID, file("*.log") into log_map
    tuple sampleID, file("*.dedup.lc") into log_map_dedup
    tuple sampleID, file("*.gz") into log_map_gz
    tuple sampleID, file("*.stat.xls") into log_map_stat
    tuple sampleID, file("*.BE*") into log_map_BE
    tuple sampleID, file("*.2.log") into log_map_pe

    script:
    log.info "----- Mapping on ${sampleID} -----"
    """
    cp ${logtxt} ${params.run}.2.log
    echo "--- Mapping Starts ---" >> ${params.run}.2.log
    echo "`date`" >> ${params.run}.2.log
  
    echo "START  ${sampleID} cpu memaln .." >> ${params.run}.2.log
    echo  "Mapping paired tags .." >> ${params.run}.2.log

    cpu memaln -T ${pair_map_qual} \
        -t $thread \
        ${params.bwaIndex}.fa \
        ${params.run}.${pairlabel}.fastq.gz 1>${params.run}.${pairlabel}.sam \
        2>> ${params.run}.2.log

    pigz -p \
        $thread \
        ${params.run}.${pairlabel}.sam \
        >> ${params.run}.2.log \
        2>&1

    echo "ENDED pair mapping" >> ${params.run}.2.log

    echo  "STARTED ${sampleID} cpu pair .. >> ${params.run}.2.log"
    echo  "Pairing paired tags .. >> ${params.run}.2.log"

    cpu pair -S -q 30 -s ${selfbp} \
        -t $thread \
        ${params.run}.${pairlabel}.sam.gz \
        1>${params.run}.${pairlabel}.stat.xls \
        2>> ${params.run}.2.log

    echo  "ENDED ${sampleID} cpu pair .. >> ${params.run}.2.log"
    
    echo  "STARTED ${sampleID} cpu span .. >> ${params.run}.2.log"
    echo  "Computing span of paired tags .. >> ${params.run}.2.log"

    cpu span -g \
        -t $thread -s ${selfbp} \
        ${params.run}.${pairlabel}.${pair_suffix}.bam \
        1>${params.run}.${pairlabel}.${pair_suffix}.span.xls \
        2>> ${params.run}.2.log

    echo  "ENDED ${sampleID} span pair .. >> ${params.run}.2.log"
    
    echo  "STARTED ${sampleID} cpu dedup .. >> ${params.run}.2.log"
    echo  "De-duplicating paired tags UU .. >> ${params.run}.2.log"

    cpu dedup -g \
        -t $thread -s ${selfbp} \
        ${params.run}.${pairlabel}.${pair_suffix}.bam \
        1>${params.run}.${pairlabel}.${pair_suffix}.dedup.lc \
        2>> ${params.run}.2.log

    echo  "ENDED ${sampleID} cpu dedup .. >> ${params.run}.2.log"    

    echo  "STARTED ${sampleID} cpu dedup span.. >> ${params.run}.2.log"
    echo  "Computing span of paired tags UU nr .. >> ${params.run}.2.log"

    cpu span \
        -t $thread -s ${selfbp} \
        ${params.run}.${pairlabel}.${pair_suffix}.nr.bam \
        1>${params.run}.${pairlabel}.${pair_suffix}.nr.span.xls \
        2>> ${params.run}.2.log

    echo  "ENDED ${sampleID} cpu dedup span.. >> ${params.run}.2.log"

    cpu cluster \
        -m -s ${selfbp} \
        -B 1000 -5 5,0 \
        -3 3,$extbp \
        -j -x -v 1 -g \
        -t $thread \
        -O ${params.run}.e$extbp \
        ${params.run}.${pairlabel}.${pair_suffix}.nr.bam \
        1>${params.run}.2.log \
        2>> ${params.run}.2.log

    echo  "ENDED ${sampleID} $pairlabel cpu clustering.. >> ${params.run}.2.log"

    mv ${params.run}.e500.clusters.cis.chiasig.gz ${params.run}.e500.clusters.cis.gz
    mv ${params.run}.e500.clusters.trans.chiasig.gz ${params.run}.e500.clusters.trans.gz

    zcat ${cis_file} | awk '{ if ( \$7 >= 3 ) print }' > ${be3_file}
    zcat ${cis_file} | awk '{ if ( \$7 >= 2 ) print }' > ${be2_file}

    """
}

process bam2pairs_juicer {
    tag "$sampleID"
    label 'tools'
    label "$resource"
    
    publishDir "${sample_outdir}", pattern: "*qcseq_pairs.hic", mode: 'copy'
    publishDir "${sample_outdir}", pattern: "*.px2", mode: 'copy'
    publishDir "${sample_outdir}", pattern: "*.pairs.gz", mode: 'copy'

    input:
    tuple sampleID, file(pebam) from map_pe_bam1
    tuple sampleID, file(logtxt) from log_map_pe

    output:
    tuple sampleID, file("*${sampleID}*qcseq_pairs.hic") \
          into ( juicer_hic )

    tuple sampleID, file("*.3.log") into log_juicer
    tuple sampleID, file("*qcseq_pairs.hic") into hic_juicer
    tuple sampleID, file("*.pairs.gz") into log_map_pairgz
    tuple sampleID, file("*.px2") into log_map_px2

    script:
    log.info "----- Running bam2pairs & Juicer  on ${sampleID} -----"
    """
    bam2pairs -c ${params.bwaIndex}.chrom.sizes ${pebam} ${params.run}
    cp ${logtxt} ${params.run}.3.log
    echo "--- Juicer Starts ---" >> ${params.run}.3.log
    echo "`date`" >> ${params.run}.3.log

    java -Xmx16g -jar /usr/local/bin/juicer_tools.jar pre -r \
        2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000 \
        ${params.run}.bsorted.pairs.gz ${hic_file} ${params.bwaIndex}.chrom.sizes >> ${params.run}.3.log
 
    """
} 

process map_single {
    tag "$sampleID"
    label 'cpuprg'
    label "$resource"

    publishDir "${sample_outdir}", pattern: "*.span.xls", mode: 'copy' 
    publishDir "${sample_outdir}", pattern: "*.stat.xls", mode: 'copy' 
    publishDir "${sample_outdir}", pattern: "*.dedup.lc", mode: 'copy' 

    input:
    tuple sampleID, file(fqse) from linker_se_fq
    tuple sampleID, file(logtxt) from log_juicer

    output:
    tuple sampleID, file("${sampleID}.${singlabel}.sam*") \
          into ( map_se_sam )
    tuple sampleID, file("${sampleID}.${singlabel}.${single_suffix}.nr.bam") \
          into ( map_se_bam )

    tuple sampleID, file("*.4.log") into log_map_se
    tuple sampleID, file("*.xls") into log_map_se_xls
    tuple sampleID, file("*.dedup.lc") into log_map_se_dedup

    script:
    log.info "-----Mapping Single on ${sampleID} -----"
    """
    cp ${logtxt} ${params.run}.4.log
    echo "--- Mapping Single Starts ---" >> ${params.run}.4.log
    echo "`date`" >> ${params.run}.4.log

    echo "START  ${sampleID} cpu memaln .." >> ${params.run}.4.log
    echo  "Mapping single tag .." >> ${params.run}.4.log

    cpu memaln -T ${single_map_qual} \
        -t $thread \
        ${params.bwaIndex}.fa \
        ${params.run}.${singlabel}.fastq.gz 1>${params.run}.${singlabel}.sam \
        2>> ${params.run}.4.log

    pigz -p \
        $thread \
        ${params.run}.${singlabel}.sam \
        >> ${params.run}.4.log \
        2>&1

    echo "ENDED single mapping" >> ${params.run}.4.log

    cpu pair -S -q ${single_map_qual} -s ${selfbp} \
        -t $thread \
        ${params.run}.${singlabel}.sam.gz \
        1>${params.run}.${singlabel}.stat.xls \
        2>> ${params.run}.4.log

    echo  "ENDED ${sampleID} cpu pair .. >> ${params.run}.4.log"

    echo  "STARTED ${sampleID} cpu span .. >> ${params.run}.4.log"
    echo  "Computing span of single tag .. >> ${params.run}.4.log"

    cpu span -g \
        -t $thread -s ${selfbp} \
        ${params.run}.${singlabel}.${single_suffix}.bam \
        1>${params.run}.${singlabel}.${single_suffix}.span.xls \
        2>> ${params.run}.4.log

    echo  "ENDED ${sampleID} span single .. >> ${params.run}.4.log"

    echo  "STARTED ${sampleID} $singlabel cpu dedup .. >> ${params.run}.4.log"
    echo  "De-duplicating single tag $singlabel .. >> ${params.run}.4.log"

    cpu dedup -g \
        -t $thread -s ${selfbp} \
        ${params.run}.${singlabel}.${single_suffix}.bam \
        1>${params.run}.${singlabel}.${single_suffix}.dedup.lc \
        2>> ${params.run}.4.log

    echo  "ENDED ${sampleID} $singlabel cpu dedup .. >> ${params.run}.4.log"

    echo  "STARTED ${sampleID} cpu dedup span.. >> ${params.run}.4.log"
    echo  "Computing span of single tag $singlabel nr .. >> ${params.run}.4.log"

    cpu span \
        -t $thread -s ${selfbp} \
        ${params.run}.${singlabel}.${single_suffix}.nr.bam \
        1>${params.run}.${singlabel}.${single_suffix}.nr.span.xls \
        2>> ${params.run}.4.log

    echo  "ENDED ${sampleID} $singlabel cpu dedup span.. >> ${params.run}.4.log"

    """
}

process map_none {
    tag "$sampleID"
    label 'cpuprg'
    label "$resource"

    publishDir "${sample_outdir}", pattern: "*.log", mode: 'copy'
    publishDir "${sample_outdir}", pattern: "*.span.xls", mode: 'copy' 
    publishDir "${sample_outdir}", pattern: "*.stat.xls", mode: 'copy' 
    publishDir "${sample_outdir}", pattern: "*.dedup.lc", mode: 'copy' 

    input:
    tuple sampleID, file(fqno) from none_fq 
    tuple sampleID, file(logtxt) from log_map_se

    output:
    tuple sampleID, file("${sampleID}.none.sam*") \
          into ( map_no_sam )
    tuple sampleID, file("${sampleID}.${nonelabel}.${none_suffix}.nr.bam") \
          into ( map_no_bam )

    tuple sampleID, file("*.5.log") into log_map_no
    tuple sampleID, file("*.xls") into log_map_no_xls
    tuple sampleID, file("*.dedup.lc") into log_map_no_dedup

    script:
    log.info "-----Mapping None on ${sampleID} -----"
    """
    cp ${logtxt} ${params.run}.5.log
    mv ${logtxt} ${logtxt}.old
    mv ${logtxt}.old ${params.run}.5.log
    rm -rf ${logtxt}.old 
    echo "--- Mapping None Starts ---" >> ${params.run}.5.log
    echo "`date`" >> ${params.run}.5.log

    echo "START  ${sampleID} cpu memaln .." >> ${params.run}.5.log
    echo  "Mapping None tag .." >> ${params.run}.5.log

    cpu memaln -T ${none_map_qual} \
        -t $thread \
        ${params.bwaIndex}.fa \
        ${params.run}.${nonelabel}.fastq.gz 1>${params.run}.${nonelabel}.sam \
        2>> ${params.run}.5.log

    pigz -p \
        $thread \
        ${params.run}.${nonelabel}.sam \
        >> ${params.run}.5.log \
        2>&1

    echo "ENDED single mapping" >> ${params.run}.5.log
    echo  "STARTED Pairing ${sampleID} $nonelabel tag .. >> ${params.run}.5.log"

    cpu pair -S -q ${none_map_qual} -s ${selfbp} \
        -t $thread \
        ${params.run}.${nonelabel}.sam.gz \
        1>${params.run}.${nonelabel}.stat.xls \
        2>> ${params.run}.5.log

    echo  "ENDED ${sampleID} $nonelabel cpu pair .. >> ${params.run}.5.log"

    echo  "STARTED ${sampleID} $nonelabel cpu span .. >> ${params.run}.5.log"
    echo  "Computing span of None tag .. >> ${params.run}.5.log"

    cpu span -g \
        -t $thread -s ${selfbp} \
        ${params.run}.${nonelabel}.${none_suffix}.bam \
        1>${params.run}.${nonelabel}.${none_suffix}.span.xls \
        2>> ${params.run}.5.log

    echo  "ENDED ${sampleID} $nonelabel span single .. >> ${params.run}.5.log"

    echo  "STARTED ${sampleID} $nonelabel cpu dedup .. >> ${params.run}.5.log"
    echo  "De-duplicating None tag $nonelabel .. >> ${params.run}.5.log"

    cpu dedup -g \
        -t $thread -s ${selfbp} \
        ${params.run}.${nonelabel}.${none_suffix}.bam \
        1>${params.run}.${nonelabel}.${none_suffix}.dedup.lc \
        2>> ${params.run}.5.log

    echo  "ENDED ${sampleID} $nonelabel cpu dedup .. >> ${params.run}.5.log"

    echo  "STARTED ${sampleID} $nonelabel cpu dedup span.. >> ${params.run}.5.log"
    echo  "Computing span of None tag $nonelabel nr .. >> ${params.run}.5.log"

    cpu span \
        -t $thread -s ${selfbp} \
        ${params.run}.${nonelabel}.${none_suffix}.nr.bam \
        1>${params.run}.${nonelabel}.${none_suffix}.nr.span.xls \
        2>> ${params.run}.5.log

    echo  "ENDED ${sampleID} $nonelabel cpu dedup span.. >> ${params.run}.5.log"

    """
}

process convert_format {
    tag "$sampleID"
    label 'samtools'
    label "$resource"

    publishDir "${sample_outdir}", pattern: "*.log", mode: 'copy'
    publishDir "${sample_outdir}", pattern: "*.BROWSER.bam*", mode: 'copy'

    input:
    tuple sampleID, file(pebam) from map_pe_bam
    tuple sampleID, file(sebam) from map_se_bam
    tuple sampleID, file(nobam) from map_no_bam

    output:
    tuple sampleID, file("${sampleID}.for.BROWSER.bam") \
          into ( cvt_fmt_bam, cvt_fmt_bam1 )
    tuple sampleID, file("${sampleID}.for.BROWSER.bam.bai") \
          into ( cvt_fmt_bai, cvt_fmt_bai1 )
    file("nr.tag.txt") into cvt_fmt_tag
    tuple sampleID, file("*.BROWSER.bam*") into log_cvt_fmt_bam

    script:
    log.info "-----Converting file formats on ${sampleID} -----"
    """
    samtools view \
        -F 2048 -h ${pebam} \
        | awk 'length(\$10) > 30 || \$1 ~ /^@/' \
        | samtools sort -@ 16 - -o ${params.run}.${pairlabel}.${pair_suffix}.nr.sorted.bam

    samtools view \
        -F 2048 -h ${sebam} \
        | awk 'length(\$10) > 30 || \$1 ~ /^@/' \
        | samtools sort -@ 16 - -o ${params.run}.${singlabel}.${single_suffix}.nr.sorted.bam

    samtools view \
        -F 2048 -h ${nobam} \
        | awk 'length(\$10) > 30 || \$1 ~ /^@/' \
        | samtools sort -@ 16 - -o ${params.run}.${nonelabel}.${none_suffix}.nr.sorted.bam

    samtools merge ${sampleID}.for.BROWSER.bam \
        ${params.run}.${pairlabel}.${pair_suffix}.nr.sorted.bam \
        ${params.run}.${singlabel}.${single_suffix}.nr.sorted.bam \
        ${params.run}.${nonelabel}.${none_suffix}.nr.sorted.bam

    samtools index ${sampleID}.for.BROWSER.bam ${sampleID}.for.BROWSER.bam.bai

    samtools view \
        -c ${sampleID}.for.BROWSER.bam \
        | numfmt --g >nr.tag.txt

    """
}

process make_bedgraph {
    tag "$sampleID"
    label 'cpuprg'
    label "$resource"

    publishDir "${sample_outdir}", pattern: "*.log", mode: 'copy' 
    publishDir "${sample_outdir}", pattern: "*.sorted.bedgraph", mode: 'copy' 

    input:
    tuple sampleID, file(brbam) from cvt_fmt_bam
    tuple sampleID, file(brbai) from cvt_fmt_bai 
    tuple sampleID, file(logtxt) from log_map_no

    output:
    tuple sampleID, file("${sampleID}.for.BROWSER.sorted*") \
          into ( mk_bg_bed )
    tuple sampleID, file("*.6.log") into (log_mk_bg1, log_mk_bg2)

    script:
    log.info "-----Making bedgraph on ${sampleID} -----"
    if(params.blackList != "none")
      """
      cp ${logtxt} ${params.run}.6.log
      bedtools genomecov \
         -ibam ${brbam} \
         -bg > ${params.run}.for.BROWSER.bedgraph

      echo "--- Removing blacklist ---" >> ${params.run}.6.log
      echo "`date`" >> ${params.run}.6.log
      
      mv ${params.run}.for.BROWSER.bedgraph ${params.run}.for.BROWSER.orig.bedgraph

      bedtools subtract \
         -a ${params.run}.for.BROWSER.orig.bedgraph \
         -b ${params.blackList} > ${params.run}.for.BROWSER.bedgraph

      bedSort \
         ${params.run}.for.BROWSER.bedgraph \
         ${params.run}.for.BROWSER.sorted.bedgraph
      """
     else
      """
      cp ${logtxt} ${params.run}.6.log
      bedtools genomecov \
         -ibam ${brbam} \
         -bg > ${params.run}.for.BROWSER.bedgraph

      bedSort \
         ${params.run}.for.BROWSER.bedgraph \
         ${params.run}.for.BROWSER.sorted.bedgraph
      """
}

process make_bigwig {
    tag "$sampleID"
    label 'kentUtils'
    label "$resource"

    publishDir "${sample_outdir}", pattern: "*.log", mode: 'copy'
    publishDir "${sample_outdir}", pattern: "*.bigwig", mode: 'copy' 

    input:
    tuple sampleID, file(srbam) from mk_bg_bed
    tuple sampleID, file(brbai) from cvt_fmt_bai

    output:
    tuple sampleID, file("${sampleID}.for.BROWSER.bigwig*") \
          into ( bg2bw_bw )

    script:
    log.info "-----Making bigwig on ${sampleID} -----"
    """
    bedGraphToBigWig \
       ${srbam} \
       ${params.bwaIndex}.chrom.sizes \
       ${params.run}.for.BROWSER.bigwig 
    
    """
}

process call_peak_spp {
    tag "$sampleID"
    label 'spp'
    label "$resource"

    publishDir "${sample_outdir}", pattern: "*.log", mode: 'copy'
    publishDir "${sample_outdir}", pattern: "*.broadPeak", mode: 'copy' 

    input:
    tuple sampleID, file(brbam) from cvt_fmt_bam1
    tuple sampleID, file(brbai) from cvt_fmt_bai1 
    tuple sampleID, file(logtxt) from log_mk_bg1

    output:
    tuple sampleID, file("${params.run}.log") into log_spp
    tuple sampleID, file("${sampleID}.for.BROWSER*broadPeak") \
          into ( spp_peak )

    when: 'spp' in params.peakCaller || 'SPP' in params.peakCaller

    script:
    log.info "-----Calling peaks on ${sampleID} -----"
    """
    cp ${logtxt} ${params.run}.log
    Rscript ${projectDir}/bin/spp_sumner.R \
        ${brbam} \
        ${params.inputCtrl} \
        . 6 >> ${params.run}.log
    """
}

process call_peak_macs2 {
    tag "$sampleID"
    label 'macs2'
    label "$resource"

    publishDir "${sample_outdir}", pattern: "*.log", mode: 'copy' 
    publishDir "${sample_outdir}", pattern: "*.narrowPeak", mode: 'copy' 

    input:
    tuple sampleID, file(brbam) from cvt_fmt_bam1
    tuple sampleID, file(brbai) from cvt_fmt_bai1 
    tuple sampleID, file(logtxt) from log_mk_bg2

    output:
    tuple sampleID, file("*.log") into log_macs2
    tuple sampleID, file("${sampleID}.no_input_all*") \
          into ( macs2_peak )

    when: 'macs2' in params.peakCaller

    script:
    log.info "-----Calling peaks on ${sampleID} -----"
    """
    cp ${logtxt} ${params.run}.log
    macs2 callpeak \
           --keep-dup all \
           --nomodel -t \
           ${brbam} \
           -f BAM -g hs -n \
           ${params.run}.no_input_all \
           1>> ${params.run}.log \
           2>> ${params.run}.log
    """
}

input_peak = Channel.empty()
if (params.peakCaller == 'spp'){
    input_peak = spp_peak
}
else {
    input_peak = macs2_peak
}

ch_seq = Channel.of( 2, 3, 4, 5, 6, 7, 8, 9, 10 )

process final_stats {
    tag "$sampleID"
    label "$resource"

    publishDir "${sample_outdir}", pattern: "*.log", mode: 'copy' 
    publishDir "${sample_outdir}", pattern: "*.final_stats.tsv", mode: 'copy'

    input:
    tuple sampleID, file(stat) from stat_cpuprg
    tuple sampleID, file(xls) from xls_pe_bam
    tuple sampleID, file(nrxls) from xls_nrpe_bam
    file(nrtag) from cvt_fmt_tag
    tuple sampleID, file(peak) from input_peak
    tuple sampleID, file(cluster) from cluster_gz
    file seqList from ch_seq.collect()

    output:
    tuple sampleID, file("${out_file}") \
          into ( final_st )

    script:
    log.info "----- Getting final stats on ${sampleID} -----"
    n_read_pair="\$(cat ${stat} | grep 'Total pairs' | awk -F'[ \t]' '{print \$3}')"
    read_pair_link="\$(cat ${stat} | grep 'Linker detected' | awk -F '[ \t]' '{print \$3}')"
    frac_link="\$(echo -e \"${read_pair_link}/${n_read_pair}\" | bc -l | xargs printf \"%.2f\n\")"
    n_read_pair="\$(printf \"%'.f\n\" ${n_read_pair})"
    read_pair_link="\$(printf \"%'.f\n\" ${read_pair_link})"
    one_tag="\$(grep 'Single Linker 1 tag (SL/ls)' ${stat} | cut -f2)"
    two_tag="\$(grep 'Single Linker 2 tags (SL/ls)' ${stat} | cut -f2)"
    one_tag="\$(printf \"%'.f\n\" ${one_tag})"
    two_tag="\$(printf \"%'.f\n\" ${two_tag})"
    unique="\$(cat ${xls} | grep 'Total pairs' | awk -F '[\t]' '{print \$2}')"
    nr="\$(cat ${nrxls} | grep 'Total pairs' | awk -F '[\t]' '{print \$2}')"
    redun="\$(echo \"(${unique} - ${nr})/${unique}\" | bc -l)"
    unique="\$(printf \"%'.f\" ${unique})"
    nr="\$(printf \"%'.f\" ${nr})"
    redun="\$(printf \"%'.2f\" ${redun})"
    nr_tag="\$(cat ${nrtag})"
    nr_tag="\$(printf \"%'.f\" ${nr_tag})"
    n_peak="\$(cat ${peak} | wc -l)"
    n_peak="\$(printf \"%'.f\" ${n_peak})"
    self_lig="\$(cat ${nrxls} | grep \"second/best<0.95\" -A5 | awk -F '[\t]' '{if(NR==4)print \$2}')"
    self_lig="\$(printf \"%'.f\" ${self_lig})"
    intra_chr_pet="\$(cat ${nrxls} | grep \"second/best<0.95\" -A5 | awk -F '[\t]' '{if(NR==5)print \$2}')"
    inter_chr_pet="\$(cat ${nrxls} | grep \"second/best<0.95\" -A5 | awk -F '[\t]' '{if(NR==2)print \$2}')"
    pet_ratio="\$(echo \"${intra_chr_pet}/${inter_chr_pet}\" | bc -l)"
    inter_lig_all="\$(echo \"${intra_chr_pet} + ${inter_chr_pet}\" | bc)"
    inter_lig_all="\$(printf \"%'.f\" ${inter_lig_all})"
    intra_chr_pet="\$(printf \"%'.f\" ${intra_chr_pet})"
    inter_chr_pet="\$(printf \"%'.f\" ${inter_chr_pet})"
    pet_ratio="\$(printf \"%'.2f\" ${pet_ratio})"
    singleton="\$(zcat $cluster | awk '\$7==1{print}' | wc -l)"
    singleton="\$(printf \"%'.f\" ${singleton})"
    intra_singleton="\$(zcat *cis.gz | awk '\$7==1{print}' | wc -l)"
    intra_singleton="\$(printf \"%'.f\" ${intra_singleton})"
    inter_singleton="\$(zcat *trans.gz | awk '\$7==1{print}' | wc -l)"
    inter_singleton="\$(printf \"%'.f\" ${inter_singleton})"
    total_cluster_number="\$(zcat $cluster | awk '\$7 !=1{print}' | wc -l)"
    total_cluster_number="\$(printf \"%'.f\" ${total_cluster_number})"
    intra_cluster="\$(zcat *cis.gz | awk '\$7 >=2 {print}' | wc -l)"
    inter_cluster="\$(zcat *trans.gz | awk '\$7 >=2 {print}'| wc -l)"
    cluster_ratio="\$(echo \"${intra_cluster}/${inter_cluster}\" | bc -l)"
    cluster_ratio="\$(printf \"%'.2f\" ${cluster_ratio})"
    intra_cluster="\$(printf \"%'.f\" ${intra_cluster})"
    inter_cluster="\$(printf \"%'.f\" ${inter_cluster})"
    """
    echo -e "Library_ID\t"${params.run} >> ${out_file}
    echo -e "Library_type\t"${params.runType} >> ${out_file}
    echo -e "Reference_genome\t"${params.genome} >> ${out_file}
    echo -e "Cell_type\t"${params.cellType} >> ${out_file}
    echo -e "Factor\t"${params.ipFactor} >> ${out_file}
    echo -e "Contact-map_URL\t"${url} >> ${out_file}
 
    echo -e "Total_read_pairs\t"${n_read_pair} >> ${out_file}
    echo -e "Read_pairs_with_linker\t"${read_pair_link} >> ${out_file}
    echo -e "Fraction_read_pairs_with_linker\t"${frac_link} >> ${out_file}

    echo -e "One_tag\t"${one_tag} >> ${out_file}
    echo -e "PET\t"${two_tag} >> ${out_file}   

    echo -e "Uniquely_mapped_PET\t"${unique} >> ${out_file}
    echo -e "Non-redundant_PET\t"${nr} >> ${out_file}
    echo -e "Redundancy\t"${redun} >> ${out_file}
    echo -e "Non-redundant_tag\t"${nr_tag} >> ${out_file}
    echo -e "Peak\t"$n_peak >> ${out_file}

    echo -e "Self-ligation_PET\t"${self_lig} >> ${out_file}
    echo -e "Inter-ligation_PET\t"${inter_lig_all} >> ${out_file}
    echo -e "Intra-chr_PET\t"${intra_chr_pet} >> ${out_file}
    echo -e "Inter-chr_PET\t"${inter_chr_pet} >> ${out_file}
    echo -e "ratio_of_intra/inter_PET\t"${pet_ratio} >> ${out_file}

    echo -e "Singleton\t"$singleton >> ${out_file}
    echo -e "Intra-chr_singleton\t"$intra_singleton >> ${out_file}
    echo -e "Inter-chr_singleton\t"$inter_singleton >> ${out_file}

    echo -e "PET_cluster\t"${total_cluster_number} >> ${out_file}
    echo -e "ratio_of_intra/inter_cluster\t"${cluster_ratio} >> ${out_file}
    echo -e "Intra-chr_PET_cluster\t"${intra_cluster} >> ${out_file}

    for content in ${seqList}
    do
        i=\$(cat \${content})
        intra_pets_number=\$(zcat *cis.gz | \
            awk -v cutoff=\${i} '\$7 == cutoff {print}' | wc -l | \
            xargs printf "%'.f")

        echo -e "pets_number_"\${i}"\t"\${intra_pets_number} >> ${out_file}
    done

    echo -e "pets_number>10\t"\$(zcat *cis.gz | \
    awk '\$7 >10 {print}' | wc -l | xargs printf "%'.f") >> ${out_file}

    echo -e "Inter-chr_PET_cluster\t"${inter_cluster} >> ${out_file}

    for content in ${seqList}
    do
        i=\$(cat \${content})
        inter_pets_number=\$(zcat *trans.gz | \
            awk -v cutoff=\${i} '\$7 == cutoff {print}' | wc -l | \
            xargs printf "%'.f")

        echo -e "pets_number_"\${i}"\t"\${inter_pets_number} >> ${out_file}
    done

    echo -e "pets_number>10\t"\$(zcat *trans.gz | \
    awk '\$7 >10 {print}' | wc -l | xargs printf "%'.f") >> ${out_file}
    """
}


