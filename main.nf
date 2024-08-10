#! /usr/bin/env nextflow

reference_map = [:]
ref_csv = file("./ref.csv") 
allLines = ref_csv.readLines()

for (line : allLines){
    splited_line = line.split(",") 
    reference_map[splited_line[0]] = splited_line[1] 
}


include {dada2_varcall as dada_sample1} from './modules/dada2'
include {dada2_varcall as dada_sample2} from './modules/dada2'
include {dada2_varcall as dada_sample3} from './modules/dada2'

params.input_files = "./fastq_sample/*.fastq"
params.results_dir = "./results/"


file_channel = Channel.fromPath( params.input_files, checkIfExists: true )
                      .map { it -> [it.baseName, it] }

variant_call_setup = multiMapCriteria {
    sample1: [it,['HAC01','HAC08','HAC15'].contains(it[0])]
    sample2: [it,['HAC02','HAC09','HAC15'].contains(it[0])]
    sample3: [it,['HAC03','HAC10','HAC15'].contains(it[0])]
}

process filter_fastq {

    publishDir("${params.results_dir}/after_filtering", mode: 'copy')

    input: 
        tuple val(fastq_name), path(fastq_file)
    output: 
        tuple val(fastq_name),path("*_filtered.fastq")
    script:
        template "filter.py" 
} 

process BWA { 

    publishDir("${params.results_dir}/bwa", mode: 'copy')

    input: 
        tuple val(fastq_name), path(filtered_fastq)
    output:
        tuple val(fastq_name), path("*_aligned.bam")
    script:
        ref_name = reference_map[fastq_name] + '.fa'
        out_name = fastq_name +'_aligned.bam'
        """
        bwa mem $baseDir/references/$ref_name $filtered_fastq | samtools sort -o $out_name
        
        """
}
process bam_index {
    
    publishDir("${params.results_dir}/bwa", mode: 'copy')
    
    input:
        tuple val(fastq_name), path(aligned_bam)
    output:
        tuple val(fastq_name), path("*_aligned.bam.bai")
    script: 
 
    """
    samtools index $aligned_bam 
    """
}

process qc_nano {
    
    publishDir("${params.results_dir}/QC", mode: 'copy')

    input: 
        tuple val(fastq_name), path(filtered_fastq) 
    output: 
        path('*_qc' )
    script:
        fastq_dir_name = fastq_name + '_qc' 
        
        """
        NanoPlot --fastq $filtered_fastq -o $fastq_dir_name --prefix $fastq_name
        
        """
}

process amplicon_clipping { 
    input:
        tuple val(fastq_name), path(to_clip) 
    output: 
        tuple val(fastq_name), path("*_clipped.bam")
    script: 
        out_name = fastq_name + '_clipped.bam'
    """
    samtools ampliconclip -b $baseDir/HVR+RS.bed --hard-clip --clipped --both-ends $to_clip -o $out_name
    """
}
process multi_qc { 
    publishDir("${params.results_dir}/multi_QC", mode: 'copy')
    input: 
        path(qc_dir)
    output: 
        path("multiqc_report.html")
        path("multiqc_data") 
    script: 
    
    """
    multiqc .
    """

}

process bam_foward_reads { 
    publishDir("${params.results_dir}/fastq/fov", mode: 'copy')
    input:
        tuple val(fastq_name), path(bam_file) 
    output: 
        tuple val(fastq_name), path('*_fov.fq')
    script: 
        out_name = fastq_name +'_fov.fq'
        """
        samtools view $bam_file -h -F 20 -o $out_name
        """
}

process bam_rev_reads { 
    publishDir("${params.results_dir}/fastq/rev", mode: 'copy')
    input:
        tuple val(fastq_name), path(bam_file) 
    output: 
        tuple val(fastq_name), path('*_rev.fq')
    script: 
        out_name = fastq_name +'_rev.fq'
        """
        samtools view $bam_file -h -f 16 -o $out_name
        """
}

process rev_compl_fastq{ 
    publishDir("${params.results_dir}/fastq/fov", mode: 'copy')
    input:
        tuple val(fastq_name), path(rev_file) 
    output:
        tuple val(fastq_name), path('*_rev_fov.fq') 
    script:
        template "rev_compl.py" 

}

process merge_foward_only { 
    publishDir("${params.results_dir}/fastq/fov_only", mode: 'copy')
    input: 
        tuple val(fastq_name), path(to_merged)
    output: 
        tuple val(fastq_name), path('*_fow_only.fq') 
    script:
        out_name = fastq_name + '_fow_only.fq'
        """
        cat ${to_merged[0]} ${to_merged[1]} > $out_name
        """

}


workflow {
    main: 
    filter_fastq(file_channel)|BWA|(bam_rev_reads&bam_foward_reads)
    rev_compl_fastq(bam_rev_reads.out)
    bam_foward_reads.out.mix(rev_compl_fastq.out).groupTuple()|merge_foward_only
    merge_foward_only.out.multiMap(variant_call_setup).set{all_samples}
    dada_sample3(all_samples.sample3.collect(flat:false))
    dada_sample2(all_samples.sample2.collect(flat:false))
    dada_sample1(all_samples.sample1.collect(flat:false))
    bam_index(BWA.out)
    qc_nano(filter_fastq.out)
    multi_qc(qc_nano.out.collect())
  

}
