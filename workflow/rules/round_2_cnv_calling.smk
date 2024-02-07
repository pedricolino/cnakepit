include: 'cnvkit.smk'
include: 'purecn.smk'

#--- create 2 long strings with all PON (anti-/)target coverage files -------

# Open the file and read the lines
with open('resources/data/lowpurity_samples.tsv', 'r') as file:
    names = file.readlines()

# Remove newline characters

names = [name.strip() for name in names]

# Format the sample names
targets_str = ['results/cnvkit/general/'+name+'.targetcoverage.cnn' for name in names]
antitargets_str = ['results/cnvkit/general/'+name+'.antitargetcoverage.cnn' for name in names]

# Join the formatted strings into a single string with spaces
long_str_target = ' '.join(targets_str)
long_str_antitarget = ' '.join(antitargets_str)

#--- PON reference track creation and fixing/bias correction with PON -------

rule cnvkit_create_panel_of_normals:
    input:
        fasta=config['reference']['fasta'],
        all_target_cnns = expand('results/cnvkit/general/{sample}.targetcoverage.cnn', sample=samples.index),
        all_antitarget_cnns = expand('results/cnvkit/general/{sample}.antitargetcoverage.cnn', sample=samples.index)
    output:
        pon = 'results/cnvkit'+suffix+'/general/pon.cnn',
    benchmark: 'benchmarks/cnvkit'+suffix+'/general/cnvkit_create_panel_of_normals.txt'
    params:
        extra = '--sample-sex female -c',
        amplicon = '--no-edge' if config['amplicon'] else '',
        target_cnn = long_str_target,
        antitarget_cnn = long_str_antitarget,
    resources:
        mem=lambda wildcards, attempt: '%dG' % (16 * attempt),
    log:
        'logs/cnvkit'+suffix+'/general/cnvkit_create_panel_of_normals/log',
    conda:
        '../envs/cnv_calling.yaml'
    shell:
        'cnvkit.py reference -o {output.pon} -f {input.fasta} {params.extra} {params.amplicon} {params.target_cnn} {params.antitarget_cnn} 2> {log}'

use rule cnvkit_fix as cnvkit_fix_w_pon with:
    input:
        target_coverage = 'results/cnvkit/general/{sample}.targetcoverage.cnn',
        antitarget_coverage = 'results/cnvkit/general/{sample}.antitargetcoverage.cnn',
        reference = 'results/cnvkit'+suffix+'/general/pon.cnn',
    output: 'results/cnvkit'+suffix+'/general/{sample}.cnr'
    benchmark: 'benchmarks/cnvkit'+suffix+'/general/fix/{sample}.txt'
    log: 'logs/cnvkit'+suffix+'/general/fix/{sample}.log',


#--- CNVkit segmentation and downstream rules --------------------------------

## First with CBS

use rule cnvkit_segment_cbs as cnvkit_segment_cbs_w_pon with:
    input:
        copy_ratios = 'results/cnvkit'+suffix+'/general/{sample}.cnr',
        germline_vcf = 'results/mutect2/germline/{sample}_germline.vcf'
    output: 'results/cnvkit'+suffix+'/cbs/{sample}.cns',
    benchmark: 'benchmarks/cnvkit'+suffix+'/cbs_segment/{sample}.txt'
    log: 'logs/cnvkit'+suffix+'/cbs/segment/{sample}.log',

use rule cnvkit_scatter_cbs as cnvkit_scatter_cbs_w_pon with:
    input:
        copy_ratio = 'results/cnvkit'+suffix+'/general/{sample}.cnr',
        segment = 'results/cnvkit'+suffix+'/cbs/{sample}.cns',
        germline_vcf = 'results/mutect2/germline/{sample}_germline.vcf'
    output: 'results/cnvkit'+suffix+'/cbs/{sample}_scatter.cnv.pdf'
    benchmark: 'benchmarks/cnvkit'+suffix+'/cbs/{sample}_scatter.txt'
    log: 'logs/cnvkit'+suffix+'/cbs/scatter/{sample}.log',

use rule cnvkit_diagram_cbs as cnvkit_diagram_cbs_w_pon with:
    input:
        copy_ratio = 'results/cnvkit'+suffix+'/general/{sample}.cnr',
        segment = 'results/cnvkit'+suffix+'/cbs/{sample}.cns',
        germline_vcf = 'results/mutect2/germline/{sample}_germline.vcf'
    output: 'results/cnvkit'+suffix+'/cbs/{sample}_diagram.cnv.pdf'
    benchmark: 'benchmarks/cnvkit'+suffix+'/cbs/{sample}_diagram.txt'
    log: 'logs/cnvkit'+suffix+'/cbs/diagram/{sample}.log',

use rule export_seg_cbs as export_seg_cbs_w_pon with:
    input: cns = 'results/cnvkit'+suffix+'/cbs/{sample}.cns',
    output: 'results/cnvkit'+suffix+'/cbs/{sample}.seg',
    benchmark: 'benchmarks/cnvkit'+suffix+'/cbs/export_seg/{sample}.txt'
    log: 'logs/cnvkit'+suffix+'/cbs/export_seg/{sample}.log',

use rule cnvkit_call_cbs as cnvkit_call_cbs_w_pon with:
    input:
        cns = 'results/cnvkit'+suffix+'/cbs/{sample}.cns',
        germline_vcf = 'results/mutect2/germline/{sample}_germline.vcf'
    output: 'results/cnvkit'+suffix+'/cbs/{sample}.icns',
    benchmark: 'benchmarks/cnvkit'+suffix+'/cbs/{sample}.call.txt'
    log: 'logs/cnvkit'+suffix+'/cbs/call/{sample}.log',


## Same with HMM
        
use rule cnvkit_segment_cbs_w_pon as cnvkit_segment_hmm_w_pon with:
    output: 'results/cnvkit'+suffix+'/hmm/{sample}.cns',
    benchmark: 'benchmarks/cnvkit'+suffix+'/hmm_segment/{sample}.txt'
    params: extra = '-m hmm --drop-low-coverage',
    log: 'logs/cnvkit'+suffix+'/hmm/segment/{sample}.log',

use rule cnvkit_scatter_hmm as cnvkit_scatter_hmm_w_pon with:
    input:
        copy_ratio = 'results/cnvkit'+suffix+'/general/{sample}.cnr',
        segment = 'results/cnvkit'+suffix+'/hmm/{sample}.cns',
        germline_vcf = 'results/mutect2/germline/{sample}_germline.vcf'
    output: 'results/cnvkit'+suffix+'/hmm/{sample}_scatter.cnv.pdf'
    benchmark: 'benchmarks/cnvkit'+suffix+'/hmm/{sample}_scatter.txt'
    log: 'logs/cnvkit'+suffix+'/hmm/scatter/{sample}.log',

use rule cnvkit_diagram_hmm as cnvkit_diagram_hmm_w_pon with:
    input:
        copy_ratio = 'results/cnvkit'+suffix+'/general/{sample}.cnr',
        segment = 'results/cnvkit'+suffix+'/hmm/{sample}.cns',
        germline_vcf = 'results/mutect2/germline/{sample}_germline.vcf'
    output: 'results/cnvkit'+suffix+'/hmm/{sample}_diagram.cnv.pdf'
    benchmark: 'benchmarks/cnvkit'+suffix+'/hmm/{sample}_diagram.txt'
    log: 'logs/cnvkit'+suffix+'/hmm/diagram/{sample}.log',

use rule export_seg_cbs as export_seg_hmm_w_pon with:
    input: cns = 'results/cnvkit'+suffix+'/hmm/{sample}.cns',
    output: 'results/cnvkit'+suffix+'/hmm/{sample}.seg',
    benchmark: 'benchmarks/cnvkit'+suffix+'/hmm/export_seg/{sample}.txt'
    log: 'logs/cnvkit'+suffix+'/hmm/export_seg/{sample}.log',

use rule cnvkit_call_cbs_w_pon as cnvkit_call_hmm_w_pon with:
    input:
        cns = 'results/cnvkit'+suffix+'/hmm/{sample}.cns',
        germline_vcf = 'results/mutect2/germline/{sample}_germline.vcf'
    output: 'results/cnvkit'+suffix+'/hmm/{sample}.icns',
    benchmark: 'benchmarks/cnvkit'+suffix+'/hmm/{sample}.call.txt'
    log: 'logs/cnvkit'+suffix+'/hmm/call/{sample}.log',


#--- PureCN rules ------------------------------------------------------------

use rule purecn_cbs_Hclust as purecn_cbs_Hclust_w_pon with:
    input:
        vcf_filt='results/mutect2/filtered/{sample}_filtered.vcf.gz',
        copy_ratios='results/cnvkit'+suffix+'/general/{sample}.cnr',
        seg='results/cnvkit'+suffix+'/cbs/{sample}.seg',
    output: 'results/purecn'+suffix+'/cbs_Hclust/{sample}/{sample}.csv'
    benchmark: 'benchmarks/purecn'+suffix+'/cbs_Hclust/{sample}.txt'
    log: 'logs/purecn'+suffix+'/cbs_Hclust/{sample}.log',
    params:
        cnvkit_method='cbs',
        purecn_method='Hclust',
        sampleid='{sample}',
        random_nb=randint(1,1000),

use rule purecn_hmm_Hclust as purecn_hmm_Hclust_w_pon with:
    input:
        vcf_filt='results/mutect2/filtered/{sample}_filtered.vcf.gz',
        copy_ratios='results/cnvkit'+suffix+'/general/{sample}.cnr',
        seg='results/cnvkit'+suffix+'/hmm/{sample}.seg',
    output: 'results/purecn'+suffix+'/hmm_Hclust/{sample}/{sample}.csv'
    benchmark: 'benchmarks/purecn'+suffix+'/hmm_Hclust/{sample}.txt'
    log: 'logs/purecn'+suffix+'/hmm_Hclust/{sample}.log',
    params:
        cnvkit_method='hmm',
        purecn_method='Hclust',
        sampleid='{sample}',
        random_nb=randint(1,1000),

use rule purecn_cbs_PSCBS as purecn_cbs_PSCBS_w_pon with:
    input:
        vcf_filt='results/mutect2/filtered/{sample}_filtered.vcf.gz',
        copy_ratios='results/cnvkit'+suffix+'/general/{sample}.cnr',
        seg='results/cnvkit'+suffix+'/cbs/{sample}.seg',
        install='results/purecn/PSCBS_check'
    output: 'results/purecn'+suffix+'/cbs_PSCBS/{sample}/{sample}.csv'
    benchmark: 'benchmarks/purecn'+suffix+'/cbs_PSCBS/{sample}.txt'
    log: 'logs/purecn'+suffix+'/cbs_PSCBS/{sample}.log',
    params:
        cnvkit_method='cbs',
        purecn_method='PSCBS',
        sampleid='{sample}',
        random_nb=randint(1,1000),

use rule purecn_hmm_PSCBS as purecn_hmm_PSCBS_w_pon with:
    input:
        vcf_filt='results/mutect2/filtered/{sample}_filtered.vcf.gz',
        copy_ratios='results/cnvkit'+suffix+'/general/{sample}.cnr',
        seg='results/cnvkit'+suffix+'/hmm/{sample}.seg',
        install='results/purecn/PSCBS_check'
    output: 'results/purecn'+suffix+'/hmm_PSCBS/{sample}/{sample}.csv'
    benchmark: 'benchmarks/purecn'+suffix+'/hmm_PSCBS/{sample}.txt'
    log: 'logs/purecn'+suffix+'/hmm_PSCBS/{sample}.log',
    params:
        cnvkit_method='hmm',
        purecn_method='PSCBS',
        sampleid='{sample}',
        random_nb=randint(1,1000),