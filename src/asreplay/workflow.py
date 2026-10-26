"""Extract SCARLET's methylation branches without unrelated analyses."""
import re
from pathlib import Path
import shutil
import subprocess
import sys
from .rapidcns2 import repair



def create(scarlet, human_variation, out):
    out=Path(out)
    source=Path(scarlet).resolve();text=(source/'main.nf').read_text()
    names=['modkit_adjust_mods','index_input_bam','index_merged_bam','human_variation_mods',
           'STURGEON_modkit_extract','STURGEON_inputtobed','STURGEON_predict',
           'nanoDx_modkit','gzip_nanoDx_modkit','nanoDx_annotate_bedMethyl','nanoDx_NN_classify',
           'rapid_cns2_meth_classification']
    blocks=[]
    asset=Path(human_variation).resolve()
    revision=subprocess.check_output(['git','-C',str(asset),'rev-parse','HEAD'],text=True).strip()
    out.mkdir(parents=True,exist_ok=False)
    # Snapshot the upstream source without its Git database, preserving the code used.
    if not (out/'wf-human-variation').exists():
        shutil.copytree(asset,out/'wf-human-variation',ignore=shutil.ignore_patterns('.git','work','.nextflow*'))
    (out/'upstream_revision.txt').write_text(revision+'\n')
    for name in names:
        m=re.search(r'^process '+name+r' \{\n',text,re.M)
        if not m: raise ValueError(name)
        nextproc=re.search(r'^process \w+ \{',text[m.end():],re.M)
        end=m.end()+nextproc.start() if nextproc else text.index('///////////////////////////\n// MAIN WORKFLOW',m.end())
        block=text[m.start():end].rstrip()
        # Trim comments between processes after the final closing brace.
        block=block[:block.rfind('}')+1]
        if name=='human_variation_mods':
            block=block.replace('nextflow run epi2me-labs/wf-human-variation -r master',
                                'nextflow -c ${projectDir}/nested.config run ${projectDir}/wf-human-variation')
        if name in ['STURGEON_predict','nanoDx_NN_classify','rapid_cns2_meth_classification']:
            block=block.replace('    input:', '    publishDir(params.outdir, mode: "copy")\n    input:',1)
        if name=='nanoDx_NN_classify':
            block=block.replace('    output:', '    output:\n        path "*nanoDx_output.txt", emit: summary',1)
        if name=='rapid_cns2_meth_classification':
            block=block.replace('    output:', '    output:\n        path "*calibrated_classification.tsv", emit: calibrated',1)
        blocks.append(block)
    # Copy the R classifier; ranger prediction otherwise defaults to all host cores.
    r=(source/'bin/methylation_classification_nanodx_v0.3.R').read_text()
    r=r.replace('predict(rf, ts, predict.all = F)', 'predict(rf, ts, predict.all = F, num.threads=opt$threads)')
    r=r.replace('predict(rf, rbind(df[,cols],df[,cols]), predict.all=F)',
                'predict(rf, rbind(df[,cols],df[,cols]), predict.all=F, num.threads=opt$threads)')
    (out/'rapidcns2.R').write_text(repair(r))
    if not (out/'bin').exists(): (out/'bin').symlink_to(source/'bin',target_is_directory=True)
    workflow=r'''
    workflow {
        input_bam = Channel.fromPath(params.bam, checkIfExists: true)
        reference = Channel.fromPath(params.reference, checkIfExists: true)
        targets = Channel.fromPath("${projectDir}/bin/NPHD_panel_hg38_clean.bed")
        mapping = Channel.fromPath("${projectDir}/bin/hglft_genome_260e9_91a970_clean.bed")
        model = Channel.fromPath("${projectDir}/bin/Capper_et_al_NN.pkl")
        classifier = Channel.fromPath("${projectDir}/bin/classify_NN_bedMethyl.py")
        rapid_script = Channel.fromPath("${projectDir}/rapidcns2.R")
        probes = Channel.fromPath("${projectDir}/bin/top_probes_hm450.Rdata")
        training = Channel.fromPath("${projectDir}/bin/capper_top_100k_betas_binarised.Rdata")
        arrayfile = Channel.fromPath("${projectDir}/bin/HM450.hg38.manifest.gencode.v22.Rdata")
        merged = modkit_adjust_mods(input_bam, params.sample, params.threads)
        indexed = index_input_bam(input_bam, params.threads)
        merged_indexed = index_merged_bam(merged.modkit_merged_bam, params.threads)
        mods = human_variation_mods(merged_indexed.indexed_bam, targets, reference,
                                   params.sample, params.outdir, 0, params.threads)
        rapid_cns2_meth_classification(rapid_script, params.sample, probes, training,
                                      arrayfile, params.threads, mods.bedmethyl_gz)
        extracted = STURGEON_modkit_extract(merged.modkit_merged_bam, params.sample, params.threads)
        converted = STURGEON_inputtobed(extracted.modkit_extract_output)
        STURGEON_predict(converted.sturgeon_bed_convert)
        pileup = nanoDx_modkit(indexed.indexed_bam_tuple, reference, mapping, params.sample, params.threads)
        zipped = gzip_nanoDx_modkit(pileup.nanoDx_modkit_out, params.sample)
        annotated = nanoDx_annotate_bedMethyl(pileup.nanoDx_modkit_out, mapping,
                                            params.sample, zipped.nanoDx_modkit_out_gz)
        nanoDx_NN_classify(classifier, model, annotated.nanoDx_annotate_out, params.sample)
    }
    '''
    (out/'main.nf').write_text('nextflow.enable.dsl=2\n'+ '\n\n'.join(blocks)+workflow)
    config='''params.threads = 16
    process {
        executor = 'local'
        cpus = params.threads
        memory = '128 GB'
        container = 'graefox/scarlet:latest'
        withName: human_variation_mods { container = null }
        env.OMP_NUM_THREADS = '1'
        env.OPENBLAS_NUM_THREADS = '1'
        env.MKL_NUM_THREADS = '1'
        env.NUMEXPR_NUM_THREADS = '1'
    }
    executor.queueSize = 1
    executor.cpus = 16
    executor.memory = '128 GB'
    docker.enabled = true
    docker.runOptions = '-u $(id -u):$(id -g) --cpus=16 --memory=120g'
    report.enabled = true
    report.overwrite = true
    '''
    # Environment limits belong at config scope, not process directives.
    config=config.replace("    env.OMP_NUM_THREADS = '1'\n    env.OPENBLAS_NUM_THREADS = '1'\n    env.MKL_NUM_THREADS = '1'\n    env.NUMEXPR_NUM_THREADS = '1'\n",'')
    config+="\nenv { OMP_NUM_THREADS='1'; OPENBLAS_NUM_THREADS='1'; MKL_NUM_THREADS='1'; NUMEXPR_NUM_THREADS='1' }\n"
    (out/'nextflow.config').write_text(config)
    (out/'nested.config').write_text("process.cpus=16\nprocess.memory='32 GB'\nexecutor.queueSize=1\nexecutor.cpus=16\ndocker.runOptions='-u $(id -u):$(id -g) --cpus=16 --memory=32g'\nenv { OMP_NUM_THREADS='1'; OPENBLAS_NUM_THREADS='1'; MKL_NUM_THREADS='1' }\n")
    return out/'main.nf'
