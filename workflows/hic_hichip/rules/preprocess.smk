ruleorder: fastp_pe > prepare_hicpro_input

rule fastp_pe:
    input:
        r1=raw_read1,
        r2=raw_read2
    output:
        r1=temp(PATHS.clean_r1("{sample}")),
        r2=temp(PATHS.clean_r2("{sample}"))
    wildcard_constraints:
        sample=CTX.sample_pattern
    params:
        html=PATHS.fastp_html("{sample}"),
        json=PATHS.fastp_json("{sample}"),
        extra=lambda wildcards: (config.get("preprocessing", {}) or {}).get("fastp_extra", "")
    log:
        PATHS.log("fastp/{sample}")
    benchmark:
        PATHS.benchmark("fastp/{sample}")
    threads:
        workflow_threads("fastp", 8)
    conda:
        ENV
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.r1:q}) $(dirname {params.html:q}) $(dirname {log:q}) $(dirname {benchmark:q})
fastp \
    -i {input.r1:q} \
    -I {input.r2:q} \
    -o {output.r1:q} \
    -O {output.r2:q} \
    -w {threads} \
    -h {params.html:q} \
    -j {params.json:q} \
    {params.extra} \
    &> {log:q}
test -s {output.r1:q}
test -s {output.r2:q}
        """

rule prepare_hicpro_input:
    input:
        reads=hicpro_reads
    output:
        r1=PATHS.hicpro_input_r1("{sample}"),
        r2=PATHS.hicpro_input_r2("{sample}"),
        done=PATHS.hicpro_input_done("{sample}")
    wildcard_constraints:
        sample=CTX.sample_pattern
    log:
        PATHS.log("prepare_hicpro_input/{sample}")
    shell:
        r"""
set -euo pipefail
mkdir -p $(dirname {output.r1:q}) $(dirname {log:q})
ln -sf $(realpath {input.reads[0]:q}) {output.r1:q}
ln -sf $(realpath {input.reads[1]:q}) {output.r2:q}
test -e {output.r1:q}
test -e {output.r2:q}
date > {output.done:q}
        """
