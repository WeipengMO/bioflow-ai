from pathlib import Path
from types import SimpleNamespace
import csv
import os
import re
import shlex
import yaml

FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
REQUIRED_SAMPLE_COLUMNS = ["sample", "fq1", "fq2", "assay", "condition", "replicate", "group"]
REQUIRED_COMPARISON_COLUMNS = ["comparison", "group1", "group2"]


def as_list(value):
    if value is None:
        return []
    if isinstance(value, str):
        return [item.strip() for item in re.split(r"[, ]+", value) if item.strip()]
    return [str(item) for item in value]


def as_bool(value, default=False):
    if value is None:
        return bool(default)
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        if value in {0, 1}:
            return bool(value)
        raise ValueError(f"Expected boolean-like value, got integer: {value}")
    normalized = str(value).strip().lower()
    if normalized in {"true", "t", "yes", "y", "1", "on"}:
        return True
    if normalized in {"false", "f", "no", "n", "0", "off", ""}:
        return False
    raise ValueError(f"Expected boolean-like value, got: {value}")


def q(path):
    return shlex.quote(str(path))


def unique_list(values):
    return list(dict.fromkeys([str(v) for v in values]))


def sample_pattern(samples):
    if not samples:
        return r"__NO_SAMPLE__"
    return "|".join(re.escape(sample) for sample in samples)


def load_tsv(path, required_columns, table_name):
    path = Path(path)
    if not path.exists():
        raise ValueError(f"{table_name} file does not exist: {path}")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{table_name} is empty: {path}")
        missing = [col for col in required_columns if col not in reader.fieldnames]
        if missing:
            raise ValueError(f"{table_name} is missing required column(s): {', '.join(missing)}")
        rows = []
        for row in reader:
            clean = {k: (v.strip() if isinstance(v, str) else v) for k, v in row.items()}
            if any(clean.values()):
                rows.append(clean)
    if not rows:
        raise ValueError(f"{table_name} has no data rows: {path}")
    return rows


class WorkflowPaths:
    def __init__(self, outdir):
        outdir = str(outdir or "results").strip() or "results"
        self.outdir = outdir if outdir == "/" else outdir.rstrip("/")
        self.clean_data = self.result("clean_data")
        self.fastp_qc = self.result("qc", "fastp")
        self.hicpro_input = self.result("hicpro_input")
        self.hicpro = self.result("hicpro")
        self.hicpro_generated_config = self.result("config", "hicpro.generated.config")
        self.valid_pairs = self.result("valid_pairs")
        self.matrix = self.result("matrix")
        self.cool = self.result("cool")
        self.bam = self.result("bam")
        self.peaks = self.result("peaks")
        self.loops = self.result("loops")
        self.consensus_loops = self.result("loops", "consensus")
        self.diffloop = self.result("diffloop")
        self.annotation = self.result("annotation")
        self.qc = self.result("qc")
        self.reports = self.result("reports")
        self.logs = self.result("logs")
        self.benchmarks = self.result("benchmarks")

    def result(self, *parts):
        return str(Path(self.outdir, *[str(part) for part in parts if str(part)]))

    def log(self, name):
        return f"{self.logs}/{name}.log"

    def benchmark(self, name):
        return f"{self.benchmarks}/{name}.txt"

    def clean_r1(self, sample):
        return f"{self.clean_data}/{sample}.R1.clean.fq.gz"

    def clean_r2(self, sample):
        return f"{self.clean_data}/{sample}.R2.clean.fq.gz"

    def fastp_html(self, sample):
        return f"{self.fastp_qc}/{sample}.html"

    def fastp_json(self, sample):
        return f"{self.fastp_qc}/{sample}.json"

    def hicpro_input_r1(self, sample):
        return f"{self.hicpro_input}/{sample}/{sample}_R1.fastq.gz"

    def hicpro_input_r2(self, sample):
        return f"{self.hicpro_input}/{sample}/{sample}_R2.fastq.gz"

    def hicpro_input_done(self, sample):
        return f"{self.hicpro_input}/{sample}/.done"

    def hicpro_done(self):
        return f"{self.hicpro}/.hicpro.done"

    def sample_valid_pairs(self, sample):
        return f"{self.valid_pairs}/{sample}.allValidPairs"

    def raw_matrix(self, sample, resolution):
        return f"{self.matrix}/{sample}/{resolution}/raw/{sample}_{resolution}.matrix"

    def raw_bins(self, sample, resolution):
        return f"{self.matrix}/{sample}/{resolution}/raw/{sample}_{resolution}_abs.bed"

    def iced_matrix(self, sample, resolution):
        return f"{self.matrix}/{sample}/{resolution}/iced/{sample}_{resolution}_iced.matrix"

    def iced_bins(self, sample, resolution):
        return f"{self.matrix}/{sample}/{resolution}/iced/{sample}_{resolution}_abs.bed"

    def cool_file(self, sample, resolution):
        return f"{self.cool}/{sample}.{resolution}.cool"

    def mcool_file(self, sample):
        return f"{self.cool}/{sample}.mcool"

    def hicpro_mapped_bam(self, sample):
        return f"{self.bam}/{sample}.hicpro.mapped.bam"

    def auto_peak_narrowpeak(self, sample):
        return f"{self.peaks}/{sample}/{sample}.auto_peaks.narrowPeak"

    def auto_peak_bed(self, sample):
        return f"{self.peaks}/{sample}/{sample}.auto_peaks.bed"

    def auto_peak_log(self, sample, caller):
        return f"{self.peaks}/{sample}/{sample}.{caller}.log"

    def peak_source_table(self):
        return f"{self.peaks}/peak_sources.tsv"

    def hicpro_qc_summary(self):
        return f"{self.qc}/hicpro_qc_summary.tsv"

    def distance_decay_pdf(self):
        return f"{self.qc}/contact_distance_decay.pdf"

    def distance_decay_tsv(self):
        return f"{self.qc}/contact_distance_decay.tsv"

    def sample_correlation_pdf(self):
        return f"{self.qc}/sample_correlation.pdf"

    def multiqc_report(self):
        return f"{self.reports}/multiqc_report.html"

    def sample_loops_bedpe(self, sample):
        return f"{self.loops}/{sample}/{sample}.loops.bedpe"

    def sample_loops_tsv(self, sample):
        return f"{self.loops}/{sample}/{sample}.loops.tsv"

    def sample_anchors_bed(self, sample):
        return f"{self.loops}/{sample}/{sample}.anchors.bed"

    def loop_universe(self):
        return f"{self.consensus_loops}/loop_universe.bedpe"

    def group_consensus_loops(self):
        return f"{self.consensus_loops}/group_consensus_loops.bedpe"

    def loop_annotation(self):
        return f"{self.consensus_loops}/loop_annotation.tsv"

    def loop_counts(self):
        return f"{self.diffloop}/counts/loop_counts.tsv"

    def loop_metadata(self):
        return f"{self.diffloop}/counts/loop_metadata.tsv"

    def sample_metadata(self):
        return f"{self.diffloop}/counts/sample_metadata.tsv"

    def diffloops(self, comparison):
        return f"{self.diffloop}/{comparison}/diffloops.tsv"

    def diffloops_significant(self, comparison):
        return f"{self.diffloop}/{comparison}/diffloops.significant.tsv"

    def diffloops_up(self, comparison):
        return f"{self.diffloop}/{comparison}/up_loops.bedpe"

    def diffloops_down(self, comparison):
        return f"{self.diffloop}/{comparison}/down_loops.bedpe"

    def diffloop_volcano(self, comparison):
        return f"{self.diffloop}/{comparison}/volcano.pdf"

    def diffloop_ma(self, comparison):
        return f"{self.diffloop}/{comparison}/ma_plot.pdf"

    def diffloop_heatmap(self, comparison):
        return f"{self.diffloop}/{comparison}/heatmap_top_diffloops.pdf"

    def loops_to_genes(self):
        return f"{self.annotation}/loops_to_genes.tsv"

    def promoter_enhancer_loops(self):
        return f"{self.annotation}/promoter_enhancer_loops.tsv"

    def workflow_report(self):
        return f"{self.reports}/hic_hichip_report.html"


PATHS = WorkflowPaths(config.get("outdir", "results"))


def resolve_path(path, base_dir=None):
    if path is None or str(path).strip() == "":
        return ""
    path = Path(str(path))
    if path.is_absolute():
        return str(path)
    if base_dir:
        return str(Path(base_dir) / path)
    return str(path)


def parse_samples(config):
    samples_path = config.get("samples", "config/samples.tsv")
    rows = load_tsv(samples_path, REQUIRED_SAMPLE_COLUMNS, "samples.tsv")
    raw_data_dir = config.get("raw_data_dir", "raw_data")
    validate_files = as_bool(config.get("validate_input_files", False), False)
    seen = set()
    records = {}
    for row in rows:
        sample = row["sample"]
        if not re.match(r"^[A-Za-z0-9_.-]+$", sample):
            raise ValueError(f"Invalid sample name '{sample}'. Use letters, numbers, dot, underscore, or dash only.")
        if sample in seen:
            raise ValueError(f"Duplicate sample name in samples.tsv: {sample}")
        seen.add(sample)
        assay = row["assay"].lower()
        if assay not in {"hic", "hichip"}:
            raise ValueError(f"Sample {sample}: assay must be hic or hichip, got {assay}")
        for key in ["condition", "replicate", "group"]:
            if not row.get(key):
                raise ValueError(f"Sample {sample}: required field '{key}' is empty")
        fq1 = resolve_path(row["fq1"], raw_data_dir)
        fq2 = resolve_path(row["fq2"], raw_data_dir)
        if validate_files:
            for fq in [fq1, fq2]:
                if not Path(fq).exists():
                    raise ValueError(f"FASTQ file does not exist for sample {sample}: {fq}")
        peak_bed = row.get("peak_bed", "")
        if peak_bed:
            peak_bed = resolve_path(peak_bed, raw_data_dir)
            if validate_files and not Path(peak_bed).exists():
                raise ValueError(f"Peak BED file does not exist for sample {sample}: {peak_bed}")
        records[sample] = {
            "sample": sample,
            "fq1": fq1,
            "fq2": fq2,
            "assay": assay,
            "condition": row["condition"],
            "replicate": row["replicate"],
            "group": row["group"],
            "target": row.get("target", ""),
            "peak_bed": peak_bed,
            "control_bam": row.get("control_bam", ""),
        }
    return records


def parse_comparisons(config, samples):
    comparisons_path = config.get("comparisons", "config/comparisons.tsv")
    if not comparisons_path or not Path(comparisons_path).exists():
        return []
    rows = load_tsv(comparisons_path, REQUIRED_COMPARISON_COLUMNS, "comparisons.tsv")
    groups = {record["group"] for record in samples.values()}
    seen = set()
    comparisons = []
    for row in rows:
        name = row["comparison"]
        if not re.match(r"^[A-Za-z0-9_.-]+$", name):
            raise ValueError(f"Invalid comparison name '{name}'. Use letters, numbers, dot, underscore, or dash only.")
        if name in seen:
            raise ValueError(f"Duplicate comparison name: {name}")
        seen.add(name)
        if row["group1"] not in groups:
            raise ValueError(f"Comparison {name}: group1 not found in samples.tsv: {row['group1']}")
        if row["group2"] not in groups:
            raise ValueError(f"Comparison {name}: group2 not found in samples.tsv: {row['group2']}")
        comparisons.append({"comparison": name, "group1": row["group1"], "group2": row["group2"]})
    return comparisons


def workflow_threads(name, default):
    return int((config.get("threads", {}) or {}).get(name, default))


def build_workflow_context(config):
    samples = parse_samples(config)
    sample_names = unique_list(samples.keys())
    comparisons = parse_comparisons(config, samples)
    hicpro_cfg = config.get("hicpro", {}) or {}
    loop_cfg = config.get("loop_calling", {}) or {}
    diff_cfg = config.get("diffloop", {}) or {}
    qc_cfg = config.get("qc", {}) or {}
    preprocessing_cfg = config.get("preprocessing", {}) or {}
    annotation_cfg = config.get("annotation", {}) or {}
    peak_cfg = config.get("peak_calling", {}) or {}
    resolutions = [str(int(x)) for x in as_list(hicpro_cfg.get("resolutions", [10000, 25000, 50000]))]
    if not resolutions:
        raise ValueError("hicpro.resolutions must contain at least one resolution")
    primary_resolution = str(int(loop_cfg.get("resolution") or min(int(r) for r in resolutions)))
    if primary_resolution not in resolutions:
        raise ValueError("loop_calling.resolution must be one of hicpro.resolutions")
    enable_loops = as_bool(loop_cfg.get("enable", True), True)
    enable_diffloop = as_bool(diff_cfg.get("enable", True), True) and len(comparisons) > 0
    enable_fastp = as_bool(preprocessing_cfg.get("enable_fastp", False), False)
    enable_multiqc = as_bool(qc_cfg.get("enable_multiqc", True), True)
    enable_distance_decay = as_bool(qc_cfg.get("enable_distance_decay", True), True)
    enable_correlation = as_bool(qc_cfg.get("enable_correlation", True), True)
    enable_annotation = as_bool(annotation_cfg.get("enable", False), False)
    enable_auto_peak = as_bool(peak_cfg.get("enable_auto_peak", True), True)
    peak_caller = str(peak_cfg.get("caller", "macs2")).lower()
    if peak_caller not in {"macs2", "macs3"}:
        raise ValueError("peak_calling.caller must be macs2 or macs3")
    loop_caller = str(loop_cfg.get("caller", "fithichip")).lower()
    if loop_caller not in {"fithichip", "mustache", "cooltools", "precomputed"}:
        raise ValueError("loop_calling.caller must be one of fithichip, mustache, cooltools, precomputed")
    if loop_caller == "fithichip":
        hic_samples = [s for s, rec in samples.items() if rec["assay"] == "hic"]
        if hic_samples:
            raise ValueError("loop_calling.caller=fithichip is for HiChIP. Hi-C sample(s) cannot use FitHiChIP: " + ", ".join(hic_samples))
        missing_peak = []
        for sample, rec in samples.items():
            has_sample_peak = bool(rec.get("peak_bed"))
            has_global_peak = bool(loop_cfg.get("external_peaks", ""))
            can_auto = rec["assay"] == "hichip" and enable_auto_peak
            if not (has_sample_peak or has_global_peak or can_auto):
                missing_peak.append(sample)
        if missing_peak:
            raise ValueError(
                "FitHiChIP requires peaks. Provide samples.tsv peak_bed, set loop_calling.external_peaks, "
                "or enable peak_calling.enable_auto_peak for sample(s): " + ", ".join(missing_peak)
            )
    for sample, rec in samples.items():
        if rec["assay"] == "hichip" and not rec.get("target"):
            print(f"[BioFlowAI hic_hichip] WARNING: HiChIP sample {sample} has empty target.")
    bind_paths = unique_list(as_list(hicpro_cfg.get("extra_bind_paths", [])))
    for key in ["fasta", "chrom_sizes", "restriction_fragments", "bowtie2_index"]:
        value = (config.get("genome", {}) or {}).get(key)
        if value:
            bind_paths.append(str(Path(value).parent if Path(str(value)).suffix else Path(value).parent))
    bind_paths.append(str(Path.cwd()))
    bind_paths = unique_list([p for p in bind_paths if p and p != "."])
    return SimpleNamespace(
        samples=samples,
        sample_names=sample_names,
        comparisons=comparisons,
        comparison_names=unique_list([c["comparison"] for c in comparisons]),
        sample_pattern=sample_pattern(sample_names),
        comparison_pattern=sample_pattern([c["comparison"] for c in comparisons]),
        resolutions=resolutions,
        primary_resolution=primary_resolution,
        enable_fastp=enable_fastp,
        enable_loops=enable_loops,
        enable_diffloop=enable_diffloop,
        enable_multiqc=enable_multiqc,
        enable_distance_decay=enable_distance_decay,
        enable_correlation=enable_correlation,
        enable_annotation=enable_annotation,
        hicpro_container=str(hicpro_cfg.get("container", "/home/mowp/software/HiC-Pro/hicpro3.sif")),
        hicpro_template=str(hicpro_cfg.get("config_template", "config/hicpro.config.template")),
        hicpro_bind_paths=bind_paths,
        hicpro_singularity_args=str(hicpro_cfg.get("singularity_args", "")),
        base_resolution=str(min(int(r) for r in resolutions)),
        loop_caller=loop_caller,
        loop_fdr=str(loop_cfg.get("fdr", 0.01)),
        peak_mode=str(loop_cfg.get("peak_mode", "external_peak")),
        enable_auto_peak=enable_auto_peak,
        peak_caller=peak_caller,
        peak_genome_size=str(peak_cfg.get("genome_size", "hs")),
        peak_qvalue=str(peak_cfg.get("qvalue", 0.01)),
        peak_mode_macs=str(peak_cfg.get("mode", "BAMPE")),
        peak_extra=str(peak_cfg.get("extra", "")),
        peak_allow_empty=as_bool(peak_cfg.get("allow_empty", False), False),
        quantification_source=str(diff_cfg.get("quantification_source", "cool")).lower(),
        use_balanced_counts=as_bool(diff_cfg.get("use_balanced_counts", False), False),
        loop_universe_mode=str(diff_cfg.get("loop_universe", "union")),
        anchor_slop=int(diff_cfg.get("anchor_slop", 0)),
        group_consensus_min_replicates=int(diff_cfg.get("group_consensus_min_replicates", 1)),
        group_consensus_min_fraction=float(diff_cfg.get("group_consensus_min_fraction", 0.5)),
        diff_fdr=float(diff_cfg.get("fdr", 0.05)),
        diff_lfc_cutoff=float(diff_cfg.get("lfc_cutoff", 0.0)),
        diff_design_formula=str(diff_cfg.get("design_formula", diff_cfg.get("design", "~ group"))),
    )


def print_workflow_summary(ctx):
    print("[BioFlowAI hic_hichip] samples:", ", ".join(ctx.sample_names))
    print("[BioFlowAI hic_hichip] resolutions:", ", ".join(ctx.resolutions))
    print("[BioFlowAI hic_hichip] loop calling:", ctx.enable_loops, ctx.loop_caller)
    print("[BioFlowAI hic_hichip] diffloop:", ctx.enable_diffloop)
    print("[BioFlowAI hic_hichip] HiC-Pro container:", ctx.hicpro_container)


def raw_read1(wildcards):
    return CTX.samples[wildcards.sample]["fq1"]


def raw_read2(wildcards):
    return CTX.samples[wildcards.sample]["fq2"]


def hicpro_reads(wildcards):
    if CTX.enable_fastp:
        return [PATHS.clean_r1(wildcards.sample), PATHS.clean_r2(wildcards.sample)]
    return [CTX.samples[wildcards.sample]["fq1"], CTX.samples[wildcards.sample]["fq2"]]


def sample_peak_bed(wildcards):
    return sample_peak_path(wildcards.sample)


def sample_external_peak(sample):
    rec = CTX.samples[sample]
    if rec.get("peak_bed"):
        return rec["peak_bed"]
    global_peak = (config.get("loop_calling", {}) or {}).get("external_peaks", "")
    return str(global_peak or "")


def needs_auto_peak(sample):
    rec = CTX.samples[sample]
    return bool(
        CTX.enable_auto_peak
        and CTX.loop_caller == "fithichip"
        and rec["assay"] == "hichip"
        and not sample_external_peak(sample)
    )


def sample_peak_path(sample):
    external = sample_external_peak(sample)
    if external:
        return external
    if needs_auto_peak(sample):
        return PATHS.auto_peak_bed(sample)
    return ""


def sample_peak_source(sample):
    rec = CTX.samples[sample]
    if rec.get("peak_bed"):
        return "sample_peak_bed"
    if (config.get("loop_calling", {}) or {}).get("external_peaks", ""):
        return "global_external_peak"
    if needs_auto_peak(sample):
        return "hichip_auto_peak"
    return "none"


def sample_peak_input(wildcards):
    if needs_auto_peak(wildcards.sample):
        return PATHS.auto_peak_bed(wildcards.sample)
    peak = sample_external_peak(wildcards.sample)
    return peak if peak else []


def comparison_group1(wildcards):
    return next(c["group1"] for c in CTX.comparisons if c["comparison"] == wildcards.comparison)


def comparison_group2(wildcards):
    return next(c["group2"] for c in CTX.comparisons if c["comparison"] == wildcards.comparison)


def all_outputs(ctx):
    outputs = []
    outputs.append(PATHS.hicpro_done())
    outputs.extend([PATHS.sample_valid_pairs(sample) for sample in ctx.sample_names])
    outputs.extend([PATHS.mcool_file(sample) for sample in ctx.sample_names])
    outputs.append(PATHS.hicpro_qc_summary())
    if ctx.enable_distance_decay:
        outputs.extend([PATHS.distance_decay_pdf(), PATHS.distance_decay_tsv()])
    if ctx.enable_correlation:
        outputs.append(PATHS.sample_correlation_pdf())
    if ctx.enable_multiqc:
        outputs.append(PATHS.multiqc_report())
    if ctx.enable_loops:
        auto_peak_samples = [sample for sample in ctx.sample_names if needs_auto_peak(sample)]
        outputs.extend([PATHS.hicpro_mapped_bam(sample) for sample in auto_peak_samples])
        outputs.extend([PATHS.auto_peak_bed(sample) for sample in auto_peak_samples])
        outputs.extend([PATHS.auto_peak_narrowpeak(sample) for sample in auto_peak_samples])
        outputs.append(PATHS.peak_source_table())
        outputs.extend([PATHS.sample_loops_bedpe(sample) for sample in ctx.sample_names])
        outputs.extend([PATHS.sample_loops_tsv(sample) for sample in ctx.sample_names])
        outputs.extend([PATHS.sample_anchors_bed(sample) for sample in ctx.sample_names])
        outputs.extend([PATHS.loop_universe(), PATHS.group_consensus_loops(), PATHS.loop_annotation()])
    if ctx.enable_diffloop:
        outputs.extend([PATHS.loop_counts(), PATHS.loop_metadata(), PATHS.sample_metadata()])
        for comp in ctx.comparison_names:
            outputs.extend([
                PATHS.diffloops(comp),
                PATHS.diffloops_significant(comp),
                PATHS.diffloops_up(comp),
                PATHS.diffloops_down(comp),
                PATHS.diffloop_volcano(comp),
                PATHS.diffloop_ma(comp),
                PATHS.diffloop_heatmap(comp),
            ])
    if ctx.enable_annotation:
        outputs.extend([PATHS.loops_to_genes(), PATHS.promoter_enhancer_loops()])
    outputs.append(PATHS.workflow_report())
    return outputs
