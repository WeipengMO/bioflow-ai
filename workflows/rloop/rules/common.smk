from pathlib import Path
from types import SimpleNamespace
import csv
import re
import shlex
import yaml


FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
PAIRED_FASTQ_RE = re.compile(
    r"^(?P<sample>.+?)(?:[._-]R?|[._-])(?P<read>[12])(?:_001)?(?P<suffix>\.f(?:ast)?q(?:\.gz)?)$",
    re.IGNORECASE,
)
ALLOWED_ROLES = {"treatment", "control", "rnaseh_control"}


class WorkflowPaths:
    def __init__(self, outdir):
        outdir = str(outdir or "results").strip() or "results"
        self.outdir = outdir if outdir == "/" else outdir.rstrip("/")

        self.clean_data = self.result("clean_data")
        self.aligned_data = self.result("aligned_data")
        self.fastp_qc = self.result("qc", "fastp")
        self.mark_duplicates_qc = self.result("qc", "mark_duplicates")
        self.insert_size_qc = self.result("qc", "insert_size")
        self.frip_qc = self.result("qc", "frip")
        self.bigwig_scale_qc = self.result("qc", "bigwig_scale")
        self.rnaseh_sensitivity_qc = self.result("qc", "rnaseh_sensitivity")
        self.bigwig = self.result("bigwig")
        self.bigwig_common_scale = self.result("bigwig_common_scale")
        self.bigwig_debug_raw = self.result("bigwig_debug", "raw_scale")
        self.bigwig_debug_dedup = self.result("bigwig_debug", "dedup_effective_fragments")
        self.macs3_results = self.result("macs3_results")
        self.rnaseh_no_overlap = self.result("rnaseh_no_overlap")
        self.rnaseh_filtered = self.result("rnaseh_filtered")
        self.rnaseh_signal = self.result("rnaseh_signal")
        self.replicate_consensus = self.result("replicate_consensus")
        self.reports = self.result("reports")
        self.logs = self.result("logs")

    def result(self, *parts):
        return str(Path(self.outdir, *[str(part) for part in parts if str(part)]))

    def clean_r1(self, sample, mode):
        suffix = ".R1.clean.fq.gz" if mode == "pe" else ".clean.fq.gz"
        return f"{self.clean_data}/{sample}{suffix}"

    def clean_r2(self, sample):
        return f"{self.clean_data}/{sample}.R2.clean.fq.gz"

    def sorted_bam(self, sample):
        return f"{self.aligned_data}/{sample}.sorted.bam"

    def sorted_bai(self, sample):
        return f"{self.aligned_data}/{sample}.sorted.bam.bai"

    def alignment_filtered_bam(self, sample):
        return f"{self.aligned_data}/{sample}.sorted.filtered.bam"

    def alignment_filtered_bai(self, sample):
        return f"{self.aligned_data}/{sample}.sorted.filtered.bam.bai"

    def markdup_bam(self, sample):
        return f"{self.aligned_data}/{sample}.sorted.markdup.bam"

    def markdup_bai(self, sample):
        return f"{self.aligned_data}/{sample}.sorted.markdup.bam.bai"

    def dedup_bam(self, sample):
        return f"{self.aligned_data}/{sample}.sorted.rmdup.bam"

    def dedup_bai(self, sample):
        return f"{self.aligned_data}/{sample}.sorted.rmdup.bam.bai"

    def signal_bam(self, sample):
        return f"{self.aligned_data}/{sample}.sorted.markdup.filtered.bam"

    def signal_bai(self, sample):
        return f"{self.aligned_data}/{sample}.sorted.markdup.filtered.bam.bai"

    def filtered_bam(self, sample):
        return f"{self.aligned_data}/{sample}.sorted.rmdup.filtered.bam"

    def filtered_bai(self, sample):
        return f"{self.aligned_data}/{sample}.sorted.rmdup.filtered.bam.bai"

    def peak_bam(self, sample, peak_duplicate_mode):
        """Return the BAM for peak calling: markdup when keep_marked/auto, rmdup when remove."""
        if peak_duplicate_mode in ("keep_marked", "auto"):
            return self.signal_bam(sample)
        return self.filtered_bam(sample)

    def peak_bai(self, sample, peak_duplicate_mode):
        if peak_duplicate_mode in ("keep_marked", "auto"):
            return self.signal_bai(sample)
        return self.filtered_bai(sample)

    def spikein_bam(self, sample):
        return f"{self.aligned_data}/{sample}.spikein.sorted.bam"

    def spikein_bai(self, sample):
        return f"{self.aligned_data}/{sample}.spikein.sorted.bam.bai"

    def spikein_counts(self, sample):
        return f"{self.bigwig_scale_qc}/{sample}.spikein_counts.tsv"

    def fastp_html(self, sample):
        return f"{self.fastp_qc}/{sample}.html"

    def fastp_json(self, sample):
        return f"{self.fastp_qc}/{sample}.json"

    def mark_duplicates_metrics(self, sample):
        return f"{self.mark_duplicates_qc}/{sample}.metrics.txt"

    def insert_size_metrics(self, sample):
        return f"{self.insert_size_qc}/{sample}.insert_size_metrics.txt"

    def insert_size_histogram(self, sample):
        return f"{self.insert_size_qc}/{sample}.insert_size_histogram.pdf"

    def frip(self, sample):
        return f"{self.frip_qc}/{sample}.frip.txt"

    def bigwig_track(self, sample):
        return f"{self.bigwig}/{sample}.sorted.markdup.filtered.CPM.bw"

    def common_scale_bigwig_track(self, sample):
        return f"{self.bigwig_common_scale}/{sample}.sorted.markdup.filtered.scaled.bw"

    def debug_raw_bigwig_track(self, sample):
        return f"{self.bigwig_debug_raw}/{sample}.sorted.markdup.filtered.raw.bw"

    def debug_dedup_bigwig_track(self, sample):
        return f"{self.bigwig_debug_dedup}/{sample}.sorted.rmdup.filtered.scaled.bw"

    def signal_scale_factors(self):
        return f"{self.bigwig_scale_qc}/signal_scale_factors.{CTX.signal_scale_factor_method}.tsv"

    def dedup_scale_factors(self):
        return f"{self.bigwig_scale_qc}/dedup_scale_factors.effective_fragments.tsv"

    def bigwig_header_qc(self, sample):
        return f"{self.bigwig_scale_qc}/{sample}.bigwig_header.tsv"

    def peak(self, sample):
        peak_type = str(config.get("peak_type", "broad")).lower()
        if peak_type == "narrow":
            return f"{self.macs3_results}/narrow/{sample}_peaks.narrowPeak"
        return f"{self.macs3_results}/broad/{sample}_peaks.broadPeak"

    def rnaseh_no_overlap_peak(self, sample):
        peak_type = str(config.get("peak_type", "broad")).lower()
        suffix = "narrowPeak" if peak_type == "narrow" else "broadPeak"
        return f"{self.rnaseh_no_overlap}/{sample}_rnaseh_no_overlap.{suffix}"

    def rnaseh_sensitive_peak(self, sample):
        peak_type = str(config.get("peak_type", "broad")).lower()
        suffix = "narrowPeak" if peak_type == "narrow" else "broadPeak"
        return f"{self.rnaseh_filtered}/{sample}_rnaseh_sensitive.{suffix}"

    def rnaseh_signal_table(self, sample):
        return f"{self.rnaseh_signal}/{sample}_rnaseh_signal.tsv"

    def rnaseh_depleted_peak(self, sample):
        return f"{self.rnaseh_signal}/{sample}_rnaseh_depleted.bed"

    def rnaseh_sensitivity_summary(self, sample):
        return f"{self.rnaseh_sensitivity_qc}/{sample}.summary.tsv"

    def final_peak(self, sample):
        if CTX.enable_rnaseh_subtraction and CTX.rnaseh_controls.get(sample):
            return self.rnaseh_no_overlap_peak(sample)
        return self.peak(sample)

    def replicate_consensus_peak(self, group):
        return f"{self.replicate_consensus}/{group}_consensus.bed"

    def multiqc_report(self):
        return f"{self.reports}/multiqc_report.html"

    def log(self, name):
        return f"{self.logs}/{name}.log"


PATHS = WorkflowPaths(config.get("outdir", "results"))


def as_list(value):
    if value is None:
        return []
    if isinstance(value, str):
        return [item.strip() for item in value.split(",") if item.strip()]
    return [str(item) for item in value]


def unique_list(values):
    return list(dict.fromkeys(str(value) for value in values))


def validate_choices(values, allowed, key):
    invalid = sorted(set(values) - set(allowed))
    if invalid:
        raise ValueError(f"config['{key}'] has invalid value(s): {', '.join(invalid)}")


def config_bool(key, default=False):
    value = config.get(key, default)
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    if isinstance(value, int):
        if value in {0, 1}:
            return bool(value)
        raise ValueError(f"config['{key}'] must be a boolean value.")
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "t", "yes", "y", "1", "on"}:
            return True
        if normalized in {"false", "f", "no", "n", "0", "off", ""}:
            return False
    raise ValueError(f"config['{key}'] must be a boolean value.")


def workflow_threads(name, default):
    threads = config.get("threads", {}) or {}
    return int(threads.get(name, default))


def sample_pattern(samples):
    return "|".join(re.escape(sample) for sample in samples) if samples else r"__no_samples__"


def load_yaml(path):
    path = Path(path)
    if not path.exists():
        raise ValueError(f"Required YAML file does not exist: {path}")
    with path.open() as handle:
        return yaml.safe_load(handle) or {}


def load_sample_data(samples):
    if not isinstance(samples, dict):
        raise ValueError("sample_config must contain a mapping of sample name to sample metadata.")

    records = {}
    groups = {}
    controls = {}
    rnaseh_controls = {}
    treatments = []

    for sample, value in samples.items():
        if value is None:
            value = {}
        if not isinstance(value, dict):
            raise ValueError(f"sample_config entry '{sample}' must be a mapping.")
        sample = str(sample)
        role = str(value.get("role", "treatment")).strip().lower()
        if role not in ALLOWED_ROLES:
            raise ValueError(f"sample_config entry '{sample}' has invalid role '{role}'.")
        group = value.get("group")
        records[sample] = {"role": role, "group": str(group) if group else None}
        if role == "treatment":
            treatments.append(sample)
            if group:
                groups.setdefault(str(group), []).append(sample)
            if value.get("control"):
                controls[sample] = str(value["control"])
            if value.get("rnaseh_control"):
                rnaseh_controls[sample] = str(value["rnaseh_control"])

    missing_controls = sorted(set(controls.values()) - set(records))
    if missing_controls:
        raise ValueError("Treatment control sample(s) are not declared: " + ", ".join(missing_controls))
    missing_rnaseh = sorted(set(rnaseh_controls.values()) - set(records))
    if missing_rnaseh:
        raise ValueError("RNase H control sample(s) are not declared: " + ", ".join(missing_rnaseh))

    return SimpleNamespace(
        records=records,
        samples=unique_list(records.keys()),
        treatments=unique_list(treatments),
        groups={group: unique_list(members) for group, members in groups.items()},
        controls=controls,
        rnaseh_controls=rnaseh_controls,
    )


def empty_fastq_record():
    return {"r1": None, "r2": None}


def set_fastq_record(index, sample, read, path):
    record = index.setdefault(sample, empty_fastq_record())
    key = f"r{read}"
    if record.get(key) is not None:
        raise ValueError(
            f"Duplicate FASTQ for sample '{sample}' read {read}: {record[key]} and {path}. "
            "Use config['fastq_manifest'] to disambiguate irregular file names."
        )
    record[key] = str(path)


def load_fastq_manifest(path, mode):
    index = {}
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "read1"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError("FASTQ manifest is missing column(s): " + ", ".join(sorted(missing)))

        for row in reader:
            sample = row["sample"].strip()
            if not sample:
                continue
            read1 = row["read1"].strip()
            read2 = (row.get("read2") or "").strip()
            if mode == "pe" and not read2:
                raise ValueError(f"FASTQ manifest sample '{sample}' is missing read2 in paired-end mode.")
            index[sample] = {"r1": read1, "r2": read2 or None}
    return index


def strip_fastq_suffix(filename):
    for suffix in sorted(FASTQ_SUFFIXES, key=len, reverse=True):
        if filename.lower().endswith(suffix):
            return filename[: -len(suffix)]
    return filename


def discover_fastqs(raw_dir, mode):
    raw_dir = Path(raw_dir)
    index = {}
    if not raw_dir.exists():
        return index

    for path in sorted(p for p in raw_dir.iterdir() if p.is_file()):
        if not path.name.lower().endswith(FASTQ_SUFFIXES):
            continue
        if mode == "pe":
            match = PAIRED_FASTQ_RE.match(path.name)
            if not match:
                continue
            set_fastq_record(index, match.group("sample"), match.group("read"), path)
        else:
            index[strip_fastq_suffix(path.name)] = {"r1": str(path), "r2": None}
    return index


def load_fastq_index(config, mode):
    manifest = config.get("fastq_manifest", "")
    index = load_fastq_manifest(manifest, mode) if manifest else discover_fastqs(config.get("raw_data_dir", "raw_data"), mode)

    for sample, reads in sorted(index.items()):
        if not reads.get("r1"):
            raise ValueError(f"Sample '{sample}' is missing read1 FASTQ.")
        if mode == "pe" and not reads.get("r2"):
            raise ValueError(f"Sample '{sample}' is missing read2 FASTQ.")
    return index


def validate_sample_fastqs(samples, fastqs, mode):
    missing = sorted(sample for sample in samples if sample not in fastqs)
    if missing:
        raise ValueError(
            "Missing FASTQ files for sample(s): "
            + ", ".join(missing)
            + ". Check raw_data_dir, file naming, sample_config, or fastq_manifest."
        )

    if mode == "pe":
        missing_r2 = sorted(sample for sample in samples if not fastqs[sample].get("r2"))
        if missing_r2:
            raise ValueError("Missing read2 FASTQ for sample(s): " + ", ".join(missing_r2))


def valid_groups(group_data, samples):
    sample_set = set(samples)
    return sorted(
        group
        for group, members in group_data.items()
        if len(members) >= 2 and all(sample in sample_set for sample in members)
    )


def parse_positive_int_config(key, default):
    try:
        value = int(config.get(key, default))
    except (TypeError, ValueError):
        raise ValueError(f"config['{key}'] must be an integer.")
    if value < 1:
        raise ValueError(f"config['{key}'] must be >= 1.")
    return value


def parse_nonnegative_float_config(key, default):
    try:
        value = float(config.get(key, default))
    except (TypeError, ValueError):
        raise ValueError(f"config['{key}'] must be a number.")
    if value < 0:
        raise ValueError(f"config['{key}'] must be >= 0.")
    return value


def parse_positive_float_config(key, default):
    value = parse_nonnegative_float_config(key, default)
    if value <= 0:
        raise ValueError(f"config['{key}'] must be > 0.")
    return value


def normalize_common_scale_normalization(value):
    raw = str(value or "None").strip()
    if raw.lower() in {"none", "raw"}:
        return "None"
    if raw.lower() == "rpgc":
        return "RPGC"
    raise ValueError("config['common_scale_normalization'] must be one of: None, raw, RPGC.")


def parse_scale_factors(path, samples):
    path = str(path or "").strip()
    if not path:
        raise ValueError(
            "enable_common_scale_bigwig requires config['signal_scale_factors_tsv'] with columns: sample, scale_factor."
        )
    tsv = Path(path)
    if not tsv.exists():
        raise ValueError(f"signal_scale_factors_tsv does not exist: {tsv}")

    factors = {}
    with tsv.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = {"sample", "scale_factor"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "signal_scale_factors_tsv is missing column(s): " + ", ".join(sorted(missing))
            )
        for row in reader:
            sample = (row.get("sample") or "").strip()
            if not sample:
                continue
            try:
                scale_factor = float(row.get("scale_factor", ""))
            except ValueError:
                raise ValueError(f"Invalid scale_factor for sample '{sample}' in {tsv}.")
            if scale_factor <= 0:
                raise ValueError(f"scale_factor must be > 0 for sample '{sample}' in {tsv}.")
            factors[sample] = scale_factor

    missing_samples = sorted(set(samples) - set(factors))
    if missing_samples:
        raise ValueError(
            "signal_scale_factors_tsv is missing sample(s): " + ", ".join(missing_samples)
        )
    return factors




def normalize_scale_factor_method(value):
    method = str(value or "effective_fragments").strip().lower()
    validate_choices([method], {"effective_fragments", "raw", "tsv", "spikein"}, "signal_scale_factor_method")
    return method


def normalize_duplicate_mode(key, default, allowed):
    mode = str(config.get(key, default) or default).strip().lower()
    validate_choices([mode], allowed, key)
    return mode


def validate_mark_duplicates_extra():
    extra = str(config.get("mark_duplicates_extra", "") or "")
    if "REMOVE_DUPLICATES" in extra.upper():
        raise ValueError(
            "config['mark_duplicates_extra'] must not contain REMOVE_DUPLICATES after duplicate branch split. "
            "Use peak_duplicate_mode and signal_duplicate_mode instead."
        )
    return extra

def build_workflow_context(config):
    mode = str(config.get("mode", "pe")).lower()
    validate_choices([mode], {"pe", "se"}, "mode")
    peak_type = str(config.get("peak_type", "broad")).lower()
    validate_choices([peak_type], {"broad", "narrow"}, "peak_type")

    sample_metadata = load_yaml(config.get("sample_config", "config/samples.yml"))
    sample_data = load_sample_data(sample_metadata)
    if not sample_data.treatments:
        raise ValueError("sample_config must contain at least one role: treatment sample.")

    fastqs = load_fastq_index(config, mode)
    validate_sample_fastqs(sample_data.samples, fastqs, mode)

    groups = valid_groups(sample_data.groups, sample_data.treatments)
    consensus_min_support_value = parse_positive_int_config("consensus_min_support", 2)
    enable_common_scale_bigwig = config_bool("enable_common_scale_bigwig", False)
    common_scale_normalization = normalize_common_scale_normalization(config.get("common_scale_normalization", "None"))
    signal_scale_factor_method = normalize_scale_factor_method(config.get("signal_scale_factor_method", "effective_fragments"))
    if signal_scale_factor_method == "tsv":
        parse_scale_factors(config.get("signal_scale_factors_tsv", ""), sample_data.samples)
    spikein_method = ""
    spikein_reference_sample = ""
    if signal_scale_factor_method == "spikein":
        spikein_genome = str(config.get("spikein_genome", "") or "").strip()
        if not spikein_genome:
            raise ValueError("signal_scale_factor_method: spikein requires config['spikein_genome'].")
        spikein_method = str(config.get("spikein_method", "ratio") or "ratio").strip().lower()
        validate_choices([spikein_method], {"ratio", "rpm"}, "spikein_method")
        if spikein_method == "ratio":
            spikein_reference_sample = str(config.get("spikein_reference_sample", "") or "").strip()
            if not spikein_reference_sample:
                raise ValueError("spikein_method: ratio requires config['spikein_reference_sample'].")
            if spikein_reference_sample not in sample_data.samples:
                raise ValueError(f"spikein_reference_sample '{spikein_reference_sample}' is not a declared sample.")
    mark_duplicates_extra = validate_mark_duplicates_extra()
    signal_duplicate_mode = normalize_duplicate_mode("signal_duplicate_mode", "keep_marked", {"keep_marked"})
    peak_duplicate_mode = normalize_duplicate_mode("peak_duplicate_mode", "auto", {"remove", "keep_marked", "auto"})

    return SimpleNamespace(
        mode=mode,
        peak_type=peak_type,
        assay=str(config.get("assay", "R-loop profiling")),
        samples=sample_data.samples,
        sample_pattern=sample_pattern(sample_data.samples),
        treatments=sample_data.treatments,
        treatment_pattern=sample_pattern(sample_data.treatments),
        fastqs=fastqs,
        sample_records=sample_data.records,
        group_data=sample_data.groups,
        groups=groups,
        group_pattern=sample_pattern(groups) if groups else r"__no_valid_groups__",
        controls=sample_data.controls,
        rnaseh_controls=sample_data.rnaseh_controls,
        consensus_min_support=consensus_min_support_value,
        enable_rnaseh_subtraction=config_bool("enable_rnaseh_subtraction", True),
        enable_common_scale_bigwig=enable_common_scale_bigwig,
        common_scale_normalization=common_scale_normalization,
        signal_scale_factor_method=signal_scale_factor_method,
        spikein_method=spikein_method,
        spikein_reference_sample=spikein_reference_sample,
        dup_threshold_for_dedup=parse_nonnegative_float_config("dup_threshold_for_dedup", 0.90),
        mark_duplicates_extra=mark_duplicates_extra,
        signal_duplicate_mode=signal_duplicate_mode,
        peak_duplicate_mode=peak_duplicate_mode,
        write_deprecated_rnaseh_sensitive_alias=config_bool("write_deprecated_rnaseh_sensitive_alias", False),
        enable_bigwig_debug_tracks=config_bool("enable_bigwig_debug_tracks", False),
        rnaseh_signal_min_fold_change=parse_positive_float_config("rnaseh_signal_min_fold_change", 2.0),
        rnaseh_signal_min_treatment_signal=parse_nonnegative_float_config("rnaseh_signal_min_treatment_signal", 0.0),
        enable_insert_size_qc=config_bool("enable_insert_size_qc", mode == "pe"),
        enable_multiqc=config_bool("enable_multiqc", True),
    )


def print_workflow_summary(ctx):
    print("Assay: " + ctx.assay)
    print("Mode: " + ctx.mode)
    print("Peak type: " + ctx.peak_type)
    print("Outdir: " + PATHS.outdir)
    print("Samples: " + ",".join(ctx.samples))
    print("Treatments: " + ",".join(ctx.treatments))
    print("Replicate groups: " + (",".join(ctx.groups) if ctx.groups else "none"))
    print("RNase H no-overlap filtering: " + ("enabled" if ctx.enable_rnaseh_subtraction else "disabled"))
    print("Common-scale bigWigs: " + ("enabled (" + ctx.common_scale_normalization + ", scale=" + ctx.signal_scale_factor_method + ")" if ctx.enable_common_scale_bigwig else "disabled"))
    print("Signal duplicate mode: " + ctx.signal_duplicate_mode)
    print("Peak duplicate mode: " + ctx.peak_duplicate_mode + (" (use markdup BAM; warn if dup rate >= " + str(int(ctx.dup_threshold_for_dedup * 100)) + "%)" if ctx.peak_duplicate_mode == "auto" else ""))
    if ctx.signal_scale_factor_method == "spikein":
        print("Spike-in method: " + ctx.spikein_method + (", reference sample=" + ctx.spikein_reference_sample if ctx.spikein_method == "ratio" else ""))
    print("BigWig debug tracks: " + ("enabled" if ctx.enable_bigwig_debug_tracks else "disabled"))
    print("MultiQC: " + ("enabled" if ctx.enable_multiqc else "disabled"))


def all_outputs(ctx):
    outputs = [
        expand(PATHS.signal_bam("{sample}"), sample=ctx.samples),
        expand(PATHS.signal_bai("{sample}"), sample=ctx.samples),
        expand(PATHS.bigwig_track("{sample}"), sample=ctx.samples),
        expand(PATHS.bigwig_header_qc("{sample}"), sample=ctx.samples),
        expand(PATHS.peak("{sample}"), sample=ctx.treatments),
        expand(PATHS.frip("{sample}"), sample=ctx.treatments),
        expand(PATHS.replicate_consensus_peak("{group}"), group=ctx.groups),
    ]
    # Dedup BAMs are needed when peak calling removes duplicates or debug tracks are enabled.
    if ctx.peak_duplicate_mode == "remove" or ctx.enable_bigwig_debug_tracks:
        outputs.append(expand(PATHS.filtered_bam("{sample}"), sample=ctx.samples))
        outputs.append(expand(PATHS.filtered_bai("{sample}"), sample=ctx.samples))
    if ctx.signal_scale_factor_method == "spikein":
        outputs.append(expand(PATHS.spikein_bam("{sample}"), sample=ctx.samples))
        outputs.append(expand(PATHS.spikein_counts("{sample}"), sample=ctx.samples))
    if ctx.enable_common_scale_bigwig:
        outputs.append(PATHS.signal_scale_factors())
        outputs.append(expand(PATHS.common_scale_bigwig_track("{sample}"), sample=ctx.samples))
    if ctx.enable_bigwig_debug_tracks:
        outputs.append(PATHS.dedup_scale_factors())
        outputs.append(expand(PATHS.debug_raw_bigwig_track("{sample}"), sample=ctx.samples))
        outputs.append(expand(PATHS.debug_dedup_bigwig_track("{sample}"), sample=ctx.samples))
    if ctx.enable_rnaseh_subtraction:
        rnaseh_treatments = [sample for sample in ctx.treatments if ctx.rnaseh_controls.get(sample)]
        outputs.append(expand(PATHS.rnaseh_no_overlap_peak("{sample}"), sample=rnaseh_treatments))
        if ctx.write_deprecated_rnaseh_sensitive_alias:
            outputs.append(expand(PATHS.rnaseh_sensitive_peak("{sample}"), sample=rnaseh_treatments))
        outputs.append(expand(PATHS.rnaseh_signal_table("{sample}"), sample=rnaseh_treatments))
        outputs.append(expand(PATHS.rnaseh_depleted_peak("{sample}"), sample=rnaseh_treatments))
        outputs.append(expand(PATHS.rnaseh_sensitivity_summary("{sample}"), sample=rnaseh_treatments))
    if ctx.mode == "pe" and ctx.enable_insert_size_qc:
        outputs.append(expand(PATHS.insert_size_metrics("{sample}"), sample=ctx.samples))
    if ctx.enable_multiqc:
        outputs.append(PATHS.multiqc_report())
    return outputs


def raw_read1(wildcards):
    return CTX.fastqs[wildcards.sample]["r1"]


def raw_read2(wildcards):
    return CTX.fastqs[wildcards.sample]["r2"]


def clean_reads(wildcards):
    if MODE == "pe":
        return [PATHS.clean_r1(wildcards.sample, MODE), PATHS.clean_r2(wildcards.sample)]
    return [PATHS.clean_r1(wildcards.sample, MODE)]


def bowtie2_reads_arg(wildcards, input):
    reads = [shlex.quote(str(read)) for read in input.reads]
    if MODE == "pe":
        return f"-1 {reads[0]} -2 {reads[1]}"
    return f"-U {reads[0]}"


def bowtie2_read_group_arg(wildcards):
    sample = str(wildcards.sample)
    return " ".join(["--rg-id", shlex.quote(sample), "--rg", shlex.quote(f"SM:{sample}"), "--rg", "PL:ILLUMINA"])


def macs3_format():
    return "BAMPE" if MODE == "pe" else "BAM"


def macs3_extra():
    if config.get("macs3_extra"):
        return config.get("macs3_extra")
    if CTX.peak_type == "narrow":
        return "-q 0.01"
    return "--broad --broad-cutoff 0.1 -q 0.05"


def macs3_control_input(wildcards):
    control = CTX.controls.get(wildcards.sample)
    return PATHS.peak_bam(control, CTX.peak_duplicate_mode) if control else []


def macs3_control_arg(wildcards, input):
    if getattr(input, "control", None):
        return "-c " + shlex.quote(str(input.control))
    return ""


def rnaseh_control_peak(wildcards):
    return PATHS.peak(CTX.rnaseh_controls[wildcards.sample])


def peak_input_bam(wildcards):
    """Return the BAM for MACS3 peak calling based on peak_duplicate_mode."""
    return PATHS.peak_bam(wildcards.sample, CTX.peak_duplicate_mode)


def peak_input_bai(wildcards):
    return PATHS.peak_bai(wildcards.sample, CTX.peak_duplicate_mode)


def rnaseh_signal_bigwig(sample):
    if CTX.enable_common_scale_bigwig:
        return PATHS.common_scale_bigwig_track(sample)
    return PATHS.bigwig_track(sample)


def rnaseh_signal_treatment_bigwig(wildcards):
    return rnaseh_signal_bigwig(wildcards.sample)


def rnaseh_signal_control_bigwig(wildcards):
    return rnaseh_signal_bigwig(CTX.rnaseh_controls[wildcards.sample])


def rnaseh_signal_track_type():
    return "common_scale" if CTX.enable_common_scale_bigwig else "CPM_warning"


def scale_factor_from_table_shell(sample_expr, table_expr):
    return "$(awk -F '\t' -v sample=" + sample_expr + " 'NR==1 {next} $1==sample {print $3; found=1; exit} END {if (!found) exit 2}' " + table_expr + ")"


def common_scale_base_args():
    extra = str(config.get("common_scale_bam_coverage_extra", "--binSize 10") or "--binSize 10").strip()
    args = f"{extra} --normalizeUsing {CTX.common_scale_normalization}"
    if CTX.common_scale_normalization == "RPGC":
        effective_genome_size = str(config.get("effective_genome_size", config.get("gsize", "")) or "").strip()
        if not effective_genome_size:
            raise ValueError("common_scale_normalization: RPGC requires config['effective_genome_size'] or config['gsize'].")
        args += f" --effectiveGenomeSize {shlex.quote(effective_genome_size)}"
    return args


def consensus_min_support(wildcards):
    return CTX.consensus_min_support


def replicates_for_group(wildcards):
    replicates = CTX.group_data[wildcards.group]
    if len(replicates) < 2:
        raise ValueError(f"Replicate group '{wildcards.group}' needs at least two samples.")
    return replicates


def replicate_peak_inputs(wildcards):
    return [PATHS.final_peak(sample) for sample in replicates_for_group(wildcards)]
