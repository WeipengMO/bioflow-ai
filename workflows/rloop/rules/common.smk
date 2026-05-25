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
        self.bigwig = self.result("bigwig")
        self.macs3_results = self.result("macs3_results")
        self.rnaseh_filtered = self.result("rnaseh_filtered")
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

    def dedup_bam(self, sample):
        return f"{self.aligned_data}/{sample}.sorted.rmdup.bam"

    def dedup_bai(self, sample):
        return f"{self.aligned_data}/{sample}.sorted.rmdup.bam.bai"

    def filtered_bam(self, sample):
        return f"{self.aligned_data}/{sample}.sorted.rmdup.filtered.bam"

    def filtered_bai(self, sample):
        return f"{self.aligned_data}/{sample}.sorted.rmdup.filtered.bam.bai"

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
        return f"{self.bigwig}/{sample}.sorted.rmdup.filtered.CPM.bw"

    def peak(self, sample):
        peak_type = str(config.get("peak_type", "broad")).lower()
        if peak_type == "narrow":
            return f"{self.macs3_results}/narrow/{sample}_peaks.narrowPeak"
        return f"{self.macs3_results}/broad/{sample}_peaks.broadPeak"

    def rnaseh_sensitive_peak(self, sample):
        peak_type = str(config.get("peak_type", "broad")).lower()
        suffix = "narrowPeak" if peak_type == "narrow" else "broadPeak"
        return f"{self.rnaseh_filtered}/{sample}_rnaseh_sensitive.{suffix}"

    def final_peak(self, sample):
        if CTX.enable_rnaseh_subtraction and CTX.rnaseh_controls.get(sample):
            return self.rnaseh_sensitive_peak(sample)
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
    print("RNase H subtraction: " + ("enabled" if ctx.enable_rnaseh_subtraction else "disabled"))
    print("MultiQC: " + ("enabled" if ctx.enable_multiqc else "disabled"))


def all_outputs(ctx):
    outputs = [
        expand(PATHS.filtered_bam("{sample}"), sample=ctx.samples),
        expand(PATHS.filtered_bai("{sample}"), sample=ctx.samples),
        expand(PATHS.bigwig_track("{sample}"), sample=ctx.samples),
        expand(PATHS.peak("{sample}"), sample=ctx.treatments),
        expand(PATHS.frip("{sample}"), sample=ctx.treatments),
        expand(PATHS.replicate_consensus_peak("{group}"), group=ctx.groups),
    ]
    if ctx.enable_rnaseh_subtraction:
        rnaseh_treatments = [sample for sample in ctx.treatments if ctx.rnaseh_controls.get(sample)]
        outputs.append(expand(PATHS.rnaseh_sensitive_peak("{sample}"), sample=rnaseh_treatments))
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
        return "--keep-dup all -q 0.01"
    return "--broad --broad-cutoff 0.1 --keep-dup all -q 0.05"


def macs3_control_input(wildcards):
    control = CTX.controls.get(wildcards.sample)
    return PATHS.filtered_bam(control) if control else []


def macs3_control_arg(wildcards, input):
    if getattr(input, "control", None):
        return "-c " + shlex.quote(str(input.control))
    return ""


def rnaseh_control_peak(wildcards):
    return PATHS.peak(CTX.rnaseh_controls[wildcards.sample])


def consensus_min_support(wildcards):
    return CTX.consensus_min_support


def replicates_for_group(wildcards):
    replicates = CTX.group_data[wildcards.group]
    if len(replicates) < 2:
        raise ValueError(f"Replicate group '{wildcards.group}' needs at least two samples.")
    return replicates


def replicate_peak_inputs(wildcards):
    return [PATHS.final_peak(sample) for sample in replicates_for_group(wildcards)]
