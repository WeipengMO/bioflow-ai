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

PEAK_MODES = {"with_control", "without_control"}
PEAK_TYPES = {"narrow", "broad"}
PEAK_EXTENSIONS = {"narrow": "narrowPeak", "broad": "broadPeak"}
PEAK_DIRS = {
    "narrow": ("with_control", "narrow"),
    "broad": ("with_control", "broad"),
    "narrow_no_control": ("without_control", "narrow"),
    "broad_no_control": ("without_control", "broad"),
}
CHIPQC_COLUMNS = [
    "SampleID",
    "Tissue",
    "Factor",
    "Condition",
    "Treatment",
    "Replicate",
    "bamReads",
    "ControlID",
    "bamControl",
    "Peaks",
    "PeakCaller",
]


class WorkflowPaths:
    def __init__(self, outdir):
        outdir = str(outdir or "results").strip() or "results"
        self.outdir = outdir if outdir == "/" else outdir.rstrip("/")

        self.clean_data = self.result("clean_data")
        self.aligned_data = self.result("aligned_data")
        self.fastp_qc = self.result("qc", "fastp")
        self.mark_duplicates_qc = self.result("qc", "mark_duplicates")
        self.bigwig = self.result("bigwig")
        self.deeptools_profile = self.result("deeptools_profile")
        self.macs3_results = self.result("macs3_results")
        self.replicate_intersect = self.result("replicate_intersect")
        self.logs = self.result("logs")
        self.default_chipqc_report = self.result("reports", "chipqc")
        self.default_homer_report = self.result("reports", "homer")
        self.default_homer_input = self.result("homer_inputs")

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

    def bigwig_track(self, sample):
        return f"{self.bigwig}/{sample}.sorted.rmdup.CPM.bw"

    def matrix(self, sample):
        return f"{self.deeptools_profile}/{sample}.matrix.gz"

    def profile_png(self, sample):
        return f"{self.deeptools_profile}/{sample}.scale.png"

    def peak(self, sample, peak_mode, peak_type):
        directory = peak_directory(peak_mode, peak_type)
        extension = PEAK_EXTENSIONS[peak_type]
        return f"{self.macs3_results}/{directory}/{sample}_peaks.{extension}"

    def replicate_consensus(self, group):
        return f"{self.replicate_intersect}/{group}_intersect.bed"

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


def load_yaml(path):
    path = Path(path)
    if not path.exists():
        raise ValueError(f"Required YAML file does not exist: {path}")
    with path.open() as handle:
        return yaml.safe_load(handle) or {}


def validate_choices(values, allowed, key):
    invalid = sorted(set(values) - set(allowed))
    if invalid:
        raise ValueError(f"config['{key}'] has invalid value(s): {', '.join(invalid)}")


def workflow_threads(name, default):
    threads = config.get("threads", {}) or {}
    return int(threads.get(name, default))


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


def sample_pattern(samples):
    return "|".join(re.escape(sample) for sample in samples)


def load_sample_data(samples):
    if not isinstance(samples, dict):
        raise ValueError("sample_config must contain a mapping of sample name to sample metadata.")

    records = {}
    for sample, value in samples.items():
        if value is None:
            value = {}
        if not isinstance(value, dict):
            raise ValueError(f"sample_config entry '{sample}' must be a mapping.")

        sample = str(sample)
        role = str(value.get("role", "treatment")).lower()
        if role in {"treatment", "treat", "peak"}:
            role = "treatment"
        elif role in {"control", "input", "igg"}:
            role = "control"
        else:
            raise ValueError(f"sample_config entry '{sample}'.role must be treatment or control.")

        records[sample] = {
            "role": role,
            "group": value.get("group"),
            "control": value.get("control"),
            "pool": value.get("pool"),
        }

    treatments = [sample for sample, record in records.items() if record["role"] == "treatment"]
    control_samples = [sample for sample, record in records.items() if record["role"] == "control"]
    pooled_control_groups = {}
    for sample in control_samples:
        pool = records[sample]["pool"]
        if pool:
            pooled_control_groups.setdefault(str(pool), []).append(sample)

    pool_name_collisions = sorted(set(pooled_control_groups) & set(records))
    if pool_name_collisions:
        raise ValueError("Control pool name(s) must not duplicate sample names: " + ", ".join(pool_name_collisions))

    valid_controls = set(control_samples) | set(pooled_control_groups)
    control_map = {}
    for treatment in treatments:
        control = records[treatment]["control"]
        if not control:
            continue
        control = str(control)
        if control not in valid_controls:
            raise ValueError(
                f"sample_config entry '{treatment}'.control references unknown control or pool: {control}"
            )
        control_map[treatment] = control

    if control_map and set(control_map) != set(treatments):
        missing = sorted(set(treatments) - set(control_map))
        raise ValueError("Treatment sample(s) missing control: " + ", ".join(missing))

    pooled_names = set(pooled_control_groups)
    raw_control_targets = set(control_map.values()) & set(control_samples)
    pooled_targets = set(control_map.values()) & pooled_names
    if raw_control_targets and pooled_targets:
        raise ValueError("Do not mix raw controls and pooled controls in one sample_config file.")
    elif pooled_targets:
        strategy = "pooled"
    elif raw_control_targets:
        strategy = "matched"
    else:
        strategy = "none"

    return SimpleNamespace(
        treatments=unique_list(treatments),
        control_strategy=strategy,
        control_map=control_map,
        control_samples=unique_list(control_samples),
        pooled_control_groups={name: unique_list(samples) for name, samples in pooled_control_groups.items()},
        sample_groups={
            sample: str(record["group"])
            for sample, record in records.items()
            if record.get("group")
        },
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
            sample = strip_fastq_suffix(path.name)
            index[sample] = {"r1": str(path), "r2": None}
    return index


def strip_fastq_suffix(filename):
    for suffix in sorted(FASTQ_SUFFIXES, key=len, reverse=True):
        if filename.lower().endswith(suffix):
            return filename[: -len(suffix)]
    return filename


def load_fastq_index(config, mode):
    manifest = config.get("fastq_manifest", "")
    index = load_fastq_manifest(manifest, mode) if manifest else discover_fastqs(config.get("raw_data_dir", "raw_data"), mode)

    for sample, reads in sorted(index.items()):
        if not reads.get("r1"):
            raise ValueError(f"Sample '{sample}' is missing read1 FASTQ.")
        if mode == "pe" and not reads.get("r2"):
            raise ValueError(
                f"Sample '{sample}' is missing read2 FASTQ. "
                "Use standard names such as sample_R1.fastq.gz/sample_R2.fastq.gz, "
                "sample_1.fq.gz/sample_2.fq.gz, or provide config['fastq_manifest']."
            )
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


def peak_directory(peak_mode, peak_type):
    validate_choices([peak_mode], PEAK_MODES, "peak_mode")
    validate_choices([peak_type], PEAK_TYPES, "peak_type")
    suffix = "_no_control" if peak_mode == "without_control" else ""
    return f"{peak_type}{suffix}"


def parse_peak_dir(directory):
    if directory not in PEAK_DIRS:
        raise ValueError(f"Invalid MACS3 peak directory: {directory}")
    return PEAK_DIRS[directory]


def macs3_format(mode):
    return "BAMPE" if mode == "pe" else "BAM"


def select_default_peak_mode(call_peak_modes, has_controls):
    if has_controls and "with_control" in call_peak_modes:
        return "with_control"
    if "without_control" in call_peak_modes:
        return "without_control"
    return "with_control"


def validate_peak_selection(peak_mode, peak_type, call_peak_modes, call_peak_types, has_controls, label):
    validate_choices([peak_mode], PEAK_MODES, f"{label}_peak_mode")
    validate_choices([peak_type], PEAK_TYPES, f"{label}_peak_type")
    if peak_mode not in call_peak_modes:
        raise ValueError(f"config['{label}_peak_mode'] is '{peak_mode}', but call_peak_modes does not include it.")
    if peak_type not in call_peak_types:
        raise ValueError(f"config['{label}_peak_type'] is '{peak_type}', but call_peak_types does not include it.")
    if peak_mode == "with_control" and not has_controls:
        raise ValueError(f"config['{label}_peak_mode'] is 'with_control', but no effective controls are configured.")


def callpeak_outputs(treatments, call_peak_modes, call_peak_types, has_controls):
    outputs = []
    if "with_control" in call_peak_modes and has_controls:
        for peak_type in call_peak_types:
            outputs.extend(expand(PATHS.peak("{sample}", "with_control", peak_type), sample=treatments))
    if "without_control" in call_peak_modes:
        for peak_type in call_peak_types:
            outputs.extend(expand(PATHS.peak("{sample}", "without_control", peak_type), sample=treatments))
    return outputs


def peak_selection_outputs(treatments, peak_mode, peak_type):
    return expand(PATHS.peak("{sample}", peak_mode, peak_type), sample=treatments)


def analysis_bam(sample):
    return PATHS.filtered_bam(sample) if CTX.enable_blacklist_filter else PATHS.dedup_bam(sample)


def analysis_bai(sample):
    return PATHS.filtered_bai(sample) if CTX.enable_blacklist_filter else PATHS.dedup_bai(sample)


def default_alignment_filter_view_extra():
    return "-f 2 -F 1804" if MODE == "pe" else "-F 1796"


def alignment_filter_view_extra(wildcards):
    return config.get("alignment_filter_view_extra", "") or default_alignment_filter_view_extra()


def mark_duplicates_input_bam(wildcards):
    if CTX.enable_alignment_filter:
        return PATHS.alignment_filtered_bam(wildcards.sample)
    return PATHS.sorted_bam(wildcards.sample)


def mark_duplicates_input_bai(wildcards):
    if CTX.enable_alignment_filter:
        return PATHS.alignment_filtered_bai(wildcards.sample)
    return PATHS.sorted_bai(wildcards.sample)


def load_replicate_data(samples):
    replicate_data = {}
    for sample, value in samples.items():
        if not isinstance(value, dict):
            continue
        role = str(value.get("role", "treatment")).lower()
        if role not in {"treatment", "treat", "peak"}:
            continue
        group = value.get("group")
        if group:
            replicate_data.setdefault(str(group), []).append(str(sample))
    return {group: unique_list(samples) for group, samples in replicate_data.items()}


def valid_replicate_groups(replicate_data, treatments):
    treatment_set = set(treatments)
    return sorted(
        group
        for group, replicates in replicate_data.items()
        if len(replicates) >= 2 and all(sample in treatment_set for sample in replicates)
    )


def build_workflow_context(config):
    mode = str(config.get("mode", "pe")).lower()
    validate_choices([mode], {"pe", "se"}, "mode")

    sample_metadata = load_yaml(config.get("sample_config", "config/samples.yml"))
    sample_data = load_sample_data(sample_metadata)
    fastqs = load_fastq_index(config, mode)

    declared_samples = unique_list(sample_data.treatments + sample_data.control_samples)
    discovered_samples = sorted(fastqs)
    samples = sorted(set(declared_samples or discovered_samples))
    if not samples:
        raise ValueError("No samples were declared in sample_config and no FASTQs were discovered.")

    treatments = sample_data.treatments or samples
    validate_sample_fastqs(samples, fastqs, mode)

    control_map = sample_data.control_map
    has_controls = bool(control_map)
    call_peak_modes = as_list(config.get("call_peak_modes", ["with_control", "without_control"]))
    call_peak_types = as_list(config.get("call_peak_types", ["narrow", "broad"]))
    validate_choices(call_peak_modes, PEAK_MODES, "call_peak_modes")
    validate_choices(call_peak_types, PEAK_TYPES, "call_peak_types")

    default_peak_mode = select_default_peak_mode(call_peak_modes, has_controls)
    default_peak_type = "narrow" if "narrow" in call_peak_types else call_peak_types[0]

    replicate_data = load_replicate_data(sample_metadata)
    replicate_peak_mode = config.get("replicate_peak_mode", default_peak_mode)
    replicate_peak_type = config.get("replicate_peak_type", default_peak_type)
    validate_peak_selection(
        replicate_peak_mode,
        replicate_peak_type,
        call_peak_modes,
        call_peak_types,
        has_controls,
        "replicate",
    )
    replicate_groups = valid_replicate_groups(replicate_data, treatments)

    enable_deeptools_profile = config_bool("enable_deeptools_profile", bool(config.get("regions_bed")))
    if enable_deeptools_profile and not config.get("regions_bed"):
        raise ValueError("config['regions_bed'] is required when enable_deeptools_profile is true.")

    enable_alignment_filter = config_bool("enable_alignment_filter", True)
    enable_blacklist_filter = config_bool("enable_blacklist_filter", True)
    blacklist = str(config.get("blacklist", "") or "")

    enable_chipqc = config_bool("enable_chipqc", False)
    chipqc_peak_mode = config.get("chipqc_peak_mode", default_peak_mode)
    chipqc_peak_type = config.get("chipqc_peak_type", default_peak_type)
    if enable_chipqc:
        validate_peak_selection(
            chipqc_peak_mode,
            chipqc_peak_type,
            call_peak_modes,
            call_peak_types,
            has_controls,
            "chipqc",
        )

    enable_homer = config_bool("enable_homer", False)
    homer_peak_mode = config.get("homer_peak_mode", default_peak_mode)
    homer_peak_type = config.get("homer_peak_type", default_peak_type)
    homer_input_source = config.get("homer_input_source", "replicate_intersect")
    validate_choices([homer_input_source], {"sample_peaks", "replicate_intersect"}, "homer_input_source")
    if enable_homer:
        validate_peak_selection(
            homer_peak_mode,
            homer_peak_type,
            call_peak_modes,
            call_peak_types,
            has_controls,
            "homer",
        )
        if not config.get("homer_genome"):
            raise ValueError("config['homer_genome'] is required when enable_homer is true.")
        if not config.get("gtf"):
            raise ValueError("config['gtf'] is required when enable_homer is true.")

    homer_targets = replicate_groups if homer_input_source == "replicate_intersect" else treatments
    homer_targets = [target for target in homer_targets if target not in set(as_list(config.get("homer_exclude_targets", [])))]
    if enable_homer and homer_input_source == "replicate_intersect" and not homer_targets:
        raise ValueError(
            "HOMER is configured to use replicate_intersect, but no valid replicate groups were found. "
            "Check treatment sample group metadata."
        )

    return SimpleNamespace(
        mode=mode,
        samples=samples,
        sample_pattern=sample_pattern(samples),
        treatments=treatments,
        treatment_pattern=sample_pattern(treatments),
        fastqs=fastqs,
        control_strategy=sample_data.control_strategy,
        use_pooled_control=sample_data.control_strategy == "pooled",
        pooled_control_groups=sample_data.pooled_control_groups,
        pooled_control_names=sorted(sample_data.pooled_control_groups),
        pooled_control_pattern=sample_pattern(sorted(sample_data.pooled_control_groups)),
        treatment_to_control=control_map,
        effective_controls=sorted(set(control_map.values())),
        has_effective_control=has_controls,
        sample_groups=sample_data.sample_groups,
        call_peak_modes=call_peak_modes,
        call_peak_types=call_peak_types,
        replicate_data=replicate_data,
        replicate_peak_mode=replicate_peak_mode,
        replicate_peak_type=replicate_peak_type,
        replicate_groups=replicate_groups,
        enable_deeptools_profile=enable_deeptools_profile,
        enable_alignment_filter=enable_alignment_filter,
        enable_blacklist_filter=enable_blacklist_filter,
        blacklist=blacklist,
        enable_chipqc=enable_chipqc,
        chipqc_peak_mode=chipqc_peak_mode,
        chipqc_peak_type=chipqc_peak_type,
        enable_homer=enable_homer,
        homer_peak_mode=homer_peak_mode,
        homer_peak_type=homer_peak_type,
        homer_input_source=homer_input_source,
        homer_targets=homer_targets,
        homer_result_tag=f"{homer_input_source}_{homer_peak_mode}_{homer_peak_type}",
    )


def print_workflow_summary(ctx):
    print("Mode: " + ctx.mode)
    print("Outdir: " + PATHS.outdir)
    print("Raw samples: " + ",".join(ctx.samples))
    print("Peak treatments: " + ",".join(ctx.treatments))
    print("Control strategy: " + ctx.control_strategy)
    print("Peak outputs: " + ";".join([f"{m}:{t}" for m in ctx.call_peak_modes for t in ctx.call_peak_types]))
    print("Replicate groups: " + (",".join(ctx.replicate_groups) if ctx.replicate_groups else "none"))
    print("Pre-MarkDuplicates alignment filter: " + ("enabled" if ctx.enable_alignment_filter else "disabled"))
    print("Post-MarkDuplicates blacklist filter: " + ("enabled" if ctx.enable_blacklist_filter else "disabled"))
    print("ChIPQC: " + ("enabled" if ctx.enable_chipqc else "disabled"))
    print("HOMER: " + ("enabled" if ctx.enable_homer else "disabled"))


def all_outputs(ctx):
    outputs = [
        expand(analysis_bam("{sample}"), sample=ctx.samples),
        expand(analysis_bai("{sample}"), sample=ctx.samples),
        expand(PATHS.bigwig_track("{sample}"), sample=ctx.samples),
        callpeak_outputs(ctx.treatments, ctx.call_peak_modes, ctx.call_peak_types, ctx.has_effective_control),
        expand(PATHS.replicate_consensus("{group}"), group=ctx.replicate_groups),
    ]

    if ctx.enable_deeptools_profile:
        outputs.append(expand(PATHS.profile_png("{sample}"), sample=ctx.samples))
    if ctx.enable_chipqc:
        outputs.append(CHIPQC_DONE)
    if ctx.enable_homer:
        outputs.append(homer_outputs(ctx.homer_targets, ctx.homer_result_tag))
    return outputs


def raw_read1(wildcards):
    return CTX.fastqs[wildcards.sample]["r1"]


def raw_read2(wildcards):
    return CTX.fastqs[wildcards.sample]["r2"]


def clean_reads(wildcards):
    if MODE == "pe":
        return [
            PATHS.clean_r1(wildcards.sample, MODE),
            PATHS.clean_r2(wildcards.sample),
        ]
    return [PATHS.clean_r1(wildcards.sample, MODE)]


def replicates_for_group(wildcards):
    replicates = CTX.replicate_data[wildcards.group]
    if len(replicates) < 2:
        raise ValueError(f"Replicate group '{wildcards.group}' needs at least two samples.")
    return replicates


def replicate_peak_inputs(wildcards):
    return [PATHS.peak(sample, CTX.replicate_peak_mode, CTX.replicate_peak_type) for sample in replicates_for_group(wildcards)]


def macs3_input(wildcards):
    peak_mode, peak_type = parse_peak_dir(wildcards.peak_dir)
    expected_ext = PEAK_EXTENSIONS[peak_type]
    if wildcards.peak_ext != expected_ext:
        raise ValueError(
            f"Inconsistent MACS3 output extension for {wildcards.peak_dir}: "
            f"expected {expected_ext}, got {wildcards.peak_ext}."
        )

    inputs = {"treatment": analysis_bam(wildcards.treatment)}
    if peak_mode == "with_control":
        if wildcards.treatment not in CTX.treatment_to_control:
            raise ValueError(f"No control configured for treatment: {wildcards.treatment}")
        inputs["control"] = analysis_bam(CTX.treatment_to_control[wildcards.treatment])
    return inputs


def macs3_control_arg(wildcards, input):
    peak_mode, _ = parse_peak_dir(wildcards.peak_dir)
    if peak_mode == "with_control":
        return "-c " + shlex.quote(str(input.control))
    return ""


def macs3_broad_arg(wildcards):
    _, peak_type = parse_peak_dir(wildcards.peak_dir)
    if peak_type != "broad":
        return ""
    return f"--broad --broad-cutoff {config.get('macs3_broad_cutoff', 0.1)}"


def pooled_control_inputs(wildcards):
    return expand(analysis_bam("{sample}"), sample=CTX.pooled_control_groups[str(wildcards.pool)])


def bowtie2_reads_arg(wildcards, input):
    reads = [shlex.quote(str(read)) for read in input.reads]
    if MODE == "pe":
        return f"-1 {reads[0]} -2 {reads[1]}"
    return f"-U {reads[0]}"


def bowtie2_read_group_arg(wildcards):
    sample = str(wildcards.sample)
    return " ".join(
        [
            "--rg-id",
            shlex.quote(sample),
            "--rg",
            shlex.quote(f"SM:{sample}"),
            "--rg",
            "PL:ILLUMINA",
        ]
    )


def homer_raw_peak_input(wildcards):
    if CTX.homer_input_source == "replicate_intersect":
        return PATHS.replicate_consensus(wildcards.target)
    return PATHS.peak(wildcards.target, CTX.homer_peak_mode, CTX.homer_peak_type)


def homer_outputs(targets, result_tag):
    report_dir = str(config.get("homer_report_dir") or PATHS.default_homer_report)
    return (
        expand(f"{report_dir}/{result_tag}/{{target}}/.motifs_complete", target=targets)
        + expand(f"{report_dir}/{result_tag}/{{target}}/annotatePeaks.txt", target=targets)
        + expand(f"{report_dir}/{result_tag}/{{target}}/peak_feature_distribution.pie.png", target=targets)
        + expand(f"{report_dir}/{result_tag}/{{target}}/peak_feature_distribution.tsv", target=targets)
        + [
            f"{report_dir}/{result_tag}/peak_feature_distribution.summary.tsv",
            f"{report_dir}/{result_tag}/peak_feature_distribution.stacked_bar.png",
        ]
    )


def chipqc_peak_caller(peak_type):
    if peak_type == "narrow":
        return "narrow"
    if peak_type == "broad":
        return "bed"
    raise ValueError("ChIPQC peak_type must be narrow or broad.")


def chipqc_replicate_map(replicate_data):
    mapping = {}
    for samples in replicate_data.values():
        for index, sample in enumerate(samples, start=1):
            mapping[str(sample)] = str(index)
    return mapping


def chipqc_sample_sheet_rows(ctx):
    replicates = chipqc_replicate_map(ctx.replicate_data)
    rows = []
    for treatment in ctx.treatments:
        control_id = ctx.treatment_to_control.get(treatment, "NA")
        bam_control = analysis_bam(control_id) if control_id != "NA" else "NA"
        condition = ctx.sample_groups.get(treatment, "NA")
        rows.append(
            {
                "SampleID": treatment,
                "Tissue": "NA",
                "Factor": "NA",
                "Condition": condition,
                "Treatment": "NA",
                "Replicate": replicates.get(treatment, "NA"),
                "bamReads": analysis_bam(treatment),
                "ControlID": control_id,
                "bamControl": bam_control,
                "Peaks": PATHS.peak(treatment, ctx.chipqc_peak_mode, ctx.chipqc_peak_type),
                "PeakCaller": chipqc_peak_caller(ctx.chipqc_peak_type),
            }
        )
    return rows


def write_chipqc_sample_sheet(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CHIPQC_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
