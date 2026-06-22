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
ALLOWED_SCALE_METHODS = {"cpm", "absolute_spikein", "matched_ref_spikein"}
ALLOWED_SPIKEIN_MODES = {"mapped_reads", "proper_pair_fragments"}
ALLOWED_RNASEH_MODES = {"ratio", "deseq2", "both"}


class WorkflowPaths:
    def __init__(self, outdir):
        outdir = str(outdir or "results").strip() or "results"
        self.outdir = outdir if outdir == "/" else outdir.rstrip("/")
        self.warning_dir = str(config.get("warning_dir", "warnings") or "warnings").rstrip("/")

        self.clean_data = self.result("clean_data")
        self.aligned_data = self.result("aligned_data")
        self.fastp_qc = self.result("qc", "fastp")
        self.mark_duplicates_qc = self.result("qc", "mark_duplicates")
        self.insert_size_qc = self.result("qc", "insert_size")
        self.frip_qc = self.result("qc", "frip")
        self.bigwig_qc = self.result("qc", "bigwig")
        self.normalization_qc = self.result("qc", "normalization")
        self.rnaseh_qc = self.result("qc")
        self.bigwig = self.result("bigwig")
        self.peaks = self.result("peaks")
        self.rnaseh = self.result("rnaseh")
        self.intersect_peaks = self.result("intersect_peaks")
        self.consensus_peaks = self.result("consensus_peaks")
        self.counts = self.result("counts")
        self.reports = self.result("reports")
        self.logs = self.result("logs")

    def result(self, *parts):
        return str(Path(self.outdir, *[str(part) for part in parts if str(part)]))

    def warning(self, *parts):
        return str(Path(self.warning_dir, *[str(part) for part in parts if str(part)]))

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
        return f"{self.aligned_data}/{sample}.aligned.filtered.bam"

    def alignment_filtered_bai(self, sample):
        return f"{self.aligned_data}/{sample}.aligned.filtered.bam.bai"

    def markdup_bam(self, sample):
        return f"{self.aligned_data}/{sample}.markdup.bam"

    def markdup_bai(self, sample):
        return f"{self.aligned_data}/{sample}.markdup.bam.bai"

    def dedup_bam(self, sample):
        return f"{self.aligned_data}/{sample}.dedup.bam"

    def dedup_bai(self, sample):
        return f"{self.aligned_data}/{sample}.dedup.bam.bai"

    def signal_bam(self, sample):
        return f"{self.aligned_data}/{sample}.signal.keepdup.filtered.bam"

    def signal_bai(self, sample):
        return f"{self.aligned_data}/{sample}.signal.keepdup.filtered.bam.bai"

    def filtered_bam(self, sample):
        return f"{self.aligned_data}/{sample}.peaks.dedup.filtered.bam"

    def filtered_bai(self, sample):
        return f"{self.aligned_data}/{sample}.peaks.dedup.filtered.bam.bai"

    def peak_bam(self, sample, peak_duplicate_mode):
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

    def bigwig_track(self, method, sample):
        return f"{self.bigwig}/{method}/{sample}.{method}.bw"

    def bigwig_header_qc(self, sample):
        return f"{self.bigwig_qc}/{sample}.bigwig_header.tsv"

    def normalization_metrics(self):
        return f"{self.normalization_qc}/spikein_summary.tsv"

    def spikein_warning_tsv(self):
        return self.warning("cut_tag_rloop_spikein.warning.tsv")

    def spikein_warning_txt(self):
        return self.warning("cut_tag_rloop_spikein.warning.txt")

    def peak(self, sample):
        suffix = "narrowPeak" if CTX.peak_type == "narrow" else "broadPeak"
        return f"{self.peaks}/{sample}.cut_tag_rloop_peaks.{suffix}"

    def peak_universe(self):
        return f"{self.counts}/peak_universe.bed"

    def peak_counts_raw(self):
        return f"{self.counts}/peak_counts.raw.tsv"

    def peak_counts_featurecounts(self):
        return f"{self.counts}/peak_counts.featureCounts.txt"

    def peak_counts_cpm(self):
        return f"{self.counts}/peak_counts.cpm.tsv"

    def peak_counts_spikein_normalized(self):
        return f"{self.counts}/peak_counts.spikein_normalized.tsv"

    def peak_annotation_input(self):
        return f"{self.counts}/peak_annotation_input.bed"

    def peak_universe_saf(self):
        return f"{self.counts}/peak_universe.saf"

    def rnaseh_ratio_table(self, sample):
        return f"{self.rnaseh}/{sample}.rnaseh_sensitive_ratio.tsv"

    def rnaseh_ratio_bed(self, sample):
        return f"{self.rnaseh}/{sample}.rnaseh_sensitive_ratio.bed"

    def rnaseh_ratio_summary(self, sample):
        return f"{self.rnaseh}/{sample}.rnaseh_sensitive_ratio.summary.tsv"

    def rnaseh_deseq2_tsv(self, contrast):
        return f"{self.rnaseh}/{contrast}.rnaseh_depleted.deseq2.tsv"

    def rnaseh_deseq2_bed(self, contrast):
        return f"{self.rnaseh}/{contrast}.rnaseh_depleted.bed"

    def intersect_peak(self, group):
        return f"{self.intersect_peaks}/{group}.intersect_peaks.bed"

    def consensus_peak(self, group):
        return f"{self.consensus_peaks}/{group}.consensus_peaks.bed"

    def replicate_correlation(self):
        return f"{self.rnaseh_qc}/replicate_correlation.tsv"

    def replicate_correlation_heatmap(self):
        return f"{self.rnaseh_qc}/replicate_correlation_heatmap.pdf"

    def sample_pca(self):
        return f"{self.rnaseh_qc}/sample_pca.pdf"

    def peak_width_distribution(self):
        return f"{self.rnaseh_qc}/peak_width_distribution.pdf"

    def rnaseh_depletion_summary(self):
        return f"{self.rnaseh_qc}/rnaseh_depletion_summary.tsv"

    def rnaseh_specificity_summary(self):
        return f"{self.rnaseh_qc}/rnaseh_specificity_summary.tsv"

    def rnaseh_signal_scatter(self):
        return f"{self.rnaseh_qc}/rnaseh_signal_scatter.pdf"

    def rnaseh_depletion_fraction_plot(self):
        return f"{self.rnaseh_qc}/rnaseh_depletion_fraction.pdf"

    def spikein_summary(self):
        return f"{self.rnaseh_qc}/spikein_summary.tsv"

    def blacklist_mito_summary(self):
        return f"{self.rnaseh_qc}/blacklist_mito_summary.tsv"

    def frip_fragment_level(self):
        return f"{self.rnaseh_qc}/frip_fragment_level.tsv"

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
    return sorted(group for group, members in group_data.items() if len(members) >= 2 and all(sample in sample_set for sample in members))


def assign_spikein_sample(mapping, sample, group, anchor, source):
    existing = mapping.get(sample)
    if existing:
        # A shared control (e.g. IgG used by multiple treatments) is allowed
        # to appear in multiple spike-in groups. Keep the first assignment so
        # the normalization script sees exactly one group per sample.
        return
    mapping[sample] = {"group": group, "anchor": anchor}


def build_spikein_group_data(sample_data):
    # Each treatment sample is its own spike-in normalization unit.
    # Multiple treatments may share the same sample_config group (for consensus peaks);
    # spike-in groups are independent: one treatment per spike-in group.
    sample_to_group = {}
    group_to_anchor = {}
    for treatment in sample_data.treatments:
        group = f"spikein_{treatment}"
        group_to_anchor[group] = treatment
        assign_spikein_sample(sample_to_group, treatment, group, treatment, "treatment")
        rnaseh = sample_data.rnaseh_controls.get(treatment)
        if rnaseh:
            assign_spikein_sample(sample_to_group, rnaseh, group, treatment, "rnaseh_control")
        control = sample_data.controls.get(treatment)
        if control:
            assign_spikein_sample(sample_to_group, control, group, treatment, "control")

    missing = [sample for sample in sample_data.samples if sample not in sample_to_group]
    if missing:
        raise ValueError(
            "Sample(s) cannot be assigned to a spike-in normalization group: "
            + ", ".join(missing)
            + ". Each non-treatment sample must be linked as control or rnaseh_control to a grouped treatment sample."
        )

    return SimpleNamespace(
        sample_to_group=sample_to_group,
        group_to_anchor=group_to_anchor,
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


def section(name):
    value = config.get(name, {}) or {}
    if not isinstance(value, dict):
        raise ValueError(f"config['{name}'] must be a mapping.")
    return value


def section_bool(section_name, key, default=False):
    value = section(section_name).get(key, default)
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    if isinstance(value, int):
        if value in {0, 1}:
            return bool(value)
        raise ValueError(f"config['{section_name}']['{key}'] must be a boolean value.")
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "t", "yes", "y", "1", "on"}:
            return True
        if normalized in {"false", "f", "no", "n", "0", "off", ""}:
            return False
    raise ValueError(f"config['{section_name}']['{key}'] must be a boolean value.")


def nested_value(section_name, key, default=None):
    return section(section_name).get(key, default)


def parse_nested_positive_int(section_name, key, default):
    try:
        value = int(nested_value(section_name, key, default))
    except (TypeError, ValueError):
        raise ValueError(f"config['{section_name}']['{key}'] must be an integer.")
    if value < 1:
        raise ValueError(f"config['{section_name}']['{key}'] must be >= 1.")
    return value


def parse_nested_nonnegative_float(section_name, key, default):
    try:
        value = float(nested_value(section_name, key, default))
    except (TypeError, ValueError):
        raise ValueError(f"config['{section_name}']['{key}'] must be a number.")
    if value < 0:
        raise ValueError(f"config['{section_name}']['{key}'] must be >= 0.")
    return value


def parse_nested_positive_float(section_name, key, default):
    value = parse_nested_nonnegative_float(section_name, key, default)
    if value <= 0:
        raise ValueError(f"config['{section_name}']['{key}'] must be > 0.")
    return value


def normalize_scale_methods(value):
    raw_methods = as_list(value if value is not None else ["cpm", "absolute_spikein", "matched_ref_spikein"])
    if not raw_methods:
        raise ValueError("config['scale_methods'] must include at least one method: cpm, absolute_spikein, or matched_ref_spikein.")
    normalized = []
    for method in raw_methods:
        low = str(method).strip().lower()
        if low == "cpm":
            normalized.append("cpm")
        elif low == "matched_ref_spikein":
            normalized.append("matched_ref_spikein")
        elif low == "absolute_spikein":
            normalized.append("absolute_spikein")
        else:
            raise ValueError(
                "config['scale_methods'] only supports cpm, absolute_spikein, and matched_ref_spikein; invalid value: "
                + str(method)
            )
    return unique_list(normalized)


def normalize_duplicate_mode(key, default, allowed):
    mode = str(config.get(key, default) or default).strip().lower()
    validate_choices([mode], allowed, key)
    return mode


def validate_mark_duplicates_extra():
    extra = str(config.get("mark_duplicates_extra", "") or "")
    if "REMOVE_DUPLICATES" in extra.upper():
        raise ValueError("config['mark_duplicates_extra'] must not contain REMOVE_DUPLICATES; use peak_duplicate_mode instead.")
    return extra


def build_workflow_context(config):
    mode = str(config.get("mode", "pe")).lower()
    validate_choices([mode], {"pe", "se"}, "mode")
    peak_type = str(config.get("peak_type", "broad")).lower()
    validate_choices([peak_type], {"broad", "narrow"}, "peak_type")
    scale_methods = normalize_scale_methods(config.get("scale_methods", ["cpm", "absolute_spikein", "matched_ref_spikein"]))

    sample_metadata = load_yaml(config.get("sample_config", "config/samples.yml"))
    sample_data = load_sample_data(sample_metadata)
    if not sample_data.treatments:
        raise ValueError("sample_config must contain at least one role: treatment sample.")
    fastqs = load_fastq_index(config, mode)
    validate_sample_fastqs(sample_data.samples, fastqs, mode)

    spikein_cfg = section("spikein")
    needs_spikein = bool({"absolute_spikein", "matched_ref_spikein"} & set(scale_methods))
    spikein_enabled = section_bool("spikein", "enabled", needs_spikein)
    spikein_index = str(spikein_cfg.get("bowtie2_index", "") or "").strip()
    spikein_genome_name = str(spikein_cfg.get("genome", "ecoli") or "ecoli").strip()
    spikein_mode = str(spikein_cfg.get("mode", "proper_pair_fragments" if mode == "pe" else "mapped_reads") or "").strip()
    validate_choices([spikein_mode], ALLOWED_SPIKEIN_MODES, "spikein.mode")
    if needs_spikein and not spikein_enabled:
        raise ValueError("scale_methods includes spike-in methods, but config['spikein']['enabled'] is false.")
    if spikein_enabled:
        if not spikein_index:
            raise ValueError("config['spikein']['bowtie2_index'] is required when spike-in normalization is enabled.")
        spikein_group_data = build_spikein_group_data(sample_data)
    else:
        spikein_group_data = SimpleNamespace(sample_to_group={}, group_to_anchor={})

    count_matrix_cfg = section("count_matrix")
    count_matrix_enabled = section_bool("count_matrix", "enabled", True)
    count_matrix_universe = str(count_matrix_cfg.get("peak_universe", "consensus") or "consensus").strip().lower()
    validate_choices([count_matrix_universe], {"union", "consensus"}, "count_matrix.peak_universe")

    rnaseh_cfg = section("rnaseh_sensitive")
    rnaseh_enabled = section_bool("rnaseh_sensitive", "enabled", True)
    rnaseh_mode = str(rnaseh_cfg.get("mode", "both") or "both").strip().lower()
    validate_choices([rnaseh_mode], ALLOWED_RNASEH_MODES, "rnaseh_sensitive.mode")

    return SimpleNamespace(
        mode=mode,
        peak_type=peak_type,
        assay=str(config.get("assay", "CUT&Tag-R-loop")),
        samples=sample_data.samples,
        sample_pattern=sample_pattern(sample_data.samples),
        treatments=sample_data.treatments,
        treatment_pattern=sample_pattern(sample_data.treatments),
        fastqs=fastqs,
        sample_records=sample_data.records,
        group_data=sample_data.groups,
        groups=valid_groups(sample_data.groups, sample_data.treatments),
        group_pattern=sample_pattern(valid_groups(sample_data.groups, sample_data.treatments)) if valid_groups(sample_data.groups, sample_data.treatments) else r"__no_valid_groups__",
        controls=sample_data.controls,
        rnaseh_controls=sample_data.rnaseh_controls,
        scale_methods=scale_methods,
        scale_method_pattern=sample_pattern(scale_methods),
        has_spikein=spikein_enabled,
        spikein_genome=spikein_genome_name,
        spikein_index=spikein_index,
        spikein_counting_mode=spikein_mode,
        spikein_sample_to_group=spikein_group_data.sample_to_group,
        spikein_group_to_anchor=spikein_group_data.group_to_anchor,
        spikein_min_mapped_reads=parse_nested_positive_int("spikein", "min_spikein_reads", 1000),
        spikein_min_fraction=parse_nested_nonnegative_float("spikein", "warn_low_fraction", 0.001),
        count_matrix_enabled=count_matrix_enabled,
        count_matrix_peak_universe=count_matrix_universe,
        count_matrix_min_peak_width=parse_nested_positive_int("count_matrix", "min_peak_width", 50),
        count_matrix_merge_distance=parse_nested_nonnegative_float("count_matrix", "merge_distance", 100),
        count_matrix_min_mapq=parse_nested_nonnegative_float("count_matrix", "min_mapq", 30),
        count_matrix_featurecounts_extra=str(count_matrix_cfg.get("featurecounts_extra", "") or ""),
        rnaseh_enabled=rnaseh_enabled,
        rnaseh_mode=rnaseh_mode,
        rnaseh_pseudocount=parse_nested_nonnegative_float("rnaseh_sensitive", "pseudocount", 0.1),
        rnaseh_signal_min_fold_change=parse_nested_positive_float("rnaseh_sensitive", "min_fold_change", 2.0),
        rnaseh_signal_min_treatment_signal=parse_nested_nonnegative_float("rnaseh_sensitive", "min_treatment_signal", 0.5),
        rnaseh_min_abs_signal_diff=parse_nested_nonnegative_float("rnaseh_sensitive", "min_abs_signal_diff", 0.2),
        rnaseh_fdr_threshold=parse_nested_positive_float("rnaseh_sensitive", "fdr_threshold", 0.05),
        rnaseh_log2fc_threshold=parse_nested_nonnegative_float("rnaseh_sensitive", "log2fc_threshold", 1.0),
        consensus_min_support=parse_positive_int_config("consensus_min_support", 2),
        mark_duplicates_extra=validate_mark_duplicates_extra(),
        peak_duplicate_mode=normalize_duplicate_mode("peak_duplicate_mode", "auto", {"remove", "keep_marked", "auto"}),
        enable_insert_size_qc=config_bool("enable_insert_size_qc", mode == "pe"),
        enable_multiqc=config_bool("enable_multiqc", True),
    )


def print_workflow_summary(ctx):
    print("Assay: " + ctx.assay)
    print("Mode: " + ctx.mode)
    print("Peak type: " + ctx.peak_type)
    print("Outdir: " + PATHS.outdir)
    print("Warning dir: " + PATHS.warning_dir)
    print("Samples: " + ",".join(ctx.samples))
    print("Treatments: " + ",".join(ctx.treatments))
    print("Replicate peak groups: " + (",".join(ctx.groups) if ctx.groups else "none"))
    print("Scale methods: " + ",".join(ctx.scale_methods))
    print("Peak duplicate mode: " + ctx.peak_duplicate_mode)
    if ctx.has_spikein:
        print("Spike-in genome: " + ctx.spikein_genome)
        print("Spike-in counting mode: " + ctx.spikein_counting_mode)
        anchors = [f"{group}:{anchor}" for group, anchor in sorted(ctx.spikein_group_to_anchor.items())]
        print("Spike-in anchor groups: " + ",".join(anchors))
    print("Count matrix: " + ("enabled" if ctx.count_matrix_enabled else "disabled"))
    print("RNaseH-sensitive analysis: " + (ctx.rnaseh_mode if ctx.rnaseh_enabled else "disabled"))
    print("MultiQC: " + ("enabled" if ctx.enable_multiqc else "disabled"))


def rnaseh_treatments(ctx):
    return [sample for sample in ctx.treatments if ctx.rnaseh_controls.get(sample)]


def rnaseh_contrasts(ctx):
    contrasts = []
    for group, members in ctx.group_data.items():
        paired = [sample for sample in members if sample in ctx.treatments and ctx.rnaseh_controls.get(sample)]
        if len(paired) >= 2:
            contrasts.append(group)
    return sorted(contrasts)


def all_outputs(ctx):
    rnaseh_samples = rnaseh_treatments(ctx)
    outputs = [
        expand(PATHS.signal_bam("{sample}"), sample=ctx.samples),
        expand(PATHS.signal_bai("{sample}"), sample=ctx.samples),
        expand(PATHS.peak("{sample}"), sample=ctx.treatments),
        expand(PATHS.frip("{sample}"), sample=ctx.treatments),
        expand(PATHS.intersect_peak("{group}"), group=ctx.groups),
        expand(PATHS.consensus_peak("{group}"), group=ctx.groups),
    ]
    if ctx.count_matrix_enabled:
        outputs.extend([
            PATHS.peak_universe(),
            PATHS.peak_counts_raw(),
            PATHS.peak_counts_featurecounts(),
            PATHS.peak_counts_cpm(),
            PATHS.peak_annotation_input(),
            PATHS.peak_universe_saf(),
        ])
        if ctx.has_spikein:
            outputs.append(PATHS.peak_counts_spikein_normalized())
    if ctx.rnaseh_enabled and ctx.rnaseh_mode in {"ratio", "both"}:
        outputs.extend([
            expand(PATHS.rnaseh_ratio_table("{sample}"), sample=rnaseh_samples),
            expand(PATHS.rnaseh_ratio_bed("{sample}"), sample=rnaseh_samples),
            expand(PATHS.rnaseh_ratio_summary("{sample}"), sample=rnaseh_samples),
        ])
    if ctx.rnaseh_enabled and ctx.rnaseh_mode in {"deseq2", "both"} and ctx.count_matrix_enabled:
        outputs.extend([
            expand(PATHS.rnaseh_deseq2_tsv("{contrast}"), contrast=rnaseh_contrasts(ctx)),
            expand(PATHS.rnaseh_deseq2_bed("{contrast}"), contrast=rnaseh_contrasts(ctx)),
        ])
    if ctx.count_matrix_enabled:
        outputs.extend([
            PATHS.replicate_correlation(),
            PATHS.replicate_correlation_heatmap(),
            PATHS.sample_pca(),
            PATHS.peak_width_distribution(),
            PATHS.rnaseh_depletion_summary(),
            PATHS.rnaseh_specificity_summary(),
            PATHS.rnaseh_signal_scatter(),
            PATHS.rnaseh_depletion_fraction_plot(),
            PATHS.spikein_summary(),
            PATHS.blacklist_mito_summary(),
            PATHS.frip_fragment_level(),
        ])
    if "cpm" in ctx.scale_methods:
        outputs.append(expand(PATHS.bigwig_track("cpm", "{sample}"), sample=ctx.samples))
        outputs.append(expand(PATHS.bigwig_header_qc("{sample}"), sample=ctx.samples))
    if ctx.has_spikein:
        outputs.append(PATHS.normalization_metrics())
        outputs.append(PATHS.spikein_warning_tsv())
        outputs.append(PATHS.spikein_warning_txt())
        if "matched_ref_spikein" in ctx.scale_methods:
            outputs.append(expand(PATHS.bigwig_track("matched_ref_spikein", "{sample}"), sample=ctx.samples))
        if "absolute_spikein" in ctx.scale_methods:
            outputs.append(expand(PATHS.bigwig_track("absolute_spikein", "{sample}"), sample=ctx.samples))
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


def spikein_bowtie2_extra():
    if config.get("spikein_bowtie2_extra"):
        return config.get("spikein_bowtie2_extra")
    if MODE == "pe":
        return "--end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant -I 10 -X 700"
    return "--end-to-end --very-sensitive"


def macs3_control_input(wildcards):
    control = CTX.controls.get(wildcards.sample)
    return PATHS.peak_bam(control, CTX.peak_duplicate_mode) if control else []


def macs3_control_arg(wildcards, input):
    if getattr(input, "control", None):
        return "-c " + shlex.quote(str(input.control))
    return ""


def peak_input_bam(wildcards):
    return PATHS.peak_bam(wildcards.sample, CTX.peak_duplicate_mode)


def peak_input_bai(wildcards):
    return PATHS.peak_bai(wildcards.sample, CTX.peak_duplicate_mode)


def scaled_bigwig_input(wildcards):
    return PATHS.bigwig_track(wildcards.scale_method, wildcards.sample)


def rnaseh_sensitive_treatment_bigwig(wildcards):
    method = "matched_ref_spikein" if CTX.has_spikein and "matched_ref_spikein" in CTX.scale_methods else "cpm"
    return PATHS.bigwig_track(method, wildcards.sample)


def rnaseh_sensitive_control_bigwig(wildcards):
    method = "matched_ref_spikein" if CTX.has_spikein and "matched_ref_spikein" in CTX.scale_methods else "cpm"
    return PATHS.bigwig_track(method, CTX.rnaseh_controls[wildcards.sample])


def consensus_min_support(wildcards):
    return CTX.consensus_min_support


def replicates_for_group(wildcards):
    replicates = CTX.group_data[wildcards.group]
    if len(replicates) < 2:
        raise ValueError(f"Replicate group '{wildcards.group}' needs at least two samples.")
    return replicates


def replicate_peak_inputs(wildcards):
    return [PATHS.peak(sample) for sample in replicates_for_group(wildcards)]


def replicate_sensitive_inputs(wildcards):
    return [PATHS.rnaseh_ratio_bed(sample) for sample in replicates_for_group(wildcards) if CTX.rnaseh_controls.get(sample)]


def count_matrix_peak_inputs(wildcards):
    return [PATHS.peak(sample) for sample in CTX.treatments]


def count_matrix_bam_inputs(wildcards):
    return [PATHS.signal_bam(sample) for sample in CTX.samples]


def rnaseh_contrast_pairs_lines(wildcards):
    lines = ["sample\tcondition\tpair"]
    for treatment in CTX.group_data[wildcards.contrast]:
        if treatment not in CTX.treatments:
            continue
        rnaseh = CTX.rnaseh_controls.get(treatment)
        if not rnaseh:
            continue
        lines.append("\t".join([treatment, "NoRNaseH", treatment]))
        lines.append("\t".join([rnaseh, "RNaseH", treatment]))
    return "\n".join(lines)


def fastq_metrics_lines():
    lines = ["sample\tread1\tread2"]
    for sample in CTX.samples:
        reads = CTX.fastqs[sample]
        lines.append("\t".join([sample, reads["r1"], reads.get("r2") or ""]))
    return "\n".join(lines)


def spikein_group_lines():
    lines = ["sample\tspikein_group\tspikein_anchor_sample"]
    for sample in CTX.samples:
        assignment = CTX.spikein_sample_to_group[sample]
        lines.append("\t".join([sample, assignment["group"], assignment["anchor"]]))
    return "\n".join(lines)
