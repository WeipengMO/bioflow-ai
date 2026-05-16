from pathlib import Path
import csv
import re
import yaml


FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
PAIRED_FASTQ_RE = re.compile(
    r"^(?P<sample>.+?)(?:[._-]R?|[._-])(?P<read>[12])(?:_001)?(?P<suffix>\.f(?:ast)?q(?:\.gz)?)$",
    re.IGNORECASE,
)

PEAK_DIRS = {
    "narrow": ("with_control", "narrow"),
    "broad": ("with_control", "broad"),
    "narrow_no_control": ("without_control", "narrow"),
    "broad_no_control": ("without_control", "broad"),
}
PEAK_EXTENSIONS = {"narrow": "narrowPeak", "broad": "broadPeak"}


def as_list(value):
    if value is None:
        return []
    if isinstance(value, str):
        return [item.strip() for item in value.split(",") if item.strip()]
    return list(value)


def unique_list(values):
    return list(dict.fromkeys(str(value) for value in values))


def load_yaml(path):
    path = Path(path)
    if not path.exists():
        return {}
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
    return bool(config.get(key, default))


def load_sample_data(path):
    data = load_yaml(path)
    treatments = unique_list(as_list(data.get("treatments", [])))
    controls = data.get("controls", {}) or {}
    strategy = str(controls.get("strategy", "none"))
    validate_choices([strategy], {"none", "pooled", "matched"}, "controls.strategy")

    pooled_name = str(controls.get("pooled_name", "pooled_input"))
    pooled_samples = unique_list(as_list(controls.get("pooled_samples", [])))
    matched_controls = {
        str(treatment): str(control)
        for treatment, control in (controls.get("matched", {}) or {}).items()
    }

    if strategy == "pooled":
        if not treatments:
            raise ValueError("controls.strategy is 'pooled', but treatments is empty.")
        if not pooled_samples:
            raise ValueError("controls.strategy is 'pooled', but controls.pooled_samples is empty.")
        control_map = {treatment: pooled_name for treatment in treatments}
        control_samples = pooled_samples
    elif strategy == "matched":
        if not treatments:
            raise ValueError("controls.strategy is 'matched', but treatments is empty.")
        missing = sorted(set(treatments) - set(matched_controls))
        if missing:
            raise ValueError("controls.matched is missing treatment(s): " + ", ".join(missing))
        control_map = {treatment: matched_controls[treatment] for treatment in treatments}
        control_samples = unique_list(control_map.values())
    else:
        control_map = {}
        control_samples = []

    return {
        "treatments": treatments,
        "control_strategy": strategy,
        "control_map": control_map,
        "control_samples": control_samples,
        "pooled_control_name": pooled_name,
        "pooled_control_samples": pooled_samples,
    }


def _empty_fastq_record():
    return {"r1": None, "r2": None}


def _set_fastq_record(index, sample, read, path):
    record = index.setdefault(sample, _empty_fastq_record())
    key = f"r{read}"
    if record.get(key) is not None:
        raise ValueError(
            f"Duplicate FASTQ for sample '{sample}' read {read}: {record[key]} and {path}. "
            "Use config['fastq_manifest'] to disambiguate irregular file names."
        )
    record[key] = str(path)


def _load_fastq_manifest(path, mode):
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


def _discover_fastqs(raw_dir, mode):
    raw_dir = Path(raw_dir)
    index = {}
    if not raw_dir.exists():
        return index

    for path in sorted(p for p in raw_dir.iterdir() if p.is_file()):
        lower_name = path.name.lower()
        if not lower_name.endswith(FASTQ_SUFFIXES):
            continue

        if mode == "pe":
            match = PAIRED_FASTQ_RE.match(path.name)
            if not match:
                continue
            sample = match.group("sample")
            read = match.group("read")
            _set_fastq_record(index, sample, read, path)
        else:
            sample = path.name
            for suffix in sorted(FASTQ_SUFFIXES, key=len, reverse=True):
                if sample.lower().endswith(suffix):
                    sample = sample[: -len(suffix)]
                    break
            index[sample] = {"r1": str(path), "r2": None}
    return index


def load_fastq_index(config, mode):
    manifest = config.get("fastq_manifest", "")
    if manifest:
        index = _load_fastq_manifest(manifest, mode)
    else:
        index = _discover_fastqs(config.get("raw_data_dir", "raw_data"), mode)

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
            + ". Check raw_data_dir, file naming, samples.yml, or fastq_manifest."
        )
    if mode == "pe":
        missing_r2 = sorted(sample for sample in samples if not fastqs[sample].get("r2"))
        if missing_r2:
            raise ValueError("Missing read2 FASTQ for sample(s): " + ", ".join(missing_r2))


def raw_read1(wildcards):
    return FASTQS[wildcards.sample]["r1"]


def raw_read2(wildcards):
    return FASTQS[wildcards.sample]["r2"]


def clean_read1(wildcards):
    if MODE == "pe":
        return f"clean_data/{wildcards.sample}.R1.clean.fq.gz"
    return f"clean_data/{wildcards.sample}.clean.fq.gz"


def clean_read2(wildcards):
    return f"clean_data/{wildcards.sample}.R2.clean.fq.gz"


def macs3_format():
    return "BAMPE" if MODE == "pe" else "BAM"


def peak_dir(peak_mode, peak_type):
    validate_choices([peak_mode], {"with_control", "without_control"}, "peak_mode")
    validate_choices([peak_type], {"narrow", "broad"}, "peak_type")
    suffix = "_no_control" if peak_mode == "without_control" else ""
    return f"{peak_type}{suffix}"


def peak_path(sample, peak_mode, peak_type):
    directory = peak_dir(peak_mode, peak_type)
    extension = PEAK_EXTENSIONS[peak_type]
    return f"macs3_results/{directory}/{sample}_peaks.{extension}"


def parse_peak_dir(directory):
    if directory not in PEAK_DIRS:
        raise ValueError(f"Invalid MACS3 peak directory: {directory}")
    return PEAK_DIRS[directory]


def expected_peak_extension(directory):
    _, peak_type = parse_peak_dir(directory)
    return PEAK_EXTENSIONS[peak_type]


def callpeak_outputs(treatments, call_peak_modes, call_peak_types, has_with_control):
    outputs = []
    if "with_control" in call_peak_modes and has_with_control:
        for peak_type in call_peak_types:
            outputs.extend(expand(peak_path("{sample}", "with_control", peak_type), sample=treatments))
    if "without_control" in call_peak_modes:
        for peak_type in call_peak_types:
            outputs.extend(expand(peak_path("{sample}", "without_control", peak_type), sample=treatments))
    return outputs


def peak_selection_outputs(treatments, peak_mode, peak_type):
    return expand(peak_path("{sample}", peak_mode, peak_type), sample=treatments)


def select_default_peak_mode(call_peak_modes, has_with_control):
    if has_with_control and "with_control" in call_peak_modes:
        return "with_control"
    if "without_control" in call_peak_modes:
        return "without_control"
    return "with_control"


def validate_peak_selection(peak_mode, peak_type, call_peak_modes, call_peak_types, has_with_control, label):
    validate_choices([peak_mode], {"with_control", "without_control"}, f"{label}_peak_mode")
    validate_choices([peak_type], {"narrow", "broad"}, f"{label}_peak_type")
    if peak_mode not in call_peak_modes:
        raise ValueError(f"config['{label}_peak_mode'] is '{peak_mode}', but call_peak_modes does not include it.")
    if peak_type not in call_peak_types:
        raise ValueError(f"config['{label}_peak_type'] is '{peak_type}', but call_peak_types does not include it.")
    if peak_mode == "with_control" and not has_with_control:
        raise ValueError(f"config['{label}_peak_mode'] is 'with_control', but no effective controls are configured.")


def load_replicate_data(path):
    data = load_yaml(path)
    return {str(group): unique_list(samples) for group, samples in data.items()}


def valid_replicate_groups(replicate_data, treatments):
    treatment_set = set(treatments)
    return sorted(
        group
        for group, replicates in replicate_data.items()
        if len(replicates) >= 2 and all(sample in treatment_set for sample in replicates)
    )


def replicates_for_group(wildcards):
    replicates = REPLICATE_DATA[wildcards.group]
    if len(replicates) < 2:
        raise ValueError(f"Replicate group '{wildcards.group}' needs at least two samples.")
    return replicates


def replicate_peak_inputs(wildcards):
    return [peak_path(sample, REPLICATE_PEAK_MODE, REPLICATE_PEAK_TYPE) for sample in replicates_for_group(wildcards)]


def replicate_min_support(wildcards):
    return int(config.get("replicate_min_support", len(replicates_for_group(wildcards))))


def homer_outputs(targets, result_tag):
    report_dir = config.get("homer_report_dir", "homer_results")
    return expand(f"{report_dir}/{result_tag}/{{target}}/.homer_complete", target=targets)
