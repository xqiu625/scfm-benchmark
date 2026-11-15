"""Gate 2 runtime rate configuration — fail-closed.

This is the ONLY route by which an `expected_doublet_rate` may reach Scrublet.
It is deliberately a different artifact from the T7 evidence layer: T7 records
what could be reasoned about, this records what a human explicitly authorized.

**No approved runtime config exists.** Every T7 row carries
`runtime_eligible=false`, no cohort has an identifiable post-Gate-1 residual
prior, and rate selection is unresolved. A production run therefore terminates
here with an error, by design. Tests build temporary fixture configs.

Why fail-closed rather than "use a sensible default": the previous
implementation had two fallbacks that both failed silently and both produced
invalid science --

    if a.prior_counts and os.path.exists(a.prior_counts):   # missing file -> no prior
    n_recovered = n_cells_recovered if ... is not None else n_cells   # -> post-QC count

Six of eight cohorts were scored against post-Gate-1 counts because of them,
and nothing raised. Every condition below therefore raises rather than
degrades.
"""
import csv
import hashlib
import math
import os

SCHEMA_VERSION = 'gate2-runtime-config/1'

#: Platforms whose rate model has been characterized. An unknown platform is an
#: error, not an occasion to reuse another platform's number.
KNOWN_PLATFORMS = frozenset({'10x', 'seqwell'})

#: The processing stage the rate must describe. A rate estimated for one stage
#: and applied to another is the central defect this contract exists to block,
#: so the two are named per row and must agree.
GATE2_INPUT_STAGE = 'post_gate1_scrublet_input'

#: How the value was arrived at. Recorded so a modelled or extrapolated number
#: can never be read back as a measurement.
VALUE_BASES = frozenset({'measured', 'modelled', 'extrapolated', 'arbitrary'})

REQUIRED_COLUMNS = (
    'schema_version',
    'cohort',
    'channel_id',
    'platform',
    'gate2_input_stage',
    'expected_doublet_rate',
    'rate_estimand',
    'rate_source',
    'decision_reference',
    'runtime_eligible',
    'value_basis',
    'n_candidate_barcodes',
    'n_called_pre_qc',
    'n_post_demux_singlets',
    'n_gate1_retained',
    'n_donors',
    'donor_proportion_stage',
    'donor_proportions_checksum',
    'source_artifact',
    'source_checksum',
    'code_provenance',
)

#: Paths under these are evidence, never runtime configuration. Checked by path
#: component so a copy or symlink elsewhere named after T7 is still refused.
FORBIDDEN_PATH_MARKERS = (
    'T7_doublet_prior_evidence',
    'literature_audit',
    'rate-sensitivity-scenarios',
    'pooled-same-donor-pairing',
    'doublet_prior_capture_evidence',
    'prior_evidence_verdicts',
    'capture_evidence_join_key',
)


class Gate2ConfigError(RuntimeError):
    """Any defect in the runtime configuration. Always fatal, never downgraded."""


def _reject_evidence_path(path):
    norm = os.path.normpath(os.path.abspath(path))
    parts = norm.replace('\\', '/').split('/')
    for marker in FORBIDDEN_PATH_MARKERS:
        if any(marker in part for part in parts):
            raise Gate2ConfigError(
                f'{path} is part of the T7 evidence layer (matched {marker!r}) and must not be '
                'used as a Gate 2 runtime configuration. Every T7 row is runtime_eligible=false; '
                'promoting a scenario to a runtime parameter requires a separate, explicit, '
                'auditable decision recorded in docs/decisions.md.')


def _require_positive_int(row, col, channel, *, optional=False):
    raw = (row.get(col) or '').strip()
    if raw == '':
        if optional:
            return None
        raise Gate2ConfigError(f'channel {channel!r}: {col} is required and empty')
    try:
        val = int(raw)
    except ValueError:
        raise Gate2ConfigError(f'channel {channel!r}: {col}={raw!r} is not an integer')
    if val <= 0:
        raise Gate2ConfigError(f'channel {channel!r}: {col}={val} must be positive')
    return val


def _validate_row(row, cohort, path):
    channel = (row.get('channel_id') or '').strip()
    if not channel:
        raise Gate2ConfigError(f'{path}: a row has an empty channel_id')

    if (row.get('schema_version') or '').strip() != SCHEMA_VERSION:
        raise Gate2ConfigError(
            f'channel {channel!r}: schema_version={row.get("schema_version")!r}, '
            f'expected {SCHEMA_VERSION!r}')

    row_cohort = (row.get('cohort') or '').strip()
    if row_cohort != cohort:
        raise Gate2ConfigError(
            f'channel {channel!r}: cohort={row_cohort!r} but this run is cohort {cohort!r}')

    eligible = (row.get('runtime_eligible') or '').strip().lower()
    if eligible != 'true':
        raise Gate2ConfigError(
            f'channel {channel!r}: runtime_eligible={row.get("runtime_eligible")!r}. Only rows '
            'explicitly marked true may be used at runtime; there is no default and no '
            'permissive reading of a blank, "1", or "yes".')

    platform = (row.get('platform') or '').strip()
    if platform not in KNOWN_PLATFORMS:
        raise Gate2ConfigError(
            f'channel {channel!r}: platform={platform!r} is not one of {sorted(KNOWN_PLATFORMS)}. '
            'An uncharacterized platform must not borrow another platform\'s rate model.')

    stage = (row.get('gate2_input_stage') or '').strip()
    if stage != GATE2_INPUT_STAGE:
        raise Gate2ConfigError(
            f'channel {channel!r}: gate2_input_stage={stage!r}, expected {GATE2_INPUT_STAGE!r}')

    estimand = (row.get('rate_estimand') or '').strip()
    if not estimand:
        raise Gate2ConfigError(f'channel {channel!r}: rate_estimand is required')
    if GATE2_INPUT_STAGE not in estimand:
        raise Gate2ConfigError(
            f'channel {channel!r}: rate_estimand={estimand!r} does not describe '
            f'{GATE2_INPUT_STAGE!r}. A rate estimated for one processing stage must not be '
            'applied to another -- state the estimand for the matrix Scrublet is actually given.')

    raw_rate = (row.get('expected_doublet_rate') or '').strip()
    if raw_rate == '':
        raise Gate2ConfigError(f'channel {channel!r}: expected_doublet_rate is empty')
    try:
        rate = float(raw_rate)
    except ValueError:
        raise Gate2ConfigError(f'channel {channel!r}: expected_doublet_rate={raw_rate!r} is not a number')
    if not math.isfinite(rate):
        raise Gate2ConfigError(f'channel {channel!r}: expected_doublet_rate={rate} is not finite')
    if not (0.0 < rate < 1.0):
        raise Gate2ConfigError(f'channel {channel!r}: expected_doublet_rate={rate} outside (0, 1)')

    for col in ('rate_source', 'decision_reference', 'source_artifact', 'source_checksum',
                'code_provenance', 'donor_proportion_stage'):
        if not (row.get(col) or '').strip():
            raise Gate2ConfigError(f'channel {channel!r}: {col} is required and empty')

    basis = (row.get('value_basis') or '').strip()
    if basis not in VALUE_BASES:
        raise Gate2ConfigError(
            f'channel {channel!r}: value_basis={basis!r} not in {sorted(VALUE_BASES)}')

    n_donors = _require_positive_int(row, 'n_donors', channel)
    n_gate1 = _require_positive_int(row, 'n_gate1_retained', channel)
    optional_counts = {c: _require_positive_int(row, c, channel, optional=True)
                       for c in ('n_candidate_barcodes', 'n_called_pre_qc',
                                 'n_post_demux_singlets')}

    return {
        'channel_id': channel,
        'cohort': cohort,
        'platform': platform,
        'gate2_input_stage': stage,
        'expected_doublet_rate': rate,
        'rate_estimand': estimand,
        'rate_source': row['rate_source'].strip(),
        'decision_reference': row['decision_reference'].strip(),
        'value_basis': basis,
        'n_donors': n_donors,
        'n_gate1_retained': n_gate1,
        'donor_proportion_stage': row['donor_proportion_stage'].strip(),
        'donor_proportions_checksum': (row.get('donor_proportions_checksum') or '').strip(),
        'source_artifact': row['source_artifact'].strip(),
        'source_checksum': row['source_checksum'].strip(),
        'code_provenance': row['code_provenance'].strip(),
        **optional_counts,
    }


def row_checksum(entry):
    """Stable checksum of the authorized row, bound into the checkpoint so a
    changed rate, source or decision reference invalidates restored results."""
    payload = '|'.join(f'{k}={entry[k]!r}' for k in sorted(entry))
    return hashlib.sha256(payload.encode()).hexdigest()


def load_runtime_config(path, cohort, required_channels):
    """Load and fully validate the runtime config. Raises Gate2ConfigError on any
    defect; never returns a partial or defaulted configuration.

    `required_channels` is the exact capture set about to be scored. Every one
    must be present, and any extra row is also an error -- a config that does
    not describe this run is not a config for this run.
    """
    if path is None or str(path).strip() == '':
        raise Gate2ConfigError(
            'no --runtime-config given. Gate 2 has no default expected_doublet_rate and no '
            'fallback to cell counts. A rate must come from an explicitly authorized runtime '
            'configuration; see docs/data-contracts.md.')
    _reject_evidence_path(path)
    if not os.path.exists(path):
        raise Gate2ConfigError(f'runtime config not found: {path}')

    with open(path, newline='', encoding='utf-8-sig') as fh:
        reader = csv.DictReader(fh)
        missing_cols = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
        if missing_cols:
            raise Gate2ConfigError(f'{path}: missing required columns {missing_cols}')
        rows = list(reader)
    if not rows:
        raise Gate2ConfigError(f'{path}: contains no rows')

    by_channel = {}
    for row in rows:
        entry = _validate_row(row, cohort, path)
        ch = entry['channel_id']
        if ch in by_channel:
            raise Gate2ConfigError(f'{path}: duplicate channel_id {ch!r}')
        by_channel[ch] = entry

    required = list(dict.fromkeys(str(c) for c in required_channels))
    missing = [c for c in required if c not in by_channel]
    if missing:
        raise Gate2ConfigError(
            f'{path}: no runtime config row for {len(missing)} capture(s), '
            f'first few: {missing[:5]}')
    extra = sorted(set(by_channel) - set(required))
    if extra:
        raise Gate2ConfigError(
            f'{path}: {len(extra)} row(s) describe captures not in this run, '
            f'first few: {extra[:5]}')

    for entry in by_channel.values():
        entry['config_row_checksum'] = row_checksum(
            {k: v for k, v in entry.items() if k != 'config_row_checksum'})
    return by_channel
