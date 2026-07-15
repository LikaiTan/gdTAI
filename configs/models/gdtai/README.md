# gdTAI model registry

`model_registry.csv` is the release and experiment index. A model is the
default only when its status is `promoted_default`; directory names alone do
not imply promotion.

Before loading a pickle, verify its SHA256 against this registry. Historical
pickles reference modules by their pre-package names, so use
`tnk_atlas.model_io.load_trusted_pickle` or a documented inference entry point.

The canonical Git gdTAI v3 release is Round 14 at threshold 0.936. The current
dirty workspace contains a pre-existing Round 12 override at the same path.
Both hashes are recorded so validation can warn without discarding work. Do not
promote or commit the override implicitly; reconcile model bytes, manifest,
metrics, and promotion decision together before the next release.
