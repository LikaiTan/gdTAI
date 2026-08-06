# gdTAI v4.0 Promotion Decision

Status: experimental candidate, not promoted as the default model.

Rationale:

- v4 high-F1 improves v3 Round12 full-primary recall (`0.792` vs `0.756`) and sorted cord-blood recall (`0.750` vs `0.675`).
- v4 high-F1 improves TRDC+TRDV- validation recall (`0.295` vs v3 `0.107`) without using a hard TRDV requirement.
- v4 high-F1 has lower full-atlas NK predicted burden than v3 (`835` vs `1,723`).
- v4 high-F1 has higher paired/any TCRAB warning burden than v3 (`3,619` vs `2,748`) but lower than v2 high-purity (`5,184`).
- v4 still does not recover sorted cord-blood gdT cells as well as v2 high-purity (`0.750` vs `0.861`), so it should be treated as an iteration result, not a release replacement.
