# Fix Sampling Plan

## Current Status

The current full-space sampler does **not** look reliable enough to run a 50-sample convexity-metric experiment on `MATPOWER/pglib_opf_case39_epri.m`.

Observed probe behavior:

- control-space dimension: `19`
- interiorised seed found, but it still has `8` active box bounds
- pilot OPF cloud spread: about `10.26`
- hit-and-run probe spread after sampling: about `5.85e-7`
- all sampled points remain feasible, but the chain is effectively not moving
- `Power flow did not converge! Using flat start as the operating point` appears while loading the case

Interpretation:

- The sampler is no longer obviously sampling outside feasibility.
- The main problem for `case39_epri` is still **severe lack of movement / mixing**.
- A `C(S)` value computed from this chain would not be trustworthy.

## Likely Root Causes

1. The chain is starting from a point that is still too close to the boundary.

   Even after the current seed search, the best feasible seed I found still has many active box bounds. In high dimension, that makes most random ambient directions nearly useless.

2. The current bounding box is too crude to support efficient hit-and-run on difficult cases.

   The sampler uses box bounds taken from the selected control variables only. If the true feasible set inside that box is thin, curved, or strongly correlated, random ambient lines will usually have extremely short feasible segments.

3. The feasibility oracle is still numerically fragile on harder cases.

   `case39_epri` is already showing a flat-start PF fallback in the case parsing flow. That is a warning sign that the PF side may be making line scans noisier than they should be.

4. The current line-interval reconstruction is still resolution-limited.

   The full-line interval-union kernel is structurally better than the earlier connected-segment kernel, but it still depends on `max_step` and repeated PF solves. On a numerically delicate case, that can still under-resolve line geometry and suppress movement.

## Plan

### 1. Build a deliberately interior feasible seed

Goal:
Start hit-and-run from a point that is not sitting on many active box bounds.

Work:

- Add a seed-improvement routine that takes feasible OPF points and solves a local feasibility-centered problem.
- Objective should reward distance from:
  - control-variable box bounds
  - active OPF inequality constraints
- Use the existing feasibility oracle and PF machinery to verify the candidate.

Success condition:

- active control-box bound count is near zero on difficult cases
- small random perturbations around the seed remain feasible at a much higher rate

### 2. Reduce the sampling space to truly free coordinates

Goal:
Do not waste random directions on coordinates that are effectively fixed.

Work:

- Detect coordinates with zero or near-zero allowable range from the box data.
- Remove them from the sampler state rather than only zeroing their direction component.
- Map sampled reduced-space points back to the full control vector before feasibility checks.

Success condition:

- every random direction used by hit-and-run acts only on coordinates that can actually move

### 3. Improve the full-line interval oracle

Goal:
Make the line-union computation more trustworthy on hard cases.

Work:

- Replace the fixed `max_step` scan with an adaptive refinement strategy.
- Start coarse, then recursively refine around sign changes or numerically ambiguous feasibility transitions.
- Cache PF warm starts along the line so nearby line evaluations reuse nearby PF states.

Success condition:

- repeated scans of the same line with tighter tolerances produce stable interval unions
- interval lengths do not collapse just because of scan resolution

### 4. Strengthen the feasibility oracle for line scans

Goal:
Reduce false negatives during hit-and-run movement.

Work:

- Track and reuse the last successful PF state along each scanned line, not just globally.
- Add a second-stage PF retry strategy for difficult points:
  - previous line point
  - previous accepted sample
  - base operating point
  - flat start as last resort
- Log PF failures separately from true OPF-constraint infeasibility.

Success condition:

- lower PF-failure rate along scanned lines
- clearer distinction between numerical failure and actual infeasibility

### 5. Add sampler diagnostics before trusting a metric run

Goal:
Prevent misleading metric runs on badly mixed chains.

Work:

- Report:
  - max pairwise sample distance
  - mean accepted step length
  - fraction of random directions that yield nontrivial feasible line intervals
  - number of active bounds at the seed
- Add a simple “movement sanity check” gate before computing `C(S)`.

Success condition:

- the code refuses or warns on chains whose movement scale is far too small relative to the OPF cloud or box scale

### 6. Re-test on `pglib_opf_case39_epri.m`

Goal:
Confirm that the repaired sampler actually moves before trying larger experiments.

Work:

- run a cheap probe:
  - `5-10` samples
  - modest burn-in
- verify:
  - all sampled points are feasible
  - sample spread is meaningful
  - accepted steps are not at `1e-7` scale
- only then run the full `50`-sample convexity-metric experiment

Success condition:

- the sample spread is orders of magnitude larger than the current `5.85e-7`
- the chain explores a visibly nontrivial region of the feasible control space

## Immediate Recommendation

Do **not** use the current sampler result for `pglib_opf_case39_epri.m`.

The next engineering target should be:

1. interior seed search
2. adaptive line-interval oracle
3. stronger PF retry logic for line scans

Only after those three pieces are in place should we rerun `case39_epri` with `50` hit-and-run samples.
