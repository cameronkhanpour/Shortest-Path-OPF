# Codex Session Log

This file records the changes I personally made during this chat session, along with the reason for each change.

## Manual Code Changes

### `hit_and_run.jl`

I made three rounds of changes in this file.

1. I replaced the earlier line-search logic with a hybrid boundary search based on:
   - `_bisect_boundary`
   - `_find_boundary_by_hybrid`
   - `_segment_is_connected`

   Reason:
   The original sampling logic was not reliably respecting disconnected or holed feasible sets. On `case9dongchan`, that showed up as samples leaking into the artificial ellipsoidal hole. The new logic explicitly tries to keep the chain on the connected feasible segment containing the current point, up to the chosen scan resolution `max_step`.

2. I extended the direction generator and sampler to support an optional pilot direction basis:
   - `_random_direction(dim, direction_basis)`
   - `next_sample(...; direction_basis=...)`

   Reason:
   For standard MATPOWER cases like `pglib_opf_case14_ieee.m` and `pglib_opf_case30_ieee.m`, isotropic ambient directions produced essentially no movement. Empirically, the feasible set was varying along correlated lower-dimensional directions. The optional basis lets hit-and-run move along those pilot-estimated directions instead of stalling numerically.

3. Later in the session, I replaced the active sampling kernel again with a full-line interval-union sampler:
   - `_line_box_parameter_range`
   - `_locate_feasibility_transition`
   - `_feasible_intervals_along_line`
   - `_sample_from_interval_union`

   I also changed `next_sample(...)` so that it samples uniformly from the full resolved feasible union along a random line, rather than only from the connected segment containing the current point.

   Reason:
   The previous connected-segment kernel was not a correct full-set hit-and-run method for nonconvex AC OPF feasible sets. It prevented jumps across separated feasible intervals on the same line and therefore did not target the full feasible space. The new kernel is much closer to the mathematically correct hit-and-run transition for a general nonconvex set, subject to the line-resolution approximation.

4. I added handling for zero-width box coordinates inside `next_sample(...)`.

   Reason:
   On `pglib_opf_case60_c.m`, at least one control coordinate had zero box width. Treating that coordinate as movable caused almost every random direction to collapse the line-box intersection to a point, which effectively froze the sampler.

## `convexity_metric.jl`

I made several substantial changes here.

1. I refactored the file around a reusable feasibility oracle:
   - `_create_feasibility_oracle(...)`

   Reason:
   The original metric workflow was too case-specific and too expensive to reuse cleanly across cases. The refactor lets the same loaded case data and oracle be reused for sampling, metric evaluation, and slicing.

2. I changed the feasibility oracle to warm-start the power flow solve using the last successful PF state and the base operating point.

   Reason:
   Nearby feasible points were sometimes being mislabeled as infeasible because Newton power flow was restarted poorly. The warm start improves the numerical honesty of the feasibility test without changing the mathematical constraint set.

3. I rewrote the main metric path to compute ratios directly from the shortest-path functions already in the repository:
   - `compute_ratio(...)`
   - `compute_convexity_metric(...)`

   Reason:
   This made the metric logic explicit and inspectable inside `convexity_metric.jl` instead of depending on a less transparent wrapper flow.

4. I made the metric reporting more honest by explicitly skipping failed pair evaluations instead of silently averaging only over successful ones.

   Reason:
   A reported `C(S)` value is misleading if many pair computations failed and that fact is hidden.

5. I added a more general analysis workflow:
   - `analyze_case(...)`
   - `analyze_cases(...)`
   - updated `main(...)`

   Reason:
   This was needed to run multiple MATPOWER cases consistently and to separate full-space analysis from visualization.

6. I added 2D slice utilities and slice plotting support:
   - `_sample_mean`
   - `_choose_slice_center`
   - `_select_slice_positions`
   - `_slice_reference_u0`
   - `_slice_bounds_from_samples`
   - `_control_label`
   - coordinate-slice plotting helpers

   Reason:
   For higher-dimensional MATPOWER cases, plotting the full feasible set directly is not possible. The slice workflow gives an honest visualization of a fixed 2D affine section rather than a misleading raw projection.

7. I added pilot-OPF helpers:
   - `_opf_point_from_solution(...)`
   - `_randomized_opf_point_cloud(...)`
   - `_pilot_direction_basis(...)`

   Reason:
   On `case14`, `case30`, and later `case60_c`, plain ambient hit-and-run did not move meaningfully. I used a cloud of randomized OPF solutions to estimate the directions along which the feasible set actually varies, then used that information to seed and guide sampling.

8. I updated sampling calls so the full-space sampler can use:
   - a chosen initial point
   - custom bounds
   - an optional pilot direction basis

   Reason:
   This was necessary to support both the full-space runs and the reduced 2D slice runs with the same sampling code.

9. Later in the session, I changed the active full-space workflow again:
   - added `_project_to_box`
   - added `_active_box_bound_count`
   - added `_find_sampling_seed`
   - changed `sample_nonconvex_space(...)` to clamp seeds/samples to the box
   - changed `analyze_case(...)` and `main(...)` to focus on full-space runs only, with no slice plotting in the active path

   Reason:
   The `case60_c` debugging showed that the main problem was still the sampler, not the shortest-path routine. The OPF-based starting points were sitting on many active box bounds, so the walk barely moved. The new seed search uses randomized OPF solutions and feasible averages of them to find a more interior feasible point before hit-and-run starts. I also removed the pilot-basis and plotting workflow from the active analysis path so the code now focuses on the full-space sampler directly.

## `convexity_metric_plots.jl`

I created this file during the session.

1. I moved the 2D feasible-space and path plotting logic out of `convexity_metric.jl`.
2. I changed the plot styling to be closer to the look of `barrier_case_pdf_v2.jl` and `barrier_case_time_v2.jl`.
3. I made the plotting save both `.png` and `.pdf` outputs.

Reason:
The plotting code had become too bulky inside `convexity_metric.jl`, and the original visuals were weaker than the existing barrier-case plots in the repository. The new file keeps the plotting separate and reusable.

## Generated Outputs

I also generated new output files while running experiments in this session. These were not hand-edited, but they were created by the code I ran.

### `case9dongchan`

I generated updated feasible-space and path figures after fixing the sampler.

### `pglib_opf_case14_ieee`

I generated:
- coordinate-slice feasible-space plots
- coordinate-slice shortest-path plots

### `pglib_opf_case30_ieee`

I generated:
- coordinate-slice feasible-space plots
- coordinate-slice shortest-path plots

### `pglib_opf_case60_c`

I generated:
- a coordinate-slice feasible-space plot during the earlier slice-based workflow
- later, a full-space `case60_c` metric run using the redesigned sampler

## Files I Intentionally Did Not Edit

I did not modify `shortest_path_v2.jl`, per your instruction.

I also did not intentionally edit any other pre-existing Julia source file besides:
- `convexity_metric.jl`
- `hit_and_run.jl`

The new Julia source file I added was:
- `convexity_metric_plots.jl`
