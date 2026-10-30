# Shared campaign objectives

In a saved Studio graph, enable **Use a shared campaign objective across design branches**, select one or more fresh compatible branches, and choose the count, timer, or passing-count objective. The branches must resolve to the same Target Structure node. For passing-count goals, configure independent prediction and explicit scientific criteria; graph validation rejects mismatched predictors, criteria, or target mappings.

The graph-level ledger tracks the selected objective. Passing-count goals count unique binder sequences meeting the configured criteria; count and timer goals do not require filtering. Each branch still runs through its normal campaign controller, native workflow, scheduler jobs, seeds, attempt ledger, receipts, and recovery. Native acceptance, shared downstream acceptance, prediction provenance, duplicate sequences, and rejection reasons remain separate in candidate Results. The shared panel reports the objective, GPU compute, attempts, and per-branch progress.

The shared controls bound total attempts, overall concurrent jobs, generation and validation concurrency, optional per-branch attempts/concurrency, and optional total GPU compute. Generator placement defaults to automatic compatible workstations and can be restricted to the current workstation. Generation remains one native attempt per execution job. Reaching the target or a shared budget prevents further dispatch; already active work is allowed to reconcile so its receipts and partial Results are retained.

Independent per-generator objectives remain available when the graph-level shared objective is disabled. Replenishing passing-count and timed objectives currently require the verified BindCraft2 or RFD4 campaign adapters. Generated-count objectives also support registered fixed-output native generator branches. Scientific acceptance criteria apply only to passing objectives. Native multi-attempt packing and shared structural diversity objectives are not implemented.

## Count, timer, and optional acceptance scopes (October 2026)

A campaign is an execution scope around the selected workflow/pipeline branches.
Its stopping objective is separate from its scientific nodes:

- `objective_type: generated`: stop dispatch at the generated-design target, subject to the attempt/time budgets. Filtering can be omitted and applied to retained designs later.
- `objective_type: time`: produce until `max_runtime_seconds` or `deadline_at` is reached. It has no count target and requires no validation stages or acceptance filter.
- `objective_type: accepted`: count unique designs meeting explicitly configured criteria. Optionally set `acceptance_node` to the ID of the pipeline's existing Shared scientific filters / Validation & Acceptance node. Its binder format, target type, mandatory filters and informational/ranking fields are snapshotted through the same acceptance backend. The existing independent prediction stage supplies the measurements. Missing evidence cannot pass. The resolved decision and prediction provenance are retained and reused at that node.

Count and passing objectives can also carry an optional duration/deadline. The first target/limit reached stops new dispatch. Active work finishes and remains in Results, so batch outputs can exceed a count target. A timed objective reports time completion separately from generated, accepted, rejected and unresolved counts. Attempt/GPU budgets remain safety bounds; exhausting one before the goal does not imply scientific success. No further campaign is launched automatically after completion.

The input tile/node campaign editor and the shared pipeline panel use the same duration/deadline controls. UI deadlines use Europe/Rome; workflow API `deadline_at` requires an ISO timestamp with an explicit timezone offset. Timers start at the durable campaign creation time and do not reset during recovery. For example, a native workflow can use:

```json
{"enabled": true, "objective_type": "time", "max_runtime_seconds": 1800,
 "deadline_at": "2026-10-31T08:00:00+01:00", "max_attempts": 1000}
```

Shared pipeline scopes use `graph.shared_campaign` with the same time fields and selected `generator_nodes`. RFD4's standalone design form also exposes these optional settings, including sequence/IDR condition-spec inputs. The workflow schema documents the nested campaign fields and the normal workflow validation API rejects invalid timers before submission.

Replenishing timed/accepted native campaigns remain limited to registered verified campaign adapters (currently BindCraft2 and RFD4). Other generator branches retain their supported fixed-output native routes. This change reuses the existing graph controller, execution ledger and scheduler; it does not establish new native multi-design-per-attempt support. Existing saved campaign-stage criteria remain supported for compatibility.

## Independent processing batches

Studio batches mean independent items per execution job: 1,000 inputs at 10 per batch produce 100 batches. Compatible generator, co-folding, sequence-design and per-candidate analysis nodes share the existing batch controller, concurrency, placement and recovery controls. Independent analysis settings are available in the pipeline processing panel. Optional dispatch duration/deadline limits stop new batches for each partitioned step; already submitted work is reconciled, and unprocessed input counts stay visible. Recovery retains the original step start time. Whole-scope design objectives remain separate from these per-step limits.

Docking retains its existing distributed controller and worker batches, with the same items-per-batch terminology. Existing clustering (including Tanimoto) and neural-network training are unchanged. Training uses one selected GPU when available, otherwise CPU; its training minibatch is a separate concept.

This is not universal standalone-tile campaign support: standalone workflows retain their native collection behavior. The workflow schema reports which pipeline nodes support partitioning and dispatch limits. Unsupported objectives are not implied by the presence of a native batch-size field. Filters remain optional separate operations.
