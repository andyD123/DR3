# Concurrency contract

DR3 supports concurrent calls only under the ownership rules below. A clean
stress or sanitizer run is supporting evidence, not a guarantee that arbitrary
sharing is safe.

| API or state | Sharing class | Mutation and ownership rule | Teardown rule |
|---|---|---|---|
| `Vec` and `VecD` values | Independent instances | Different threads may own different values. Concurrent mutation of one value is unsupported. Pool allocation/free is synchronized. | A value may be destroyed on another joined worker because pool return is synchronized. |
| Allocator registry and pools | Internally synchronized | Allocation, return, policy creation, and cleanup hold the registry mutex. Callers must not retain raw pool pointers after return. | `freeAll()` is permitted only after all pool-backed values are destroyed and workers are joined. |
| Immutable curves | Concurrent read-only | Concurrent evaluation is allowed only while no thread resets or mutates curve data. | The curve owner outlives all readers. |
| Tree pricers | Independent calls | Inputs and outputs are caller-owned; temporary `Vec` allocations use the synchronized pool. | All calls finish before allocator-global cleanup. |
| European PDE solver | Independent calls | Each call owns its grid, result, and `PdeWorkspace`. Do not share mutable result vectors. | Results may outlive the call and use standard container ownership. |
| American PDE solver | Not implemented | Concurrency classification is required with its implementation. | To be specified. |
| Batched PDE solver | Not implemented | Concurrency classification is required with its implementation. | To be specified. |
| Parallel executor | Not implemented | Workers will own disjoint output ranges and one workspace each. | The executor joins every worker before returning or rethrowing. |

The allocator mutex is acquired only when storage is obtained or returned. It
is not acquired by SIMD arithmetic, tree rollback, or PDE time stepping.

## Unsupported use

- Concurrent reads and writes to the same `Vec`, `VecD`, result vector, or
  workspace are unsupported.
- Curve mutation while another thread reads the curve is unsupported.
- Calling allocator-global cleanup while workers or pool-backed values are
  alive is unsupported.
- Sharing a future solver object is unsupported until its contract explicitly
  says otherwise.

Debug assertions may be added for ownership violations that can be detected
without adding state to numerical hot loops. Tests use independent values,
fixed work partitions, synchronized starts, and joined workers.
