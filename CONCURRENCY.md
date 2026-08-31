# Concurrency contract

DR3 supports concurrent allocation, use, and destruction of independently
owned vector values under the rules below. A clean stress or sanitizer run is
supporting evidence, not a guarantee that arbitrary sharing is safe.

| API or state | Sharing class | Mutation and ownership rule | Teardown rule |
|---|---|---|---|
| `Vec` and `VecD` values | Independent instances | Different threads may own different values. Concurrent mutation of one value is unsupported. Pool allocation/free is synchronized. | A value may be destroyed on another joined worker because pool return is synchronized. |
| Allocator registry and pools | Internally synchronized | Allocation, return, policy creation, and cleanup hold the registry mutex. Callers must not retain raw pool pointers after return. | Join workers and destroy all pool-backed values before allocator-global cleanup. Cleanup rejects live blocks with `std::logic_error`. |

The allocator mutex is acquired only when storage is obtained or returned. It
is not acquired by SIMD arithmetic on an already-owned value.

The registry map and its mutex are intentionally heap-allocated
process-lifetime objects. This avoids cross-translation-unit static destruction
ordering hazards. Explicit cleanup still releases all pools once their live
block counts reach zero. The `AllAllocatorsGuard` teardown path is nonthrowing;
if a static value is still live during process exit, its pool is left for the
operating system to reclaim instead of terminating the process.

## Unsupported use

- Concurrent mutation of the same `Vec` or `VecD` is unsupported.
- Retaining or using a raw pool pointer after returning it is unsupported.
- Calling `freeAllAllocators(...)` while workers or pool-backed values are
  alive is unsupported and throws when a live block can be detected.
- Removing a size-specific policy while one of its blocks is live is
  unsupported and throws. Libraries built in Debug mode also diagnose unknown
  and duplicate pool returns; Release libraries preserve the legacy no-op
  behavior even when consumed by a Debug application.

Tests use independent values, synchronized starts, and joined workers. The
focused ThreadSanitizer gate exercises allocator and independently owned
vector operations only.
