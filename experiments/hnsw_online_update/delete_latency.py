"""Query latency and recall of hnswlib as the tombstone (mark_deleted) ratio grows.

Deleted nodes are still expanded during search, so a growing tombstone ratio
should cost latency even when recall holds up.
"""
import time

import hnswlib
import numpy as np

DIM, N, NQ, K, M, EF_C, EF = 64, 20_000, 1_000, 10, 16, 200, 20
rng = np.random.default_rng(1)
centers = rng.normal(size=(100, DIM)).astype(np.float32)


def make(n):
    return (centers[rng.integers(0, 100, n)] + rng.normal(size=(n, DIM))).astype(np.float32)


data, queries = make(N), make(NQ)
idx = hnswlib.Index(space="l2", dim=DIM)
idx.init_index(max_elements=N, M=M, ef_construction=EF_C, random_seed=1)
idx.add_items(data, np.arange(N), num_threads=4)
idx.set_ef(EF)

order = rng.permutation(N)
deleted = 0
print("deleted_ratio  recall@10  us/query(1 thread)")
for ratio in (0.0, 0.1, 0.3, 0.5, 0.7, 0.9):
    target = int(ratio * N)
    for v in order[deleted:target]:
        idx.mark_deleted(int(v))
    deleted = target
    live = np.sort(order[deleted:])
    lv = data[live]
    d = (queries ** 2).sum(1)[:, None] - 2 * queries @ lv.T + (lv ** 2).sum(1)[None, :]
    gt = live[np.argpartition(d, K, axis=1)[:, :K]]
    idx.knn_query(queries, k=K, num_threads=1)  # warm-up
    times = []
    for _ in range(5):
        t = time.perf_counter()
        labels, _ = idx.knn_query(queries, k=K, num_threads=1)
        times.append(time.perf_counter() - t)
    us = np.median(times) / NQ * 1e6
    rec = np.mean([len(set(a) & set(b)) / K for a, b in zip(labels, gt)])
    print(f"{ratio:>13.1f}  {rec:9.4f}  {us:8.1f}")
