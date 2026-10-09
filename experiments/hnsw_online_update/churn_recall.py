"""hnswlib online-update experiment.

Measures recall@10 of hnswlib under online updates (mark-delete, replace-deleted,
in-place vector update) and compares against a freshly built index on the same
live set. Ground truth is exact brute force over the live vectors.
"""
import time

import hnswlib
import numpy as np

DIM = 64
N = 20_000
NQ = 500
K = 10
M = 16
EF_C = 200
EFS = (20, 50, 100)
THREADS = 4
SEED = 0

rng = np.random.default_rng(SEED)


# Gaussian mixture shared by base, queries and inserts (closer to real
# embeddings than isotropic noise, and keeps queries in-distribution).
CENTERS = rng.normal(size=(100, DIM)).astype(np.float32)


def make_data(n, centers=CENTERS):
    assign = rng.integers(0, len(centers), size=n)
    return (centers[assign] + rng.normal(size=(n, DIM))).astype(np.float32)


def exact_knn(live_ids, live_vecs, queries):
    d = (queries ** 2).sum(1)[:, None] - 2 * queries @ live_vecs.T + (live_vecs ** 2).sum(1)[None, :]
    idx = np.argpartition(d, K, axis=1)[:, :K]
    return [set(live_ids[row]) for row in idx]


def recall(index, gt, queries):
    out = {}
    for ef in EFS:
        index.set_ef(ef)
        labels, _ = index.knn_query(queries, k=K, num_threads=THREADS)
        out[ef] = np.mean([len(set(l) & g) / K for l, g in zip(labels, gt)])
    return out


def fresh_index(ids, vecs):
    idx = hnswlib.Index(space="l2", dim=DIM)
    idx.init_index(max_elements=len(ids), M=M, ef_construction=EF_C, random_seed=SEED)
    idx.add_items(vecs, ids, num_threads=THREADS)
    return idx


def fmt(r):
    return "  ".join(f"ef={ef}:{v:.4f}" for ef, v in r.items())


def main():
    queries = make_data(NQ)
    store = {}  # label -> vector (live set)
    next_label = 0

    base = make_data(N)
    labels = np.arange(N)
    store.update(zip(labels.tolist(), base))
    next_label = N

    idx = hnswlib.Index(space="l2", dim=DIM)
    idx.init_index(max_elements=N, M=M, ef_construction=EF_C, random_seed=SEED,
                   allow_replace_deleted=True)
    t = time.perf_counter()
    idx.add_items(base, labels, num_threads=THREADS)
    print(f"build {N} pts: {time.perf_counter() - t:.2f}s ({THREADS} threads)")

    def live():
        ids = np.fromiter(store.keys(), dtype=np.int64)
        return ids, np.stack([store[i] for i in ids])

    ids, vecs = live()
    print(f"[round 0] online      {fmt(recall(idx, exact_knn(ids, vecs, queries), queries))}")

    # --- 1) mark-delete only: graph keeps tombstones ---
    frac = 0.3
    victims = rng.choice(ids, size=int(frac * len(ids)), replace=False)
    for v in victims:
        idx.mark_deleted(int(v))
        del store[int(v)]
    ids, vecs = live()
    gt = exact_knn(ids, vecs, queries)
    print(f"[delete 30%] tombstone {fmt(recall(idx, gt, queries))}")
    print(f"[delete 30%] fresh     {fmt(recall(fresh_index(ids, vecs), gt, queries))}")

    # --- 2) churn rounds: delete 20% then refill slots via replace_deleted ---
    for rnd in range(1, 6):
        if rnd > 1:
            ids, _ = live()
            victims = rng.choice(ids, size=int(0.2 * N), replace=False)
            for v in victims:
                idx.mark_deleted(int(v))
                del store[int(v)]
        n_new = N - len(store)
        new_vecs = make_data(n_new)
        new_labels = np.arange(next_label, next_label + n_new)
        next_label += n_new
        t = time.perf_counter()
        idx.add_items(new_vecs, new_labels, num_threads=THREADS, replace_deleted=True)
        dt = time.perf_counter() - t
        store.update(zip(new_labels.tolist(), new_vecs))
        ids, vecs = live()
        gt = exact_knn(ids, vecs, queries)
        print(f"[churn {rnd}] replace {n_new} pts in {dt:.2f}s  online {fmt(recall(idx, gt, queries))}")
    print(f"[churn 5] fresh                           {fmt(recall(fresh_index(ids, vecs), gt, queries))}")

    # --- 3) in-place update of existing labels (vector drift) ---
    upd = rng.choice(ids, size=int(0.2 * N), replace=False)
    new_vecs = make_data(len(upd))
    t = time.perf_counter()
    idx.add_items(new_vecs, upd, num_threads=THREADS)
    dt = time.perf_counter() - t
    for l, v in zip(upd.tolist(), new_vecs):
        store[l] = v
    ids, vecs = live()
    gt = exact_knn(ids, vecs, queries)
    print(f"[update 20%] in-place {dt:.2f}s  online {fmt(recall(idx, gt, queries))}")
    print(f"[update 20%] fresh                {fmt(recall(fresh_index(ids, vecs), gt, queries))}")


if __name__ == "__main__":
    main()
