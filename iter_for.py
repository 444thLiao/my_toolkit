def batch_iter(iter, batch_size):
    # generating batch according batch_size
    n_iter = []
    batch_d = 0
    for batch_u in range(0, len(iter), batch_size):
        if batch_u != 0:
            n_iter.append(iter[batch_d:batch_u])
        batch_d = batch_u
    n_iter.append(iter[batch_d: len(iter) + 1])
    return n_iter


def batch_iter2(iter, num_batch):
    # generating batch according num_batch
    n_iter = []
    idx = []
    num_iter = len(iter)
    batch_size = num_iter // num_batch
    batch_d = 0
    for batch_u in range(0, num_iter + batch_size, batch_size):
        if batch_u != 0:
            n_iter.append(iter[batch_d:batch_u])
            idx.append((batch_d, batch_u))
        batch_d = batch_u
    if len(n_iter) > num_batch:
        final_v = [_ for v in n_iter[num_batch - 1:] for _ in v]
        n_iter = n_iter[:num_batch]
        n_iter[-1] = final_v
    return n_iter
