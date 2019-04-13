import time
from multiprocessing import Manager, Process, cpu_count, Pool

from tqdm import tqdm

from iter_for import batch_iter2


def assign_work(func, args, n_iter, num_thread, verbose=1):
    """
    parallel func which doesn't use different params. Each run of func is independent.
    With progress bar(tqdm)
    :param func: last two params must be sub_i, q. first is the iteration times, second is the shared list to stodge result.
    :param args: params passed to func.
    :param n_iter: total iteration times.
    :param num_thread: number of threads.
    :param verbose: verbose of this function
    :return:
    """
    pbar = None
    if verbose:
        pbar = tqdm(total=n_iter)
    # initiate a progress bar according verbose
    if num_thread == 0:
        # if number of thread equal to 0, then use all threads of computer.
        num_thread = cpu_count()
    chunks = batch_iter2(range(n_iter), num_thread)
    # split total iteration times into a list which contains `num_thread` list contains almost same number.(could be simple)
    # but for more robust `batch_iter2`, it make a little complicated.
    manager = Manager()
    q = manager.list()
    # create a shared list between processes used to stodge results
    for _ in chunks:
        p = Process(target=func,
                    args=(*args,
                          len(_),
                          q))
        p.daemon = True
        p.start()

    while 1:
        # check all processes has been complete. like .join()
        # and also used to update pbar
        if pbar:
            pbar.update(len(q) - pbar.n)
        if len(q) == n_iter:
            if pbar:
                pbar.update(n_iter - pbar.n)
            break
            # jump out the while loop until length of q equal to the total iteration times.
    return q


def assign_work_pool(func, differ_args, num_thread):
    if num_thread == 0:
        # if number of thread equal to 0, then use all threads of computer.
        num_thread = cpu_count()
    with Pool(processes=num_thread) as pool:
        results = list(tqdm(pool.imap(func, differ_args), total=len(differ_args)))

    return results


if __name__ == '__main__':

    def cal_for(x, y, sub_i, q):
        for _ in range(sub_i):
            t = x * y
            time.sleep(0.1)
            q.append(t)


    # t1 = time.time()
    # results = assign_work(cal_for,
    #                       (10, 20),
    #                       n_iter=1000,
    #                       num_thread=7)
    # print(len(results))
    # print(time.time() - t1)
    ############################################################
    def cal_for2(x, y):
        t = x * y
        time.sleep(0.1)
        return x, t


    def accesory(arg):
        return cal_for2(*arg)


    t1 = time.time()
    results = assign_work_pool(cal_for2,
                               (*[(x, y) for x, y in zip(range(100),
                                                         range(78, 178))],),
                               num_thread=7)
    print(len(results))
    print(time.time() - t1)
    for _1, _2 in results:
        print(_1, _2)
