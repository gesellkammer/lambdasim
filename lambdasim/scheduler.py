from concurrent import futures
import subprocess
import threading

_SCHEDULER = None

class _Scheduler(object):
    def __init__(self):
        self.threadpool = futures.ThreadPoolExecutor(max_workers=5)

    def check_output(self, *args, **kws):
        """
        like subprocess.check_output, returns a Future
        """
        return self.threadpool.submit(subprocess.check_output, *args, **kws)

    def popen(self, *args, **kws):
        return self.threadpool.submit(subprocess.Popen, *args, **kws)

    def subproc_call(self, *args, **kws):
        return self.threadpool.submit(subprocess.call, *args, **kws)

    def call_now(self, fn, args=[], kws={}):
        """
        Call fn now, returns a Future
        """
        return self.threadpool.submit(fn, *args, **kws)

def wrap_subproc(subproc):
    fut = futures.Future()
    def _func():
        subproc.wait()
        fut.set_result(subproc.poll())
    get_sched().call_now(_func)
    return fut

def call_now(fn, args=(), kws={}):
    return get_sched().call_now(fn, args, kws)

def popen(*args, **kws):
    return get_sched().popen(*args, **kws)

def subproc_call(*args, **kws):
    return get_sched().subproc_call(*args, **kws)

def check_output(*args, **kws):
    return get_sched().check_output(*args, **kws)

def call_later(seconds, fn, args=(), kws={}):
    """
    Returns --> a Future
    """
    fut = futures.Future()
    def _func(*args, **kws):
        out = fn(*args, **kws)
        fut.set_result(out)
    threading.Timer(seconds, _func, args=args, kwargs=kws).start()
    return fut


def get_sched():
    global _SCHEDULER
    if _SCHEDULER is None:
        _SCHEDULER = _Scheduler()
    return _SCHEDULER
