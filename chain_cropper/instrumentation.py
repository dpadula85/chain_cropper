"""instrumentation.py - logging, debugging, and timing helpers for chain_cropper.

Ported from polymer_couplings/instrumentation.py, which was written with
zero dependency on that package's own domain logic specifically so it
could be reused here.
"""

import sys
import pdb
import time
import logging
import argparse as arg
import functools
from pathlib import Path
from typing import Any, Callable, Optional

log = logging.getLogger(__name__)

_LOG_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
_LOG_DATEFMT = "%d-%m-%Y %H:%M:%S"


#: Used when a caller doesn't name a log file explicitly -- see
#: `add_instrumentation_args`. Log text always has somewhere to go: the
#: console shows only tqdm progress bars (see `cli.py`'s per-frame loops,
#: each wrapped in `tqdm.contrib.logging.logging_redirect_tqdm`), never
#: log records, so a log file is not optional the way it used to be.
DEFAULT_LOG_FILE = Path("chain-cropper.log")


def setup_logging(verbose: bool = False, log_file: Optional[Path] = None) -> None:
    """
    Configure the root logger for a chain_cropper run.

    There is deliberately NO console handler: the console is reserved for
    tqdm progress bars, and a log record printed there would either
    scroll past unread or corrupt a bar's line. All log text goes to
    `log_file` instead, which defaults to `DEFAULT_LOG_FILE` in the
    working directory rather than being silently dropped when the caller
    doesn't name one.

    Parameters
    ----------
    verbose: bool.
        If True, the log file gets DEBUG-level detail; otherwise INFO.
    log_file: Optional[Path].
        Where to write the log. Defaults to `DEFAULT_LOG_FILE` -- pass
        this explicitly (e.g. from `--log-file`) to choose another path,
        not `None` to disable logging altogether, since there is nowhere
        else for the text to go.
    """
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)

    # Remove any handlers a previous call installed, so repeated calls
    # (e.g. across tests) don't duplicate output.
    for handler in list(root.handlers):
        root.removeHandler(handler)

    formatter = logging.Formatter(_LOG_FORMAT, datefmt=_LOG_DATEFMT)

    # Verbose means "tell me what THIS pipeline is doing" in the log file.
    # Without this the root logger at DEBUG also turns on MDAnalysis'
    # internal tracing, which emits a line per topology attribute per
    # Universe and buries everything worth reading.
    for noisy in ("MDAnalysis", "matplotlib", "PIL", "joblib"):
        logging.getLogger(noisy).setLevel(logging.WARNING)

    log_file = Path(log_file) if log_file is not None else DEFAULT_LOG_FILE
    log_file.parent.mkdir(parents=True, exist_ok=True)
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    file_handler.setFormatter(formatter)
    root.addHandler(file_handler)


def install_debug_excepthook(enabled: bool = True) -> None:
    """
    Drop into pdb.post_mortem() on an unhandled exception instead of just
    printing a traceback and exiting.

    Parameters
    ----------
    enabled: bool.
        If False, restore the default excepthook (no-op if it was never
        replaced).
    """
    if not enabled:
        sys.excepthook = sys.__excepthook__
        return

    def _debug_excepthook(exc_type, exc_value, exc_traceback):
        if issubclass(exc_type, KeyboardInterrupt):
            sys.__excepthook__(exc_type, exc_value, exc_traceback)
            return
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        log.debug("Unhandled exception -- dropping into pdb.post_mortem()")
        pdb.post_mortem(exc_traceback)

    sys.excepthook = _debug_excepthook


def add_instrumentation_args(parser: arg.ArgumentParser) -> arg.ArgumentParser:
    """
    Add the standard -v/--verbose, --pdb/--debug, and --log-file flags to
    an argparse parser. Call `apply_instrumentation_args(args)` on the
    parsed namespace to actually wire them up.
    """
    parser.add_argument("-v", "--verbose", action="store_true",
                         help="enable DEBUG-level detail in the log file")
    parser.add_argument("--pdb", "--debug", dest="pdb", action="store_true",
                         help="drop into pdb.post_mortem() on an unhandled exception")
    parser.add_argument("--log-file", type=Path, default=DEFAULT_LOG_FILE,
                         help="where to write the run log (default: %(default)s). "
                              "The console shows only progress bars, never log text.")
    return parser


def apply_instrumentation_args(args: arg.Namespace) -> None:
    """Wire up logging/debugging from a namespace populated via `add_instrumentation_args`."""
    setup_logging(verbose=args.verbose, log_file=args.log_file)
    install_debug_excepthook(enabled=args.pdb)


class Timer:
    """
    Context manager that logs how long a block took.

    Usage
    -----
    with Timer("crop_chains"):
        ...
    """

    def __init__(self, label: str, level: Optional[int] = logging.INFO):
        self.label = label
        #: None times the block without logging it -- see TimingSummary.measure.
        self.level = level
        self.elapsed: Optional[float] = None
        self._start: Optional[float] = None

    def __enter__(self) -> "Timer":
        self._start = time.perf_counter()
        if self.level is not None:
            log.log(self.level, "Starting %s", self.label)
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback) -> None:
        self.elapsed = time.perf_counter() - self._start
        if self.level is None:
            return
        if exc_type is None:
            log.log(self.level, "Finished %s in %.3fs", self.label, self.elapsed)
        else:
            log.log(logging.ERROR, "%s failed after %.3fs", self.label, self.elapsed)


def timed(label: Optional[str] = None, level: int = logging.INFO) -> Callable:
    """
    Decorator equivalent of `Timer` for a whole function call.

    Usage
    -----
    @timed()
    def run_couplings(...):
        ...
    """
    def decorator(func: Callable) -> Callable:
        step_label = label or func.__name__

        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            with Timer(step_label, level=level):
                return func(*args, **kwargs)

        return wrapper

    return decorator


class TimingSummary:
    """
    Accumulates named durations across a run and logs a summary table at
    the end -- useful for seeing where time actually goes across a
    multi-step pipeline.
    """

    def __init__(self):
        self._durations: dict = {}

    def record(self, label: str, seconds: float) -> None:
        self._durations[label] = self._durations.get(label, 0.0) + seconds

    def timer(self, label: str, level: int = logging.INFO) -> "_RecordingTimer":
        """Time a block, log start and finish, and add it to the summary."""
        return _RecordingTimer(self, label, level=level)

    def measure(self, label: str) -> "_RecordingTimer":
        """
        Time a block for the summary without logging anything itself.

        Use this when the code inside is already wrapped in `@timed`, so the
        start and finish are announced once rather than twice under the same
        label -- which reads as the pipeline having run the step twice.
        """
        return _RecordingTimer(self, label, level=None)

    def log_summary(self, level: int = logging.INFO) -> None:
        if not self._durations:
            log.log(level, "Timing summary: no steps recorded")
            return
        total = sum(self._durations.values())
        log.log(level, "Timing summary (total %.3fs):", total)
        for label, seconds in sorted(self._durations.items(), key=lambda kv: -kv[1]):
            pct = 100.0 * seconds / total if total > 0 else 0.0
            log.log(level, "  %-30s %10.3fs (%5.1f%%)", label, seconds, pct)


class _RecordingTimer(Timer):
    """A `Timer` that also feeds its elapsed time into a `TimingSummary`."""

    def __init__(self, summary: TimingSummary, label: str,
                 level: Optional[int] = logging.INFO):
        super().__init__(label, level=level)
        self._summary = summary

    def __exit__(self, exc_type, exc_value, exc_traceback) -> None:
        super().__exit__(exc_type, exc_value, exc_traceback)
        if self.elapsed is not None:
            self._summary.record(self.label, self.elapsed)
