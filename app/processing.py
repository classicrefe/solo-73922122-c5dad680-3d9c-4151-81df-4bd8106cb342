"""Core demultiplexing: streaming paired gzip FASTQ split with crash safety.

Memory profile: every record is handled one at a time using line-based IO;
no whole-file buffers are allocated, so a 512 MiB lane fits in a 128 MiB
container together with the Python runtime.
"""
from __future__ import annotations

import errno
import fcntl
import gzip
import os
import shutil
import time
import zlib
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterator

from . import config, db

QUARANTINE = "__quarantine__"
_VALID_BASES = frozenset(b"ACGTN")


class FastqFormatError(Exception):
    """Raised when either mate is malformed; the whole lane then fails."""


class SizeExceededError(Exception):
    """Raised when decompressed input exceeds the 512 MiB budget."""


@dataclass
class Record:
    name: bytes   # full header line including '@', without newline
    seq: bytes    # bases only
    plus: bytes   # header line without newline
    qual: bytes   # qualities only, raw Phred+33 bytes


def _read_record(gf: gzip.GzipFile, mate: str, index: int) -> Record | None:
    """Read and validate one FASTQ record. Returns None at clean EOF."""
    name = gf.readline()
    if name == b"":
        return None
    seq_line = gf.readline()
    plus_line = gf.readline()
    qual_line = gf.readline()
    if not name.startswith(b"@") or seq_line == b"" or not plus_line.startswith(b"+"):
        raise FastqFormatError(
            f"{mate}: record {index} is not a well-formed FASTQ record")
    if qual_line == b"" or seq_line.count(b"\n") == 0:
        raise FastqFormatError(f"{mate}: record {index} is truncated")

    seq = seq_line.rstrip(b"\r\n")
    qual = qual_line.rstrip(b"\r\n")
    name_token = name.rstrip(b"\r\n")
    if len(seq) == 0:
        raise FastqFormatError(f"{mate}: record {index} has an empty sequence")
    if any(b not in _VALID_BASES for b in seq):
        raise FastqFormatError(
            f"{mate}: record {index} contains characters outside A/C/G/T/N")
    if len(seq) != len(qual):
        raise FastqFormatError(
            f"{mate}: record {index} sequence/quality length mismatch")
    if any(q < 33 or q > 126 for q in qual):
        raise FastqFormatError(
            f"{mate}: record {index} has quality bytes outside Phred+33 range")
    if len(name_token) < 2:
        raise FastqFormatError(f"{mate}: record {index} has an empty read name")
    return Record(name_token, seq, plus_line.rstrip(b"\r\n"), qual)


def _paired_records(r1_path: Path, r2_path: Path,
                    size_budget: int) -> Iterator[tuple[Record, Record]]:
    """Yield validated R1/R2 pairs, enforcing name/order equality and budget."""
    used = 0
    with open(r1_path, "rb") as f1, open(r2_path, "rb") as f2:
        try:
            with gzip.GzipFile(fileobj=f1) as g1, gzip.GzipFile(fileobj=f2) as g2:
                index = 0
                while True:
                    try:
                        r1 = _read_record(g1, "R1", index)
                        r2 = _read_record(g2, "R2", index)
                    except (OSError, EOFError, zlib.error) as exc:
                        raise FastqFormatError(f"gzip stream corrupted: {exc}") from exc
                    if r1 is None and r2 is None:
                        return
                    if r1 is None or r2 is None:
                        raise FastqFormatError(
                            "R1 and R2 contain different numbers of records")
                    n1 = r1.name[1:].split()[0]
                    n2 = r2.name[1:].split()[0]
                    if n1 != n2:
                        raise FastqFormatError(
                            f"record {index}: R1/R2 read names differ "
                            f"({n1.decode('ascii', 'replace')} vs "
                            f"{n2.decode('ascii', 'replace')})")
                    used += (len(r1.name) + len(r1.seq) + len(r1.plus)
                             + len(r1.qual) + 8)
                    used += (len(r2.name) + len(r2.seq) + len(r2.plus)
                             + len(r2.qual) + 8)
                    if used > size_budget:
                        raise SizeExceededError(
                            f"decompressed total exceeds "
                            f"{config.MAX_TOTAL_UNCOMPRESSED} bytes")
                    yield r1, r2
                    index += 1
        except FastqFormatError:
            raise
        except (OSError, EOFError, zlib.error) as exc:
            raise FastqFormatError(f"gzip stream corrupted: {exc}") from exc


@dataclass
class SampleKey:
    sample_id: str
    combined: bytes  # r1 barcode + r2 barcode, 16 bases


def _best_sample(seq1: bytes, seq2: bytes,
                 keys: list[SampleKey]) -> SampleKey | None:
    """Return unique sample within <=1 substitution over both mates, else None."""
    if len(seq1) < config.BARCODE_LEN or len(seq2) < config.BARCODE_LEN:
        return None
    probe = seq1[:config.BARCODE_LEN] + seq2[:config.BARCODE_LEN]
    best: SampleKey | None = None
    best_dist = config.BARCODE_LEN * 2 + 1
    tied = False
    for key in keys:
        ref = key.combined
        dist = 0
        for a, b in zip(probe, ref):
            if a != b:
                dist += 1
                if dist > best_dist:
                    break
        if dist < best_dist:
            best_dist = dist
            best = key
            tied = False
        elif dist == best_dist:
            tied = True
    if best is None or best_dist > 1 or tied:
        return None
    return best


class _Writers:
    """One gzip writer per (bucket, mate); files opened eagerly."""

    def __init__(self, out_dir: Path, sample_ids: list[str]):
        self._files: dict[tuple[str, str], gzip.GzipFile] = {}
        out_dir.mkdir(parents=True, exist_ok=True)
        buckets = sample_ids + [QUARANTINE]
        for bucket in buckets:
            for mate in ("r1", "r2"):
                fh = open(out_dir / f"{bucket}.{mate}.fastq.gz", "wb")
                self._files[(bucket, mate)] = gzip.GzipFile(
                    filename="", mode="wb", fileobj=fh, compresslevel=6)
        self._raw = {k: v.fileobj for k, v in self._files.items()}

    def write(self, bucket: str, mate: str, rec: Record, trimmed: bool) -> None:
        gf = self._files[(bucket, mate)]
        gf.write(b"@" + rec.name[1:] + b"\n")
        if trimmed:
            gf.write(rec.seq[config.BARCODE_LEN:] + b"\n")
            gf.write(b"+\n")
            gf.write(rec.qual[config.BARCODE_LEN:] + b"\n")
        else:
            gf.write(rec.seq + b"\n")
            gf.write(b"+\n")
            gf.write(rec.qual + b"\n")

    def close(self) -> None:
        # GzipFile does not own an externally supplied fileobj: close the
        # gzip layer (flush footer + buffered data), then the raw file.
        for key, gf in self._files.items():
            gf.close()
            self._raw[key].close()


def _work_dir(job_id: str) -> Path:
    return config.WORK_DIR / job_id


def _lock_path(job_id: str) -> Path:
    return _work_dir(job_id) / "lock"


def acquire_job_lock(job_id: str):
    """Exclusive flock held for the whole processing/publish cycle.

    Retries for a bounded time: startup recovery may briefly hold the lock
    of a dead job while it checks/deletes debris, and a requeued job can be
    claimed in that window. The work directory (or lock file) may also have
    just been removed by recovery, so both are recreated as needed.
    """
    deadline = time.monotonic() + config.LOCK_ACQUIRE_SECONDS
    while True:
        wd = _work_dir(job_id)
        try:
            wd.mkdir(parents=True, exist_ok=True)
            fh = open(_lock_path(job_id), "w")
        except OSError:
            # Directory vanished under us (recovery cleanup); recreate.
            if time.monotonic() >= deadline:
                raise
            time.sleep(0.1)
            continue
        try:
            fcntl.flock(fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            fh.close()
            if exc.errno in (errno.EACCES, errno.EAGAIN):
                if time.monotonic() >= deadline:
                    raise RuntimeError(
                        f"job {job_id} is locked by another instance") from exc
                time.sleep(0.1)
                continue
            raise
        fh.write(f"{os.getpid()}\n")
        fh.flush()
        return fh


def process_job(job_row, instance_id: str) -> None:
    """Full lifecycle: split -> stage -> atomic publish -> metadata."""
    job_id = job_row["job_id"]
    lane = job_row["lane"]
    sample_rows = db.get_samples(job_id)
    keys = [
        SampleKey(
            sample_id=r["sample_id"],
            combined=r["barcode_r1"].encode() + r["barcode_r2"].encode(),
        )
        for r in sample_rows
    ]
    sample_ids = [r["sample_id"] for r in sample_rows]
    r1_in = Path(job_row["r1_incoming"])
    r2_in = Path(job_row["r2_incoming"])

    lock = acquire_job_lock(job_id)
    work = _work_dir(job_id)
    staged = work / "output"
    tmp_publish = config.PUBLISHED_DIR / "lanes" / f".tmp.{job_id}"
    final_publish = config.PUBLISHED_DIR / "lanes" / lane
    try:
        # Prove liveness right after the claim -> lock handoff so startup
        # recovery on a sibling instance never mistakes us for a dead owner.
        db.heartbeat(job_id, instance_id)
        # Remove leftovers from a previous crashed attempt of this same job.
        if staged.exists():
            shutil.rmtree(staged)
        if tmp_publish.exists():
            shutil.rmtree(tmp_publish)
        if final_publish.exists() and db.get_conn().execute(
                "SELECT 1 FROM lane_publishes WHERE lane=?", (lane,)).fetchone() is None:
            shutil.rmtree(final_publish)

        writers = _Writers(staged, sample_ids)
        counts = {s: 0 for s in sample_ids}
        counts[QUARANTINE] = 0
        try:
            for pair_index, (r1, r2) in enumerate(
                    _paired_records(r1_in, r2_in,
                                    config.MAX_TOTAL_UNCOMPRESSED)):
                assigned = _best_sample(r1.seq, r2.seq, keys)
                if assigned is None:
                    bucket = QUARANTINE
                    writers.write(bucket, "r1", r1, trimmed=False)
                    writers.write(bucket, "r2", r2, trimmed=False)
                else:
                    bucket = assigned.sample_id
                    writers.write(bucket, "r1", r1, trimmed=True)
                    writers.write(bucket, "r2", r2, trimmed=True)
                counts[bucket] += 1
                if pair_index % config.HEARTBEAT_FLUSH_EVERY == 0:
                    db.heartbeat(job_id, instance_id)
        finally:
            writers.close()
        db.heartbeat(job_id, instance_id)

        # Every bucket (samples + quarantine), both mates, must exist and be
        # a readable gzip file, even when the bucket received no reads.
        for bucket in sample_ids + [QUARANTINE]:
            for mate in ("r1", "r2"):
                p = staged / f"{bucket}.{mate}.fastq.gz"
                if not p.exists() or p.stat().st_size == 0:
                    raise FastqFormatError(f"missing output file {p.name}")
                try:
                    with gzip.open(p, "rb") as verify:
                        while verify.read(65536):
                            pass
                except (OSError, EOFError, zlib.error) as exc:
                    raise FastqFormatError(
                        f"output {p.name} is not a valid gzip: {exc}") from exc
        db.heartbeat(job_id, instance_id)

        # Atomic same-filesystem move into the published tree.
        tmp_publish.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(staged), str(tmp_publish))
        os.rename(tmp_publish, final_publish)

        count_rows = [
            {
                "sample_id": bucket,
                "r1_file": str(final_publish / f"{bucket}.r1.fastq.gz"),
                "r2_file": str(final_publish / f"{bucket}.r2.fastq.gz"),
                "read_pairs": n,
            }
            for bucket, n in counts.items()
        ]
        db.publish_job(job_id, lane, count_rows)

        # Publish is durable; only now remove inputs and scratch space.
        for p in (r1_in, r2_in):
            try:
                p.unlink()
            except FileNotFoundError:
                pass
        incoming_dir = r1_in.parent
        try:
            incoming_dir.rmdir()
        except OSError:
            pass
        shutil.rmtree(work, ignore_errors=True)
    finally:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        finally:
            lock.close()


# ---------------- startup recovery ----------------
#
# Recovery runs on every instance start, possibly while sibling instances are
# actively claiming, processing and publishing lanes on the same volume. It
# must therefore obey three invariants, all enforced on *freshly re-read*
# state rather than cached snapshots:
#
#   1. A 'processing' job whose flock is held OR whose heartbeat is fresh is
#      live: it is never marked failed and its files are never touched. The
#      fresh-heartbeat rule covers the claim -> lock handoff, where the job
#      is already 'processing' but has not created its lock file yet.
#   2. A lane with a lane_publishes row is a published product: its directory
#      is never deleted, even if the job row or lock looks odd.
#   3. Every delete of a job-associated path is serialized against the job
#      flock: recovery holds the lock across the final state re-check and the
#      delete. A worker holds that lock for its whole process/publish cycle,
#      so a publish can neither be mid-flight nor start while recovery holds
#      it — the check and the delete are atomic with respect to publishing.
#
# Everything uses LOCK_NB (recovery never waits on a live lane) and every
# step is idempotent, so a recovery run that is itself interrupted simply
# resumes on the next start.

RecoveryHook = Callable[[str, str], None]


def _lock_is_free(path: Path) -> bool:
    """True iff nobody currently holds flock on path."""
    try:
        fh = open(path, "r")
    except OSError:
        return True  # missing (or just unlinked) lock file: nothing held
    try:
        try:
            fcntl.flock(fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            if exc.errno in (errno.EACCES, errno.EAGAIN):
                return False
            raise
        fcntl.flock(fh.fileno(), fcntl.LOCK_UN)
        return True
    finally:
        fh.close()


def _hold_job_lock(job_id: str):
    """Hold the job's existing flock for a delete decision.

    Returns False when a live owner holds the lock, None when no lock file
    exists (then nobody can be inside the locked process/publish section —
    a not-yet-locked fresh claim is caught by the heartbeat re-check), else
    an open handle the caller must release. Never creates files: recovery
    must not leave debris of its own.
    """
    path = _lock_path(job_id)
    try:
        fh = open(path, "r")
    except OSError:
        return None
    try:
        fcntl.flock(fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as exc:
        fh.close()
        if exc.errno in (errno.EACCES, errno.EAGAIN):
            return False
        raise
    return fh


def _release(fh) -> None:
    try:
        fcntl.flock(fh.fileno(), fcntl.LOCK_UN)
    finally:
        fh.close()


def _heartbeat_fresh(job_row, now: float) -> bool:
    """True when the job's last heartbeat is inside the recovery grace."""
    hb = job_row["heartbeat_at"]
    if not hb:
        return False
    try:
        ts = datetime.fromisoformat(hb).timestamp()
    except (TypeError, ValueError):
        return False
    return (now - ts) < config.RECOVERY_GRACE_SECONDS


def _job_is_live(job_row, now: float) -> bool:
    """A 'processing' job is live when its flock is held or its heartbeat is
    fresh (the claim -> lock window carries the fresh claim heartbeat)."""
    if not _lock_is_free(_lock_path(job_row["job_id"])):
        return True
    return _heartbeat_fresh(job_row, now)


def _delete_with_job_lock(path: Path, job_id: str, still_safe,
                          notes: list[str], hook: RecoveryHook,
                          note: str) -> None:
    """Delete `path` after a fresh re-check, serialized on the job lock.

    When the job's lock file exists we hold it across the `still_safe()`
    re-check and the delete; because a worker holds that same lock for its
    entire process/publish cycle, the re-check and the delete cannot
    straddle a concurrent publish of this job. When no lock file exists the
    fresh heartbeat/status re-check alone is authoritative (a live owner
    always has either the lock or a fresh claim heartbeat).
    """
    held = _hold_job_lock(job_id)
    if held is False:
        return  # live owner: leave everything alone
    try:
        if not still_safe():
            return
        hook("delete:before", str(path))
        shutil.rmtree(path, ignore_errors=True)
        notes.append(note)
    finally:
        if held is not None:
            _release(held)


def _reap_work_dir(wd: Path, now: float, notes: list[str],
                   hook: RecoveryHook) -> None:
    """Remove a work directory unless a live job owns it right now."""
    hook("work_dir:check", wd.name)
    job = db.get_job(wd.name)
    if job is None:
        # No job row at all: pure orphan. Respect a held lock and give a
        # just-created directory a short grace before reclaiming it.
        if not _lock_is_free(wd / "lock"):
            return
        try:
            age = time.time() - wd.stat().st_mtime
        except OSError:
            return  # vanished (or unstatable) mid-scan: leave it alone
        if age < 5:
            return
        hook("delete:before", str(wd))
        shutil.rmtree(wd, ignore_errors=True)
        notes.append(f"removing orphan work directory {wd.name}")
        return
    if job["status"] == "processing" and _job_is_live(job, now):
        return  # actively processed by a sibling instance

    def still_safe() -> bool:
        fresh = db.get_job(wd.name)
        if fresh is None:
            return True
        # A (re-)claim landing while we decided must stop the delete.
        return not (fresh["status"] == "processing"
                    and _job_is_live(fresh, now))

    _delete_with_job_lock(wd, wd.name, still_safe, notes, hook,
                          f"removing work directory {wd.name}")


def _reap_tmp_staging(child: Path, now: float, notes: list[str],
                      hook: RecoveryHook) -> None:
    """Remove a `.tmp.<job>` publish-staging dir with no live publisher."""
    job_id = child.name[len(".tmp."):]
    hook("tmp_dir:check", job_id)
    job = db.get_job(job_id)
    if job is not None and job["status"] == "processing" \
            and _job_is_live(job, now):
        return  # sibling instance is mid-publish for this job

    def still_safe() -> bool:
        fresh = db.get_job(job_id)
        if fresh is not None and fresh["status"] == "processing" \
                and _job_is_live(fresh, now):
            return False
        return True

    _delete_with_job_lock(child, job_id, still_safe, notes, hook,
                          f"removing unfinished publish staging {child.name}")


def _reap_lane_dir(child: Path, now: float, notes: list[str],
                   hook: RecoveryHook) -> None:
    """Remove a final lane dir only when it is not, and is not becoming,
    a published product."""
    lane = child.name
    hook("lane_dir:check", lane)
    if db.lane_is_published(lane):
        return  # published product: never touched
    job = db.get_job_by_lane(lane)
    if job is not None and job["status"] == "published":
        return  # defensive: product of a published job
    if job is not None and job["status"] == "processing" \
            and _job_is_live(job, now):
        return  # sibling instance is mid-publish for this lane
    job_id = job["job_id"] if job is not None else None

    def still_safe() -> bool:
        # Re-read everything: a publish committing between our first look
        # and the delete must stop it (this is the rolling-restart race).
        if db.lane_is_published(lane):
            return False
        fresh = db.get_job_by_lane(lane)
        if fresh is not None and fresh["status"] == "published":
            return False
        if fresh is not None and fresh["status"] == "processing" \
                and _job_is_live(fresh, now):
            return False
        return True

    if job_id is None:
        if still_safe():
            hook("delete:before", str(child))
            shutil.rmtree(child, ignore_errors=True)
            notes.append(f"removing unpublished lane directory {lane}")
        return
    _delete_with_job_lock(child, job_id, still_safe, notes, hook,
                          f"removing unpublished lane directory {lane}")


def recover_on_startup(event_hook: RecoveryHook | None = None) -> list[str]:
    """Reclaim debris from truly interrupted runs; see invariants above.

    `event_hook(phase, subject)` — when given — is invoked at defined
    handoff points so tests can inject controlled interleavings (e.g. a
    publish committing while recovery is mid-scan).
    """
    notes: list[str] = []

    def hook(phase: str, subject: str = "") -> None:
        if event_hook is not None:
            event_hook(phase, subject)

    conn = db.get_conn()
    now = time.time()
    # Same UTC ISO format as db.utcnow(), so lexicographic comparison in the
    # SQL predicate matches chronological order.
    stale_before = datetime.fromtimestamp(
        now - config.RECOVERY_GRACE_SECONDS, timezone.utc).isoformat()

    hook("start")

    # 1. Jobs stuck in 'processing' whose owner is truly gone. The flock
    #    pre-check is the first authority (a live processor always holds
    #    it); the atomic stale-heartbeat predicate is the second — a fresh
    #    claim heartbeat makes the UPDATE a no-op, so the claim -> lock
    #    handoff window can never be misjudged as a crash. Scratch removal
    #    is left to steps 2/3 so an interrupted recovery simply resumes.
    for job in conn.execute(
            "SELECT * FROM jobs WHERE status='processing'").fetchall():
        if not _lock_is_free(_lock_path(job["job_id"])):
            continue  # another instance is actively working this job
        hook("mark_interrupted:before", job["job_id"])
        if db.fail_interrupted_job(
                job["job_id"], stale_before,
                "interrupted (owner process disappeared)"):
            notes.append(f"recovering interrupted job {job['job_id']}")
    hook("mark_interrupted:after")

    # 2. Work directories whose owning job is dead, finished or unknown.
    if config.WORK_DIR.exists():
        for wd in sorted(config.WORK_DIR.iterdir()):
            if wd.is_dir():
                _reap_work_dir(wd, now, notes, hook)

    # 3. Unpublished staging/final directories in the published tree.
    lanes_root = config.PUBLISHED_DIR / "lanes"
    if lanes_root.exists():
        for child in sorted(lanes_root.iterdir()):
            if child.name.startswith(".tmp."):
                _reap_tmp_staging(child, now, notes, hook)
            elif child.is_dir():
                _reap_lane_dir(child, now, notes, hook)

    # 4. Incoming debris of jobs that no longer need their uploads.
    #    A directory without a COMPLETE marker is an upload still streaming
    #    (or a freshly aborted one younger than the grace period) and is
    #    never touched, so another instance's in-flight upload is safe.
    keep_incoming = {
        Path(r["r1_incoming"]).parent.name
        for r in conn.execute(
            "SELECT r1_incoming FROM jobs WHERE status IN "
            "('queued','processing','failed')").fetchall()
    }
    grace_seconds = 3600
    if config.INCOMING_DIR.exists():
        now_fs = time.time()
        for d in config.INCOMING_DIR.iterdir():
            if not d.is_dir() or d.name in keep_incoming:
                continue
            if not (d / "COMPLETE").exists():
                # Aborted mid-upload: only reclaim well after any plausible
                # transfer could still be running.
                try:
                    age = now_fs - d.stat().st_mtime
                except OSError:
                    continue  # vanished mid-scan
                if age < grace_seconds:
                    continue
            notes.append(f"removing stale incoming upload {d.name}")
            shutil.rmtree(d, ignore_errors=True)
    hook("done")
    return notes
