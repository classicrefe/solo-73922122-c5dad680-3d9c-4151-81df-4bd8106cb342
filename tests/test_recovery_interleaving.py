"""Regression tests for the claim / publish / startup-recovery handoff.

The two rolling-restart failure modes these guard against:

* Bug 1: B's startup recovery decided from cached snapshots of the publish
  manifest and job table. When A committed a publish between B's read and
  B's delete, B removed the freshly published lane directory — the task
  still reported 'published' but every download 404'd.
* Bug 2: B's recovery treated a freshly claimed job (status 'processing',
  lock file not yet created) as crashed and marked the live task failed.

Interleavings are made deterministic through the recovery event hook plus
direct flock/heartbeat manipulation — no real second process needed.
"""
import gzip
import io
import json
import time

BARCODES = [
    {"sample_id": "s1", "barcode_r1": "AAAAAAAA", "barcode_r2": "CCCCCCCC"},
]

STALE = "2000-01-01T00:00:00+00:00"  # ancient heartbeat: owner truly gone


def _make_job(dbmod, data_dir, job_id, lane, status="queued",
              heartbeat=None, instance=None):
    """Create a job row plus its retained incoming uploads."""
    upload = data_dir / "incoming" / job_id
    upload.mkdir(parents=True, exist_ok=True)
    r1 = upload / "R1.fastq.gz"
    r2 = upload / "R2.fastq.gz"
    with gzip.open(r1, "wb") as g:
        g.write(b"@x\nAAAAAAAATT\n+\nIIIIIIIIII\n")
    with gzip.open(r2, "wb") as g:
        g.write(b"@x\nCCCCCCCCGG\n+\nJJJJJJJJJJ\n")
    (upload / "COMPLETE").touch()
    dbmod.create_job(
        {"job_id": job_id, "lane": lane, "input_hash": "h",
         "barcode_table": json.dumps(BARCODES),
         "r1_incoming": str(r1), "r2_incoming": str(r2),
         "created_at": dbmod.utcnow()},
        BARCODES,
    )
    if status != "queued":
        dbmod.get_conn().execute(
            "UPDATE jobs SET status=?, instance_id=?, heartbeat_at=? "
            "WHERE job_id=?", (status, instance, heartbeat, job_id))
    return r1, r2


# ---------------- bug 2: live claim must not be failed ----------------

def test_freshly_claimed_job_is_not_failed(data_dir):
    """Claimed (status=processing) but lock not yet acquired: hands off."""
    from app import db as dbmod, processing
    _make_job(dbmod, data_dir, "j1", "laneJ1", status="processing",
              heartbeat=dbmod.utcnow(), instance="instA")
    notes = processing.recover_on_startup()
    assert dbmod.get_job("j1")["status"] == "processing"
    assert notes == []


def test_racing_claim_during_recovery_is_not_failed(data_dir):
    """The claim heartbeat lands while B's recovery is mid-decision; the
    atomic stale-heartbeat predicate must make the fail-update a no-op."""
    from app import db as dbmod, processing
    _make_job(dbmod, data_dir, "j2", "laneJ2", status="processing",
              heartbeat=STALE, instance="instA")

    def hook(phase, subject):
        if phase == "mark_interrupted:before" and subject == "j2":
            # Instance A claims/heartbeats the job right here.
            dbmod.heartbeat("j2", "instA")

    notes = processing.recover_on_startup(event_hook=hook)
    assert dbmod.get_job("j2")["status"] == "processing"
    assert not any("j2" in n for n in notes)


# ---------------- truly dead owners are still reclaimed ----------------

def test_truly_dead_job_is_failed_and_cleaned(data_dir):
    """A genuinely crashed owner is failed, its scratch removed, and its
    original uploads kept so the lane can be retried."""
    from app import db as dbmod, processing
    r1, r2 = _make_job(dbmod, data_dir, "j3", "laneJ3", status="processing",
                       heartbeat=STALE, instance="ghost")
    work = data_dir / "work" / "j3"
    (work / "output").mkdir(parents=True)
    (work / "lock").write_text("ghost\n")
    tmp_stage = data_dir / "published" / "lanes" / ".tmp.j3"
    tmp_stage.mkdir(parents=True)
    final = data_dir / "published" / "lanes" / "laneJ3"
    final.mkdir(parents=True)
    (final / "s1.r1.fastq.gz").write_bytes(b"partial")

    notes = processing.recover_on_startup()
    row = dbmod.get_job("j3")
    assert row["status"] == "failed"
    assert "interrupted" in row["error"]
    assert not work.exists()
    assert not tmp_stage.exists()
    assert not final.exists()
    # Original uploads survive so the lane can be retried from its input.
    assert r1.is_file() and r2.is_file()
    assert dbmod.requeue_job("j3") is True
    assert dbmod.get_job("j3")["status"] == "queued"
    # A follow-up recovery has nothing left to do.
    assert processing.recover_on_startup() == []


# ---------------- bug 1: publish handoff ----------------

def test_publish_committing_mid_recovery_preserves_lane(data_dir):
    """The rolling-restart race: A moved outputs into the final lane
    directory (still holding the job lock); B's recovery starts and scans;
    A commits the publish while B is mid-scan. DB state and files survive."""
    from app import db as dbmod, processing
    _make_job(dbmod, data_dir, "j4", "laneJ4", status="processing",
              heartbeat=STALE, instance="instA")
    lock = processing.acquire_job_lock("j4")  # A is live: holds the lock
    final = data_dir / "published" / "lanes" / "laneJ4"
    final.mkdir(parents=True)
    payload = gzip.compress(b"@x\nTT\n+\nII\n")
    (final / "s1.r1.fastq.gz").write_bytes(payload)
    (final / "s1.r2.fastq.gz").write_bytes(payload)
    counts = [{
        "sample_id": "s1",
        "r1_file": str(final / "s1.r1.fastq.gz"),
        "r2_file": str(final / "s1.r2.fastq.gz"),
        "read_pairs": 1,
    }]

    def hook(phase, subject):
        if phase == "lane_dir:check" and subject == "laneJ4":
            # A finishes the publish while B's recovery is mid-scan.
            dbmod.publish_job("j4", "laneJ4", counts)
            processing._release(lock)

    notes = processing.recover_on_startup(event_hook=hook)
    assert dbmod.get_job("j4")["status"] == "published"
    assert (final / "s1.r1.fastq.gz").read_bytes() == payload
    assert (final / "s1.r2.fastq.gz").read_bytes() == payload
    assert dbmod.lane_counts("laneJ4")[0]["read_pairs"] == 1
    assert not any("laneJ4" in n for n in notes)


def test_mid_publish_state_without_commit_is_kept(data_dir):
    """A holds the lock, files already in the final dir, commit pending:
    recovery must leave the directory (and staging dir) alone."""
    from app import db as dbmod, processing
    _make_job(dbmod, data_dir, "j5", "laneJ5", status="processing",
              heartbeat=STALE, instance="instA")
    lock = processing.acquire_job_lock("j5")
    try:
        final = data_dir / "published" / "lanes" / "laneJ5"
        final.mkdir(parents=True)
        (final / "s1.r1.fastq.gz").write_bytes(b"data")
        tmp_stage = data_dir / "published" / "lanes" / ".tmp.j5"
        tmp_stage.mkdir(parents=True)
        notes = processing.recover_on_startup()
        assert (final / "s1.r1.fastq.gz").read_bytes() == b"data"
        assert tmp_stage.exists()
        assert (data_dir / "work" / "j5").exists()
        assert dbmod.get_job("j5")["status"] == "processing"
        assert notes == []
    finally:
        processing._release(lock)


def test_published_lane_survives_recovery_without_any_lock(data_dir):
    """A finished publish (lock long released) is untouchable."""
    from app import db as dbmod, processing
    _make_job(dbmod, data_dir, "j7", "laneJ7")
    final = data_dir / "published" / "lanes" / "laneJ7"
    final.mkdir(parents=True)
    (final / "s1.r1.fastq.gz").write_bytes(
        gzip.compress(b"@x\nTT\n+\nII\n"))
    counts = [{
        "sample_id": "s1",
        "r1_file": str(final / "s1.r1.fastq.gz"),
        "r2_file": str(final / "s1.r2.fastq.gz"),
        "read_pairs": 1,
    }]
    dbmod.publish_job("j7", "laneJ7", counts)
    notes = processing.recover_on_startup()
    assert (final / "s1.r1.fastq.gz").is_file()
    assert not any("laneJ7" in n for n in notes)


# ---------------- recovery interrupted again: resume ----------------

def test_interrupted_recovery_resumes_on_next_run(data_dir):
    """Recovery died after marking the job failed but before cleaning up:
    the next run must still clear every remnant."""
    from app import db as dbmod, processing
    _make_job(dbmod, data_dir, "j6", "laneJ6", status="failed",
              heartbeat=STALE, instance=None)
    work = data_dir / "work" / "j6"
    (work / "output").mkdir(parents=True)
    (work / "lock").write_text("ghost\n")
    tmp_stage = data_dir / "published" / "lanes" / ".tmp.j6"
    tmp_stage.mkdir(parents=True)
    final = data_dir / "published" / "lanes" / "laneJ6"
    final.mkdir(parents=True)
    (final / "s1.r1.fastq.gz").write_bytes(b"partial")

    processing.recover_on_startup()
    assert not work.exists()
    assert not tmp_stage.exists()
    assert not final.exists()
    assert processing.recover_on_startup() == []


# ---------------- end to end through the API ----------------

def test_recovery_after_publish_keeps_downloads(client, tmp_path):
    """End to end: publish a lane, run startup recovery as a restarting
    sibling instance would, then downloads must still work."""
    r1 = tmp_path / "R1.fastq.gz"
    r2 = tmp_path / "R2.fastq.gz"
    with gzip.open(r1, "wb") as g:
        g.write(b"@x\nAAAAAAAATT\n+\nIIIIIIIIII\n")
    with gzip.open(r2, "wb") as g:
        g.write(b"@x\nCCCCCCCCGG\n+\nJJJJJJJJJJ\n")
    with open(r1, "rb") as f1, open(r2, "rb") as f2:
        resp = client.post(
            "/jobs",
            data={"lane": "laneRoll", "barcodes": json.dumps(BARCODES)},
            files={"r1": ("R1.fastq.gz", f1, "application/gzip"),
                   "r2": ("R2.fastq.gz", f2, "application/gzip")},
        )
    assert resp.status_code == 201, resp.text
    job_id = resp.json()["job"]["job_id"]
    deadline = time.time() + 15
    row = {}
    while time.time() < deadline:
        row = client.get(f"/jobs/{job_id}").json()
        if row["status"] in ("published", "failed"):
            break
        time.sleep(0.1)
    assert row["status"] == "published"

    # A second instance boots now and runs startup recovery on the volume.
    from app import processing
    processing.recover_on_startup()

    assert client.get("/lanes/laneRoll").status_code == 200
    dl = client.get("/lanes/laneRoll/samples/s1/r1")
    assert dl.status_code == 200
    with gzip.open(io.BytesIO(dl.content)) as g:
        assert g.read().splitlines()[1] == b"TT"  # trimmed payload intact
