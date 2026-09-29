"""Detached, resumable SPHEREx retrieval: start, status, or stop."""
from pathlib import Path
import argparse
import collections
import csv
import fcntl
import json
import os
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / 'local_data/spherex/sample'
IDS = 'F01 F02 F03 F05 F07 F08 F09 F10 F11 F12 R02 R03 R04 R05 R06 R07 R08 R09 R10 R11 R12'.split()
STATE = OUT / 'background_state.json'
STOP = OUT / 'background.stop'


def state(**values):
    values.update(pid=os.getpid(), updated_unix=time.time())
    tmp = STATE.with_suffix('.tmp')
    tmp.write_text(json.dumps(values, indent=2) + '\n')
    tmp.replace(STATE)


def is_running():
    with (OUT / 'background.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
    return False


def complete():
    try:
        with (OUT / 'extraction_status.csv').open() as stream:
            rows = {row['id']: row for row in csv.DictReader(stream)}
        return all(rows.get(target, {}).get('complete') == 'True' for target in IDS)
    except (OSError, ValueError):
        return False


def worker():
    with (OUT / 'background.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return
        attempt = 0
        while not STOP.exists():
            attempt += 1
            state(status='downloading', attempt=attempt)
            print(f'\nStarting resumable pass {attempt}', flush=True)
            child = subprocess.Popen([sys.executable, '-u', str(HERE / 'retrieve_spherex_sample.py'), '--ids', *IDS])
            while child.poll() is None:
                if STOP.exists():
                    child.terminate()  # Extractor finishes its batch and writes partial CSVs.
                    child.wait()
                    break
                time.sleep(2)
            if STOP.exists():
                break
            if child.returncode == 0 and complete():
                state(status='complete', attempt=attempt)
                return
            delay = min(60 * 2 ** min(attempt - 1, 4), 900)
            state(status='waiting_to_retry', attempt=attempt,
                  last_exit_code=child.returncode, retry_at_unix=time.time() + delay)
            print(f'Incomplete pass; retrying in {delay} seconds.', flush=True)
            deadline = time.monotonic() + delay
            while time.monotonic() < deadline and not STOP.exists():
                time.sleep(2)
        state(status='stopped', attempt=attempt)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['start', 'status', 'stop', 'worker'])
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.action == 'worker':
        worker()
    elif args.action == 'start':
        if is_running():
            print('Already running.')
            return
        STOP.unlink(missing_ok=True)
        with (OUT / 'background.log').open('a') as log:
            process = subprocess.Popen([sys.executable, '-u', str(Path(__file__).resolve()), 'worker'],
                                       stdin=subprocess.DEVNULL, stdout=log, stderr=log,
                                       start_new_session=True, close_fds=True, cwd=HERE)
        print(f'Started background worker PID {process.pid}; log: {OUT / "background.log"}')
    elif args.action == 'stop':
        STOP.touch()
        print('Stop requested; the current batch will be saved before exit.')
    else:
        counts = collections.Counter()
        for target in IDS:
            for path in (OUT / target).glob('*.json'):
                try:
                    counts[json.loads(path.read_text())['status']] += 1
                except (ValueError, KeyError):
                    counts['unreadable_retry_needed'] += 1
        print(json.dumps(dict(running=is_running(), saved_image_outcomes=dict(counts),
                              last_worker_state=json.loads(STATE.read_text()) if STATE.exists() else None), indent=2))


if __name__ == '__main__':
    main()
