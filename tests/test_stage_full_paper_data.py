import gzip
import hashlib
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import stage_full_paper_data
from stage_full_paper_data import (
    byte_ranges,
    decompress_and_characterize,
    download_archive,
    median_from_counts,
    normalized_name,
)


def test_byte_ranges_cover_input_exactly(monkeypatch):
    monkeypatch.setattr(stage_full_paper_data, "MIN_RANGE_BYTES", 1)
    ranges = byte_ranges(101, 4)
    assert ranges == [(0, 25), (26, 51), (52, 77), (78, 100)]
    assert sum(end - start + 1 for start, end in ranges) == 101


def test_multipart_download_assembles_and_verifies_archive(tmp_path, monkeypatch):
    payload = bytes(range(251)) * 20

    class RangeHandler(BaseHTTPRequestHandler):
        def do_GET(self):
            value = self.headers["Range"]
            start_text, end_text = value.removeprefix("bytes=").split("-", 1)
            start, end = int(start_text), int(end_text)
            body = payload[start : end + 1]
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {start}-{end}/{len(payload)}")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, format, *args):
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), RangeHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    monkeypatch.setattr(stage_full_paper_data, "MIN_RANGE_BYTES", 1)
    source = {
        "url": f"http://127.0.0.1:{server.server_port}/reads.fastq.gz",
        "source_bytes": len(payload),
        "source_md5": hashlib.md5(payload).hexdigest(),
    }
    destination = tmp_path / "reads.fastq.gz"
    try:
        result = download_archive(source, destination, connections=4)
    finally:
        server.shutdown()
        server.server_close()
        thread.join()

    assert destination.read_bytes() == payload
    assert result["retrieval"] == "downloaded_multipart"
    assert result["range_connections"] == 4
    assert result["sha256"] == hashlib.sha256(payload).hexdigest()
    assert not destination.with_name(destination.name + ".parts").exists()


def test_normalized_name_strips_pair_suffix_and_comment():
    assert normalized_name(b"@read42/1 comment\n") == b"read42"
    assert normalized_name(b"@read42 1:N:0:1\n") == b"read42"


def test_median_from_length_counts():
    from collections import Counter

    assert median_from_counts(Counter({10: 2, 20: 1}), 3) == 10
    assert median_from_counts(Counter({10: 1, 20: 1}), 2) == 15


def test_decompress_and_characterize_validates_and_records_provenance(tmp_path):
    records = (
        b"@a/1 extra\nAC\n+\nII\n"
        b"@b/1\nTGCA\n+\nJJJJ\n"
    )
    archive = tmp_path / "reads.fastq.gz"
    output = tmp_path / "reads.fastq"
    with gzip.open(archive, "wb") as handle:
        handle.write(records)

    result = decompress_and_characterize(archive, output)

    assert output.read_bytes() == records
    assert result["records"] == 2
    assert result["sequence_bases"] == 6
    assert result["length_summary"] == {
        "minimum": 2,
        "median": 3,
        "maximum": 4,
        "distinct_lengths": 2,
    }
    assert result["first_normalized_id"] == "a"
    assert result["last_normalized_id"] == "b"


def test_existing_fastq_requires_matching_provenance(tmp_path):
    archive = tmp_path / "reads.fastq.gz"
    output = tmp_path / "reads.fastq"
    with gzip.open(archive, "wb") as handle:
        handle.write(b"@a\nA\n+\nI\n")
    result = decompress_and_characterize(archive, output)
    assert decompress_and_characterize(archive, output, result) == result
