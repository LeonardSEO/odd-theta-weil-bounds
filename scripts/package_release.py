#!/usr/bin/env python3
"""Create deterministic source archives and verify their contents and checksums."""
from pathlib import Path
import gzip
import hashlib
import io
import json
import tarfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
VERSION = json.loads((ROOT / '.zenodo.json').read_text())['version']
PREFIX = f'odd-theta-weil-bounds-v{VERSION}'
SKIP_DIRS = {'.git', '.lake', '.venv', '__pycache__', 'dist'}
SKIP_FILES = {'MANIFEST.sha256', 'SHA256SUMS', '.DS_Store'}
SKIP_SUFFIXES = {'.pyc', '.aux', '.log', '.out', '.toc', '.gz'}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    required = ['paper.tex', 'paper.pdf', 'supplement.tex', 'supplement.pdf',
                'lean/lean-toolchain', 'lean/lakefile.toml', 'lean/lake-manifest.json',
                'lean/OddThetaWeil.lean', 'lean/Axioms.lean', 'README.md', 'NOTICE',
                'CITATION.cff', '.zenodo.json', 'verification/VERIFICATION_REPORT.md']
    for name in required:
        assert (ROOT / name).is_file(), f'Missing required file: {name}'
    files = []
    for p in sorted(ROOT.rglob('*')):
        rel = p.relative_to(ROOT)
        if any(part in SKIP_DIRS for part in rel.parts):
            continue
        assert not p.is_symlink(), f'Nonportable symlink: {rel}'
        if p.is_file() and p.name not in SKIP_FILES and p.suffix not in SKIP_SUFFIXES:
            files.append(rel.as_posix())
    assert any(n.startswith('lean/OddThetaWeil/') and n.endswith('.lean') for n in files)
    content = {name: (ROOT / name).read_bytes() for name in files}
    manifest = ''.join(f'{digest(data)}  {name}\n' for name, data in content.items())
    (ROOT / 'MANIFEST.sha256').write_text(manifest)
    content['MANIFEST.sha256'] = manifest.encode()
    dist = ROOT / 'dist'
    dist.mkdir(exist_ok=True)
    tar_path = dist / f'{PREFIX}-source.tar.gz'
    zip_path = dist / f'{PREFIX}-source.zip'
    with tar_path.open('wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode='w') as archive:
                for name, data in content.items():
                    entry = tarfile.TarInfo(f'{PREFIX}/{name}')
                    entry.size = len(data)
                    entry.mode = 0o644
                    entry.mtime = 0
                    archive.addfile(entry, io.BytesIO(data))
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, data in content.items():
            entry = zipfile.ZipInfo(f'{PREFIX}/{name}', date_time=(1980, 1, 1, 0, 0, 0))
            entry.external_attr = 0o100644 << 16
            entry.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(entry, data)
    expected = {f'{PREFIX}/{name}': data for name, data in content.items()}
    with tarfile.open(tar_path) as archive:
        actual = {member.name: archive.extractfile(member).read() for member in archive}
        assert actual == expected, 'Tar archive content mismatch'
    with zipfile.ZipFile(zip_path) as archive:
        assert archive.testzip() is None, 'ZIP CRC failure'
        assert {name: archive.read(name) for name in archive.namelist()} == expected
    assets = [ROOT / 'paper.pdf', ROOT / 'supplement.pdf', tar_path, zip_path]
    sums = ''.join(f'{digest(p.read_bytes())}  {p.relative_to(ROOT).as_posix()}\n' for p in assets)
    (ROOT / 'SHA256SUMS').write_text(sums)
    for line in sums.splitlines():
        sha, name = line.split('  ', 1)
        assert digest((ROOT / name).read_bytes()) == sha
    print(f'PASS: {len(content)} files in each source archive; Lean sources and lock included; all asset checksums verified.')
    print('Archives exclude SHA256SUMS to avoid a self-referential archive checksum; publish it alongside them.')


if __name__ == '__main__':
    main()
