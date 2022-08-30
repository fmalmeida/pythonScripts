from pathlib import Path

def find_files(start_dir, pattern):
    matches = []
    for path in Path(start_dir).rglob(pattern):
        matches.append(path.resolve())
    return matches