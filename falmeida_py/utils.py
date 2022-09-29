from pathlib import Path
import pandas as pd

def find_files(start_dir, pattern):
    matches = []
    for path in Path(start_dir).rglob(pattern):
        matches.append(path.resolve())
    return matches

def load_and_subset_gff(file, col, pattern):
    df = pd.read_csv(
        file, sep='\t', 
        names=[
            "seq", "source", "type", "start", "end", 
            "score", "strand", "frame", "attributes"
        ]
    )

    filter = df[col].str.contains(pattern)

    df = df[filter]
    df.drop_duplicates(inplace=True)

    return df