# Changelog

Started only in version 0.5

## v1.2.3

- Added ICEberg and PHAST results to bacannot json summary file

## v1.2.2

- Added a file size check, only enter bacannot summary creation if results files are not empty
- Fixed mob_suite summary parsing (added .astype(str))
- Moved MGE dict key generation for inside the integron_finder loop

## v1.2.1

- Added a quick fix for plasmid finder results parsing, as results from gram-negative have only one database, but for gram-positive it has >1.

## v0.5

### hotfix

A small fix has been performed, the blast databases now in the scripts are created without the `-parse_seqids` option because it was bringing problem in the subset module (`align2subsetgbk`).