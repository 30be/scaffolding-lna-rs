# Scaffolding-LNA (Rust)

Mostly vibecoded, so this description below too. See [./thesis.pdf] for details

---

Fast, parallelized, and database-backed implementation of the antibody scaffolding tool.

## Features

- **SQLite Storage**: Efficient management of PDB structures and metadata.
- **Parallel Processing**: Uses `rayon` for concurrent downloads and analysis.
- **Lean Architecture**: Minimal dependencies (no async runtime), fast compile times.
- **Structured Output**: JSON output for easy integration with Nushell or other tools.

## Prerequisites

- Rust (stable)
- `ANARCII` (устанавливается автоматически в `.venv` при первом запуске или вручную).

## Usage

The tool is designed to be used with a single command. It handles initialization (downloading database, processing) automatically on the first run.

```bash
# Match a target PDB against the database
cargo run -- match input.pdb
```

### Flags

- `-f`, `--force-update`: Force re-downloading and re-processing of the SAbDab database.

## Output

The output is a JSON array of matches, sorted by score (descending).

```json
[
  {
    "pdb_id": "1t66",
    "score": 0.85,
    "method": "X-RAY DIFFRACTION"
  }
]
```

## Developer Notes

See [DOCS.md](DOCS.md) for architectural details.

