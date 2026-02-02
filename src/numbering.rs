use anyhow::{Result, bail};
use log::{warn, debug};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::NamedTempFile;

pub trait NumberingStrategy {
    fn number(&self, sequence: &str, chain_type: &str) -> Result<Vec<(String, String)>>; // (Number, Residue)
}

pub struct AnarciStrategy;

impl AnarciStrategy {
    pub fn new() -> Self {
        Self
    }

    fn find_binary() -> PathBuf {
        // Check local venv first
        let venv_bin = Path::new(".venv/bin/anarcii");
        if venv_bin.exists() {
            return venv_bin.to_path_buf();
        }
        
        // Also check if we are inside a directory that has .venv (e.g. running from root)
        // This is a bit hacky, normally we expect correct env.
        if let Ok(cwd) = std::env::current_dir() {
            let venv_full = cwd.join(".venv/bin/anarcii");
            if venv_full.exists() {
                return venv_full;
            }
        }

        // Default to system PATH
        PathBuf::from("anarcii")
    }
}

impl NumberingStrategy for AnarciStrategy {
    fn number(&self, sequence: &str, _chain_type: &str) -> Result<Vec<(String, String)>> {
        // Create temp fasta
        let mut input_file = NamedTempFile::new()?;
        writeln!(input_file, ">seq\n{}", sequence)?;
        let input_path = input_file.path();
        
        // Output file
        let output_file = NamedTempFile::new()?;
        let _output_path = output_file.path().with_extension("csv"); // ANARCII likely needs suffix or implies it?
        // Actually CLI said -o FILE (must end in .csv)
        // NamedTempFile path usually doesn't end in .csv.
        // We need a temp path that ends in .csv.
        let temp_dir = std::env::temp_dir();
        let output_csv_path = temp_dir.join(format!("anarcii_{}.csv", uuid::Uuid::new_v4()));

        let binary = Self::find_binary();
        debug!("Using ANARCII binary at: {:?}", binary);

        let output = Command::new(binary)
            .arg(input_path)
            .arg("--scheme")
            .arg("martin")
            .arg("-o")
            .arg(&output_csv_path)
            .output();

        let result = match output {
            Ok(o) if o.status.success() => {
                if !output_csv_path.exists() {
                    bail!("ANARCII finished successfully but no output file found.");
                }
                
                let content = std::fs::read_to_string(&output_csv_path)?;
                // Parse CSV
                // Headers: Name,Chain,Score,Query start,Query end,1,2,...
                // Row: seq,H,31.0,0,112,E,V,...
                
                let mut reader = csv::Reader::from_reader(content.as_bytes());
                let headers = reader.headers()?.clone();
                
                // We expect only one record (or one relevant chain if we passed one seq)
                // ANARCII might split chains if it detects multiple domains.
                // For now, take the first valid chain row.
                
                let mut numbered_seq = Vec::new();
                for result in reader.records() {
                    let record = result?;
                    // Iterate columns starting from index 5 (after Name,Chain,Score,Qstart,Qend)
                    // Check headers to be sure.
                    
                    for (i, field) in record.iter().enumerate() {
                        if i < 5 { continue; } // Skip metadata
                        if field == "-" { continue; } // Gap/Missing
                        
                        let number = &headers[i];
                        numbered_seq.push((number.to_string(), field.to_string()));
                    }
                    if !numbered_seq.is_empty() {
                         break; // Found our chain
                    }
                }
                Ok(numbered_seq)
            }
            Ok(o) => {
                 let stderr = String::from_utf8_lossy(&o.stderr);
                 warn!("ANARCII failed: {}", stderr);
                 bail!("ANARCII failed: {}", stderr);
            }
            Err(e) => {
                warn!("Failed to execute ANARCII: {}", e);
                bail!("Failed to execute ANARCII: {}", e);
            }
        };

        // Cleanup
        let _ = std::fs::remove_file(&output_csv_path);
        
        result
    }
}