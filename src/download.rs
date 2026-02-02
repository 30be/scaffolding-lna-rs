use crate::db::Db;
use anyhow::{Context, Result};
use log::{info, warn, debug};
use rayon::prelude::*;
use std::fs;
use std::io::{Read, Write};
use std::path::Path;
use std::time::Duration;
use rusqlite::params;

const SUMMARY_URL: &str = "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/summary/all/";

pub fn download_summary(path: &Path) -> Result<()> {
    if path.exists() {
        info!("Summary file already exists at {:?}", path);
        return Ok(());
    }
    info!("Downloading summary from {}", SUMMARY_URL);
    let mut response = ureq::get(SUMMARY_URL).call()?.into_body().into_reader();
    let mut file = fs::File::create(path)?;
    std::io::copy(&mut response, &mut file)?;
    Ok(())
}

pub struct Record {
    pub pdb: String,
    pub h_chain: String,
    pub l_chain: String,
    pub resolution: Option<f64>,
    pub species: String,
    pub method: String,
    pub scfv: bool,
}

pub fn parse_summary(path: &Path) -> Result<Vec<Record>> {
    let content = fs::read_to_string(path)?;
    let mut records = Vec::new();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(content.as_bytes());

    for result in reader.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue, 
        };
        
        let pdb = record.get(0).unwrap_or("").to_string();
        let h_chain = record.get(1).unwrap_or("").to_string();
        let l_chain = record.get(2).unwrap_or("").to_string();
        let species = record.get(12).unwrap_or("").to_lowercase();
        let resolution = record.get(13).and_then(|s| s.parse::<f64>().ok());
        let method = record.get(14).unwrap_or("").to_uppercase();
        let scfv = record.get(17).map(|s| s == "True").unwrap_or(false); 

        if species == "homo sapiens" 
           && resolution.map(|r| r <= 3.0).unwrap_or(false)
           && (method.contains("X-RAY") || method.contains("ELECTRON MICROSCOPY"))
           && !scfv 
        {
             records.push(Record {
                pdb,
                h_chain,
                l_chain,
                resolution,
                species,
                method,
                scfv,
            });
        }
    }
    Ok(records)
}

pub fn fetch_pdb(pdb_id: &str) -> Result<String> {
    let url = format!("https://files.rcsb.org/download/{}.pdb", pdb_id);
    let mut body = String::new();
    ureq::get(&url)
        .call()
        .context("Failed to fetch PDB")?
        .into_body()
        .into_reader()
        .read_to_string(&mut body)?;
    Ok(body)
}

pub fn populate_db(db: &mut Db, summary_path: &Path) -> Result<()> {
    download_summary(summary_path)?;
    let records = parse_summary(summary_path)?;
    info!("Found {} valid records after filtering.", records.len());

    // Insert metadata first
    {
        let conn = db.get_conn();
        let mut stmt = conn.prepare(
            "INSERT OR IGNORE INTO antibodies (pdb_id, h_chain, l_chain, resolution, species, method, scfv)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)"
        )?;
        
        conn.execute("BEGIN TRANSACTION", [])?;
        for rec in &records {
            stmt.execute(params![
                rec.pdb,
                rec.h_chain,
                rec.l_chain,
                rec.resolution,
                rec.species,
                rec.method,
                rec.scfv
            ])?;
        }
        conn.execute("COMMIT", [])?;
    }

    // Cleanup: Remove the large summary file as it's now in the DB
    if summary_path.exists() {
        info!("Removing temporary summary file...");
        fs::remove_file(summary_path)?;
    }

    // Identify what needs downloading
    let mut to_download = Vec::new();
    {
        let conn = db.get_conn();
        let mut stmt = conn.prepare("SELECT pdb_id FROM antibodies WHERE pdb_blob IS NULL")?;
        let rows = stmt.query_map([], |row| row.get::<_, String>(0))?;
        for r in rows {
            to_download.push(r?);
        }
    }

    if to_download.is_empty() {
        info!("All PDBs are already downloaded.");
        return Ok(());
    }

    info!("Downloading {} PDBs...", to_download.len());
    
    let chunk_size = 50;
    for chunk in to_download.chunks(chunk_size) {
        let fetched: Vec<(String, Option<String>)> = chunk.par_iter().map(|pdb_id| {
            for _ in 0..3 {
                match fetch_pdb(pdb_id) {
                    Ok(content) => return (pdb_id.clone(), Some(content)),
                    Err(e) => {
                        debug!("Error fetching {}: {}", pdb_id, e);
                        std::thread::sleep(Duration::from_millis(500));
                    }
                }
            }
            warn!("Failed to download {}", pdb_id);
            (pdb_id.clone(), None)
        }).collect();

        let conn = db.get_conn();
        conn.execute("BEGIN TRANSACTION", [])?;
        let mut stmt = conn.prepare("UPDATE antibodies SET pdb_blob = ?1 WHERE pdb_id = ?2")?;
        for (pdb_id, content) in fetched {
            if let Some(c) = content {
                stmt.execute(params![c.as_bytes(), pdb_id])?;
            }
        }
        conn.execute("COMMIT", [])?;
        print!(".");
        std::io::stdout().flush()?;
    }
    println!();
    
    Ok(())
}