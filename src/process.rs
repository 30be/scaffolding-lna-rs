use crate::db::Db;
use crate::pdb::Pdb;
use crate::numbering::{AnarciStrategy, NumberingStrategy};
use anyhow::Result;
use log::{info, debug};
use rayon::prelude::*;
use rusqlite::params;
use serde_json::json;

pub fn process_all(db: &mut Db) -> Result<()> {
    info!("Starting processing pipeline...");
    
    // Select unprocessed PDBs
    let mut tasks = Vec::new();
    {
        let conn = db.get_conn();
        let mut stmt = conn.prepare("SELECT pdb_id, pdb_blob, h_chain, l_chain FROM antibodies WHERE processed = FALSE AND pdb_blob IS NOT NULL")?;
        let rows = stmt.query_map([], |row| {
            let id: String = row.get(0)?;
            let blob: Vec<u8> = row.get(1)?;
            let h: String = row.get(2)?;
            let l: String = row.get(3)?;
            Ok((id, blob, h, l))
        })?;
        
        for r in rows {
            tasks.push(r?);
        }
    }

    if tasks.is_empty() {
        info!("Nothing to process.");
        return Ok(());
    }

    info!("Processing {} PDBs...", tasks.len());
    
    let strategy = AnarciStrategy::new();

    let processed_results: Vec<(String, String, usize, usize, bool)> = tasks.par_iter().map(|(id, blob, h_chain, l_chain)| {
        let content = String::from_utf8_lossy(blob);
        let pdb = Pdb::from_str(&content);
        
        // 1. Validation
        let report = pdb.validate();
        let passed_qc = report.is_pass();

        // Extract sequences for chains
        // H_chain field in DB might be "H" or "H,I" etc.
        // We take the first one for MVP simplicity
        let h_id = h_chain.chars().next().unwrap_or('H');
        let l_id = l_chain.chars().next().unwrap_or('L');

        let h_seq = pdb.get_sequence(h_id);
        let l_seq = pdb.get_sequence(l_id);
        
        let mut numbered_h = Vec::new();
        let mut numbered_l = Vec::new();

        // Attempt numbering only if QC passed (optimization)
        if passed_qc {
            if !h_seq.is_empty() {
                 match strategy.number(&h_seq, "antibody") {
                     Ok(res) => numbered_h = res,
                     Err(e) => debug!("Failed to number H chain for {}: {}", id, e),
                 }
            }
            if !l_seq.is_empty() {
                 match strategy.number(&l_seq, "antibody") {
                     Ok(res) => numbered_l = res,
                     Err(e) => debug!("Failed to number L chain for {}: {}", id, e),
                 }
            }
        }

        // Store result as JSON
        let json_meta = json!({
            "status": "processed", 
            "id": id,
            "h_chain_seq": h_seq,
            "l_chain_seq": l_seq,
            "h_numbering": numbered_h,
            "l_numbering": numbered_l,
            "qc": report
        });
        
        (id.clone(), json_meta.to_string(), report.missing_backbone_residues, report.geometric_gaps + report.numbering_gaps, passed_qc)
    }).collect();

    let conn = db.get_conn();
    conn.execute("BEGIN TRANSACTION", [])?;
    let mut stmt = conn.prepare("UPDATE antibodies SET processed = TRUE, json_blob = ?1, missing_backbone = ?2, gaps = ?3, passed_qc = ?4 WHERE pdb_id = ?5")?;
    for (id, json, missing, gaps, passed) in processed_results {
        stmt.execute(params![json, missing as u32, gaps as u32, passed, id])?;
    }
    conn.execute("COMMIT", [])?;

    Ok(())
}