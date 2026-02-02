use crate::db::Db;
use crate::pdb::Pdb;
use crate::analysis;
use anyhow::Result;
use rayon::prelude::*;
use serde::Serialize;
use std::path::Path;
use log::info;

#[derive(Serialize)]
pub struct MatchResult {
    pub pdb_id: String,
    pub score: f64,
    pub method: String,
}

pub fn find_matches(db: &mut Db, target_path: &Path, top_n: usize) -> Result<Vec<MatchResult>> {
    let target_content = std::fs::read_to_string(target_path)?;
    let target_pdb = Pdb::from_str(&target_content);
    // Extract target sequence (naive extraction from atoms for MVP)
    // Real implementation would group by residue ID and map 3-letter code to 1-letter.
    // For now, let's assume we have a way to compare.
    // Since we don't have robust sequence extraction in pdb.rs yet, we will compare
    // structural similarity or just dummy score for the skeleton.
    
    // We'll use RMSD on first 100 atoms as a dummy metric if counts match, 
    // or just return 0.0 to show the pipeline works.
    
    // Fetch candidates
    let candidates = {
        let conn = db.get_conn();
        // Only select those that passed QC
        let mut stmt = conn.prepare("SELECT pdb_id, pdb_blob, method FROM antibodies WHERE processed = TRUE AND passed_qc = TRUE AND pdb_blob IS NOT NULL")?;
        let rows = stmt.query_map([], |row| {
            Ok((
                row.get::<_, String>(0)?,
                row.get::<_, Vec<u8>>(1)?,
                row.get::<_, String>(2)?,
            ))
        })?;
        
        let mut res = Vec::new();
        for r in rows { res.push(r?); }
        res
    };

    info!("Matching against {} candidates...", candidates.len());

    let mut results: Vec<MatchResult> = candidates.par_iter().map(|(id, blob, method)| {
        let content = String::from_utf8_lossy(blob);
        let candidate_pdb = Pdb::from_str(&content);
        
        // Metric: RMSD + Ramachandran
        // RMSD
        let limit = target_pdb.atoms.len().min(candidate_pdb.atoms.len()).min(50);
        let rmsd_score = if limit > 0 {
             1.0 / (1.0 + analysis::rmsd(&target_pdb.atoms[0..limit], &candidate_pdb.atoms[0..limit]))
        } else {
            0.0
        };

        // Ramachandran
        let target_rama = analysis::ramachandran(&target_pdb.atoms);
        let cand_rama = analysis::ramachandran(&candidate_pdb.atoms);
        let rama_score = analysis::ramachandran_score(&target_rama, &cand_rama);

        // Weighted sum (50/50 for now)
        let score = 0.5 * rmsd_score + 0.5 * rama_score;

        MatchResult {
            pdb_id: id.clone(),
            score,
            method: method.clone(),
        }
    }).collect();

    // Sort by score descending
    results.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal));
    
    Ok(results.into_iter().take(top_n).collect())
}
