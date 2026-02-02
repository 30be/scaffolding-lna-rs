use std::str::FromStr;
use std::fmt;
use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn distance(&self, other: &Point) -> f64 {
        ((self.x - other.x).powi(2) + (self.y - other.y).powi(2) + (self.z - other.z).powi(2)).sqrt()
    }

    // Vector operations
    pub fn sub(&self, other: &Point) -> Point {
        Point::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }

    pub fn add(&self, other: &Point) -> Point {
        Point::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }

    pub fn dot(&self, other: &Point) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn cross(&self, other: &Point) -> Point {
        Point::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    pub fn norm(&self) -> f64 {
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt()
    }

    pub fn normalize(&self) -> Point {
        let n = self.norm();
        if n == 0.0 { *self } else { Point::new(self.x / n, self.y / n, self.z / n) }
    }
}

#[derive(Debug, Clone)]
pub struct Atom {
    pub serial: i32,
    pub name: String,
    pub alt_loc: char,
    pub res_name: String,
    pub chain_id: char,
    pub res_seq: i32,
    pub i_code: char,
    pub pos: Point,
    pub occupancy: f64,
    pub temp_factor: f64,
    pub element: String,
}

impl Atom {
    // Parse a standard PDB ATOM/HETATM line
    pub fn from_line(line: &str) -> Option<Self> {
        if !line.starts_with("ATOM") && !line.starts_with("HETATM") {
            return None;
        }
        if line.len() < 54 {
            return None;
        }

        // Fixed column widths according to PDB format
        let serial = line.get(6..11)?.trim().parse().ok()?;
        let name = line.get(12..16)?.trim().to_string();
        let alt_loc = line.chars().nth(16)?;
        let res_name = line.get(17..20)?.trim().to_string();
        let chain_id = line.chars().nth(21)?;
        let res_seq = line.get(22..26)?.trim().parse().ok()?;
        let i_code = line.chars().nth(26)?;
        let x = line.get(30..38)?.trim().parse().ok()?;
        let y = line.get(38..46)?.trim().parse().ok()?;
        let z = line.get(46..54)?.trim().parse().ok()?;
        
        let occupancy = line.get(54..60).and_then(|s| s.trim().parse().ok()).unwrap_or(1.0);
        let temp_factor = line.get(60..66).and_then(|s| s.trim().parse().ok()).unwrap_or(0.0);
        let element = line.get(76..78).map(|s| s.trim().to_string()).unwrap_or_default();

        Some(Atom {
            serial,
            name,
            alt_loc,
            res_name,
            chain_id,
            res_seq,
            i_code,
            pos: Point { x, y, z },
            occupancy,
            temp_factor,
            element,
        })
    }
}

pub struct Pdb {
    pub atoms: Vec<Atom>,
}

impl Pdb {
    pub fn from_str(content: &str) -> Self {
        let atoms = content
            .lines()
            .filter_map(Atom::from_line)
            .collect();
        Self { atoms }
    }

    pub fn get_sequence(&self, chain_id: char) -> String {
        let mut seq = String::new();
        let mut seen_residues = std::collections::HashSet::new();
        
        // Filter by chain
        let mut chain_atoms: Vec<&Atom> = self.atoms.iter()
            .filter(|a| a.chain_id == chain_id)
            .collect();
        
        // Sort by residue sequence and i_code? 
        // Typically atoms are sorted, but we should be robust.
        // PDB residue ordering: res_seq asc, then i_code (A, B, ...).
        // Let's rely on simple iteration order for now if file is standard.
        
        for atom in chain_atoms {
            let key = (atom.res_seq, atom.i_code);
            if !seen_residues.contains(&key) {
                seen_residues.insert(key);
                // Simple 3to1 mapping
                seq.push(three_to_one(&atom.res_name));
            }
        }
        seq
    }

    pub fn validate(&self) -> QualityReport {
        let mut report = QualityReport::default();
        
        // Group by chain
        let mut chains: HashMap<char, Vec<&Atom>> = HashMap::new();
        for atom in &self.atoms {
            chains.entry(atom.chain_id).or_default().push(atom);
        }

        for (chain_id, atoms) in chains {
            // Group by residue
            let mut residues: Vec<Vec<&Atom>> = Vec::new();
            let mut curr_res = Vec::new();
            let mut last_key = (-999, ' ');

            for atom in atoms {
                let key = (atom.res_seq, atom.i_code);
                if key != last_key {
                    if !curr_res.is_empty() {
                        residues.push(curr_res);
                    }
                    curr_res = Vec::new();
                    last_key = key;
                }
                curr_res.push(atom);
            }
            if !curr_res.is_empty() {
                residues.push(curr_res);
            }

            // Check Residues
            for res in &residues {
                let has_n = res.iter().any(|a| a.name == "N");
                let has_ca = res.iter().any(|a| a.name == "CA");
                let has_c = res.iter().any(|a| a.name == "C");
                
                if !has_n || !has_ca || !has_c {
                    report.missing_backbone_residues += 1;
                }
            }

            // Check Gaps (Distance between C_i and N_i+1)
            for i in 0..residues.len().saturating_sub(1) {
                let c_curr = residues[i].iter().find(|a| a.name == "C");
                let n_next = residues[i+1].iter().find(|a| a.name == "N");

                if let (Some(c), Some(n)) = (c_curr, n_next) {
                    let dist = c.pos.distance(&n.pos);
                    // Peptide bond is ~1.33A. If > 2.0A (allowing for some error), it's likely a break.
                    // Or if numbering is not sequential (e.g. 10 -> 12).
                    
                    // Check numbering gap (simplified, ignores insertion codes logic for distance)
                    // If res_seq diff > 1, it's a numbering gap.
                    let seq_diff = residues[i+1][0].res_seq - residues[i][0].res_seq;
                    if seq_diff > 1 {
                        report.numbering_gaps += 1;
                    }

                    // Check geometric gap
                    if dist > 2.5 {
                        report.geometric_gaps += 1;
                    }
                }
            }
        }
        report
    }
}

use serde::Serialize;

#[derive(Debug, Default, Clone, Serialize)]
pub struct QualityReport {
    pub missing_backbone_residues: usize,
    pub numbering_gaps: usize,
    pub geometric_gaps: usize,
}

impl QualityReport {
    pub fn is_pass(&self) -> bool {
        // Strict criteria: No gaps, few missing atoms
        self.geometric_gaps == 0 && self.missing_backbone_residues < 5
    }
}

fn three_to_one(res: &str) -> char {
    match res {
        "ALA" => 'A', "CYS" => 'C', "ASP" => 'D', "GLU" => 'E', "PHE" => 'F',
        "GLY" => 'G', "HIS" => 'H', "ILE" => 'I', "LYS" => 'K', "LEU" => 'L',
        "MET" => 'M', "ASN" => 'N', "PRO" => 'P', "GLN" => 'Q', "ARG" => 'R',
        "SER" => 'S', "THR" => 'T', "VAL" => 'V', "TRP" => 'W', "TYR" => 'Y',
        _ => 'X',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_point_distance() {
        let p1 = Point { x: 0.0, y: 0.0, z: 0.0 };
        let p2 = Point { x: 3.0, y: 4.0, z: 0.0 };
        assert_eq!(p1.distance(&p2), 5.0);
    }

    #[test]
    fn test_atom_parsing() {
        let line = "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N";
        let atom = Atom::from_line(line).unwrap();
        assert_eq!(atom.serial, 1);
        assert_eq!(atom.name, "N");
        assert_eq!(atom.res_name, "ALA");
        assert_eq!(atom.chain_id, 'A');
        assert_eq!(atom.pos.x, 10.0);
    }

    #[test]
    fn test_sequence_extraction() {
        let content = "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N\n\
                       ATOM      2  N   GLY A   2      11.000  10.000  10.000  1.00  0.00           N";
        let pdb = Pdb::from_str(content);
        assert_eq!(pdb.get_sequence('A'), "AG");
    }
    
    #[test]
    fn test_three_to_one() {
        assert_eq!(three_to_one("ALA"), 'A');
        assert_eq!(three_to_one("UNK"), 'X');
    }
}