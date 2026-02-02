use crate::pdb::{Atom, Point};
use std::f64::consts::PI;

// Helper to calculate torsion angle between 4 points
fn torsion_angle(p1: Point, p2: Point, p3: Point, p4: Point) -> f64 {
    let b1 = p2.sub(&p1);
    let b2 = p3.sub(&p2);
    let b3 = p4.sub(&p3);

    let n1 = b1.cross(&b2).normalize();
    let n2 = b2.cross(&b3).normalize();
    let m1 = n1.cross(&b2.normalize());

    let x = n1.dot(&n2);
    let y = m1.dot(&n2);

    -y.atan2(x) // Returns radians [-PI, PI]
}

pub fn ramachandran(atoms: &[Atom]) -> Vec<(f64, f64)> {
    // We need to group atoms by residue to find N, CA, C
    // Then iterate sliding window of residues
    let mut residues = Vec::new();
    
    // Simple grouping by res_seq. 
    // Assumes atoms are sorted. A robust way is a HashMap or just strict iteration.
    // Let's assume standard PDB sorting.
    let mut current_res = Vec::new();
    let mut last_seq = -999;
    
    for atom in atoms {
        if atom.res_seq != last_seq {
            if !current_res.is_empty() {
                residues.push(current_res);
            }
            current_res = Vec::new();
            last_seq = atom.res_seq;
        }
        current_res.push(atom.clone());
    }
    if !current_res.is_empty() {
        residues.push(current_res);
    }

    let mut angles = Vec::new();

    for i in 1..residues.len().saturating_sub(1) {
        let prev = &residues[i-1];
        let curr = &residues[i];
        let next = &residues[i+1];

        // Find necessary atoms: C(prev), N(curr), CA(curr), C(curr), N(next)
        let c_prev = prev.iter().find(|a| a.name == "C");
        let n_curr = curr.iter().find(|a| a.name == "N");
        let ca_curr = curr.iter().find(|a| a.name == "CA");
        let c_curr = curr.iter().find(|a| a.name == "C");
        let n_next = next.iter().find(|a| a.name == "N");

        if let (Some(cp), Some(n), Some(ca), Some(c), Some(nn)) = (c_prev, n_curr, ca_curr, c_curr, n_next) {
            let phi = torsion_angle(cp.pos, n.pos, ca.pos, c.pos);
            let psi = torsion_angle(n.pos, ca.pos, c.pos, nn.pos);
            angles.push((phi, psi));
        }
    }

    angles
}

pub fn ramachandran_score(target: &[(f64, f64)], candidate: &[(f64, f64)]) -> f64 {
    // Simple metric: Mean Squared Difference of angles
    // Problem: Angles are periodic. -PI is close to PI.
    // Distance d = min(|a-b|, 2PI - |a-b|)
    
    if target.is_empty() || candidate.is_empty() {
        return 0.0;
    }

    let len = target.len().min(candidate.len());
    let mut sum_sq = 0.0;

    for i in 0..len {
        let (phi1, psi1) = target[i];
        let (phi2, psi2) = candidate[i];

        let d_phi = (phi1 - phi2).abs();
        let d_phi = if d_phi > PI { 2.0 * PI - d_phi } else { d_phi };

        let d_psi = (psi1 - psi2).abs();
        let d_psi = if d_psi > PI { 2.0 * PI - d_psi } else { d_psi };

        sum_sq += d_phi.powi(2) + d_psi.powi(2);
    }

    // Convert to similarity score [0, 1]
    let mse = sum_sq / len as f64;
    1.0 / (1.0 + mse)
}

#[allow(dead_code)]
pub fn align(s1: &[char], s2: &[char]) -> f64 {
    let gap_open = -11.0;
    let gap_extend = -1.0;
    
    let match_score = |c1: char, c2: char| -> f64 {
        if c1 == c2 { 4.0 } else { -1.0 }
    };

    let n = s1.len();
    let m = s2.len();
    let mut dp = vec![vec![0.0; m + 1]; n + 1];
    
    // Init
    for (i, row) in dp.iter_mut().enumerate().take(n + 1).skip(1) {
        row[0] = gap_open + (i as f64 - 1.0) * gap_extend;
    }
    for (j, val) in dp[0].iter_mut().enumerate().take(m + 1).skip(1) {
        *val = gap_open + (j as f64 - 1.0) * gap_extend;
    }

    for i in 1..=n {
        for j in 1..=m {
            let match_val = dp[i-1][j-1] + match_score(s1[i-1], s2[j-1]);
            let delete = dp[i-1][j] + gap_extend; 
            let insert = dp[i][j-1] + gap_extend;
            dp[i][j] = match_val.max(delete).max(insert);
        }
    }
    
    dp[n][m]
}

pub fn rmsd(atoms1: &[Atom], atoms2: &[Atom]) -> f64 {
    if atoms1.len() != atoms2.len() || atoms1.is_empty() {
        return f64::INFINITY;
    }
    let sum_sq: f64 = atoms1.iter().zip(atoms2.iter())
        .map(|(a, b)| a.pos.distance(&b.pos).powi(2))
        .sum();
    (sum_sq / atoms1.len() as f64).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pdb::Point;

    fn mock_atom(x: f64, y: f64, z: f64) -> Atom {
        Atom {
            serial: 1, name: "CA".into(), alt_loc: ' ', res_name: "ALA".into(),
            chain_id: 'A', res_seq: 1, i_code: ' ',
            pos: Point { x, y, z }, occupancy: 1.0, temp_factor: 0.0, element: "C".into()
        }
    }

    #[test]
    fn test_rmsd_identical() {
        let atoms = vec![mock_atom(0.0, 0.0, 0.0), mock_atom(1.0, 0.0, 0.0)];
        assert_eq!(rmsd(&atoms, &atoms), 0.0);
    }

    #[test]
    fn test_rmsd_offset() {
        let a = vec![mock_atom(0.0, 0.0, 0.0)];
        let b = vec![mock_atom(2.0, 0.0, 0.0)];
        assert_eq!(rmsd(&a, &b), 2.0);
    }

    #[test]
    fn test_align_identical() {
        let seq: Vec<char> = "AAAA".chars().collect();
        assert_eq!(align(&seq, &seq), 16.0);
    }

    #[test]
    fn test_align_mismatch() {
        let s1: Vec<char> = "A".chars().collect();
        let s2: Vec<char> = "B".chars().collect();
        assert_eq!(align(&s1, &s2), -1.0);
    }

    #[test]
    fn test_torsion_angle() {
        let p1 = Point::new(1.0, 0.0, 0.0);
        let p2 = Point::new(0.0, 0.0, 0.0);
        let p3 = Point::new(0.0, 1.0, 0.0);
        let p4 = Point::new(0.0, 1.0, 1.0);
        
        let angle = torsion_angle(p1, p2, p3, p4);
        assert!((angle.abs() - PI/2.0).abs() < 1e-6);
    }
}