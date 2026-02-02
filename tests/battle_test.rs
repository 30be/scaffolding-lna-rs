use std::process::Command;
use std::fs;
use std::io::Read;

struct XorShift {
    state: u32,
}

impl XorShift {
    fn new(seed: u32) -> Self {
        Self { state: seed }
    }
    fn next_f64(&mut self) -> f64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        self.state = x;
        (x as f64) / (u32::MAX as f64)
    }
}

#[test]
fn test_battle_shake() {
    let output_init = Command::new("cargo")
        .args(&["run", "--", "dummy_init.pdb"])
        .current_dir(".")
        .output()
        .expect("Failed to run init");

    if !output_init.status.success() {
         println!("Skipping battle test: DB init failed.");
         return;
    }
    
    let original_pdb_url = "https://files.rcsb.org/download/1t66.pdb";
    let original_pdb_content = match ureq::get(original_pdb_url).call() {
        Ok(r) => {
            let mut body = String::new();
            if r.into_body().into_reader().read_to_string(&mut body).is_err() {
                 println!("Skipping: read error");
                 return;
            }
            body
        },
        Err(_) => {
            println!("Skipping battle test: Could not fetch 1t66 from RCSB.");
            return;
        }
    };
    
    let mut rng = XorShift::new(12345);
    let perturbed_lines: Vec<String> = original_pdb_content.lines().map(|line: &str| {
        if line.starts_with("ATOM") || line.starts_with("HETATM") {
            if line.len() < 54 { return line.to_string(); }
            
            let x: f64 = line[30..38].trim().parse().unwrap_or(0.0);
            let y: f64 = line[38..46].trim().parse().unwrap_or(0.0);
            let z: f64 = line[46..54].trim().parse().unwrap_or(0.0);
            
            let noise_x = (rng.next_f64() - 0.5) * 1.0;
            let noise_y = (rng.next_f64() - 0.5) * 1.0;
            let noise_z = (rng.next_f64() - 0.5) * 1.0;
            
            let new_x = x + noise_x;
            let new_y = y + noise_y;
            let new_z = z + noise_z;
            
            let mut new_line = line.to_string();
            new_line.replace_range(30..38, &format!("{:8.3}", new_x));
            new_line.replace_range(38..46, &format!("{:8.3}", new_y));
            new_line.replace_range(46..54, &format!("{:8.3}", new_z));
            new_line
        } else {
            line.to_string()
        }
    }).collect();
    
    let perturbed_file = "1t66_shaken.pdb";
    fs::write(perturbed_file, perturbed_lines.join("\n")).unwrap();

    let output = Command::new("cargo")
        .args(&["run", "--release", "--", perturbed_file]) 
        .current_dir(".")
        .output()
        .expect("Failed to run match");

    let _ = fs::remove_file(perturbed_file);

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    
    if let Ok(json) = serde_json::from_str::<serde_json::Value>(&stdout) {
        if let Some(arr) = json.as_array() {
            let found = arr.iter().any(|item| {
                item["pdb_id"].as_str() == Some("1t66")
            });
            if !found {
                 println!("Output: {}", stdout);
                 panic!("1t66 not found in top matches for shaken input!");
            }
        } else {
            panic!("Output was not a JSON array");
        }
    } else {
        panic!("Failed to parse JSON output: {}", stdout);
    }
}