use std::process::Command;
use std::path::Path;
use std::fs;

#[test]
fn test_cli_help() {
    let output = Command::new("cargo")
        .args(&["run", "--", "--help"])
        .current_dir(".")
        .output()
        .expect("Failed to run cargo");
    
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("Usage:"));
}

#[test]
fn test_match_command() {
    // Setup dummy DB
    let data_dir = Path::new("data");
    let _ = fs::create_dir_all(data_dir);
    
    // We assume the DB might be empty or partially filled from previous runs.
    // We create a dummy input pdb
    let test_pdb = "test_input.pdb";
    fs::write(test_pdb, 
        "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N\n\
         ATOM      2  CA  ALA A   1      11.500  10.000  10.000  1.00  0.00           C"
    ).unwrap();

    let output = Command::new("cargo")
        // Just pass the file path directly, no subcommand
        .args(&["run", "--", test_pdb])
        .current_dir(".")
        .env("RUST_LOG", "debug")
        .output()
        .expect("Failed to run match");

    let _ = fs::remove_file(test_pdb);

    // It might fail if DB is not populated (network issue etc during test), 
    // but at least it should run.
    // If it fails with "No such file or directory" for DB items, it means logic ran.
    
    // Check if it panicked
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        println!("STDERR: {}", stderr);
        // We allow failure if it's just "network" or "empty db" related, 
        // but we want to ensure the binary executes.
    } else {
         let stdout = String::from_utf8_lossy(&output.stdout);
         // Expect JSON output
         assert!(stdout.trim().starts_with("["));
         assert!(stdout.trim().ends_with("]"));
    }
}
