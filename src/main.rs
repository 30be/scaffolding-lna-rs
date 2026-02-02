use anyhow::Result;
use clap::Parser;
use std::path::{Path, PathBuf};
use log::info;
use scaffolding_lna_rs::{db, download, process, match_ab};

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Path to the PDB file to match
    #[arg(required = true)]
    input: PathBuf,

        /// Force update of the database
        #[arg(short, long)]
        force_update: bool,
    
        /// Number of top matches to return
        #[arg(short = 'n', long, default_value_t = 5)]
        top_n: usize,
    }
    
    fn main() -> Result<()> {
        env_logger::init();
        let cli = Cli::parse();
        
        let db_path = Path::new("data/antibodies.db");
        if let Some(parent) = db_path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        
        let mut db = db::Db::open(db_path)?;
    
        // Auto-initialization
        let needs_init = !db.is_populated()? || cli.force_update;
        if needs_init {
            info!("Database needs initialization or update...");
            let summary_path = Path::new("data/sabdab_summary_all.tsv");
            download::populate_db(&mut db, summary_path)?;
            process::process_all(&mut db)?;
        }
    
        // Default mode: Match
        let matches = match_ab::find_matches(&mut db, &cli.input, cli.top_n)?;
        println!("{}", serde_json::to_string_pretty(&matches)?);
    
        Ok(())
    }
