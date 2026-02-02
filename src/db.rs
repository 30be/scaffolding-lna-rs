use rusqlite::{params, Connection, Result};
use std::path::Path;

pub struct Db {
    conn: Connection,
}

impl Db {
    pub fn open(path: impl AsRef<Path>) -> anyhow::Result<Self> {
        let conn = Connection::open(path)?;
        // Enable WAL mode for better concurrency
        conn.pragma_update(None, "journal_mode", "WAL")?;
        Self::init(&conn)?;
        Ok(Self { conn })
    }

    // For testing: in-memory DB
    #[allow(dead_code)]
    pub fn open_in_memory() -> anyhow::Result<Self> {
        let conn = Connection::open_in_memory()?;
        Self::init(&conn)?;
        Ok(Self { conn })
    }

    fn init(conn: &Connection) -> Result<()> {
        conn.execute(
            "CREATE TABLE IF NOT EXISTS meta (
                key TEXT PRIMARY KEY,
                value TEXT
            )",
            [],
        )?;
        conn.execute(
            "CREATE TABLE IF NOT EXISTS antibodies (
                pdb_id TEXT PRIMARY KEY,
                h_chain TEXT,
                l_chain TEXT,
                resolution REAL,
                species TEXT,
                method TEXT,
                scfv BOOLEAN,
                pdb_blob BLOB,
                json_blob TEXT,
                processed BOOLEAN DEFAULT FALSE,
                missing_backbone INT DEFAULT 0,
                gaps INT DEFAULT 0,
                passed_qc BOOLEAN DEFAULT FALSE
            )",
            [],
        )?;
        Ok(())
    }

    pub fn get_conn(&self) -> &Connection {
        &self.conn
    }
    
    pub fn is_populated(&self) -> Result<bool> {
        let count: i64 = self.conn.query_row(
            "SELECT COUNT(*) FROM antibodies WHERE processed = TRUE",
            [],
            |row| row.get(0),
        )?;
        Ok(count > 0)
    }

    #[allow(dead_code)]
    pub fn insert_raw(
        &self,
        pdb_id: &str,
        h_chain: &str,
        l_chain: &str,
        resolution: Option<f64>,
        species: &str,
        method: &str,
        scfv: bool,
    ) -> Result<()> {
        self.conn.execute(
            "INSERT OR IGNORE INTO antibodies (pdb_id, h_chain, l_chain, resolution, species, method, scfv)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)",
            params![pdb_id, h_chain, l_chain, resolution, species, method, scfv],
        )?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_db_init() {
        let db = Db::open_in_memory().unwrap();
        assert!(!db.is_populated().unwrap());
    }

    #[test]
    fn test_insert_and_query() {
        let mut db = Db::open_in_memory().unwrap();
        db.insert_raw("1t66", "H", "L", Some(2.8), "human", "x-ray", false).unwrap();
        
        let conn = db.get_conn();
        let mut stmt = conn.prepare("SELECT pdb_id FROM antibodies").unwrap();
        let mut rows = stmt.query([]).unwrap();
        let id: String = rows.next().unwrap().unwrap().get(0).unwrap();
        assert_eq!(id, "1t66");
    }
}
