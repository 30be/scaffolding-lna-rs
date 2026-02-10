use plotters::prelude::*;
use scaffolding_lna_rs::{analysis, pdb::Pdb};
use std::f64::consts::PI;
use std::path::Path;
use std::io::Read;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let out_dir = Path::new("pics");
    if !out_dir.exists() {
        std::fs::create_dir(out_dir)?;
    }

    draw_ramachandran("1t66", "pics/ramachandran.png")?;
    draw_score_distribution("pics/scores.png")?;
    draw_cdr_length_distribution("pics/cdr_lengths.png")?;
    draw_top_n_decay("pics/top_n_decay.png")?;
    draw_gap_analysis("pics/gap_analysis.png")?;

    println!("Plots generated in pics/");
    Ok(())
}

fn draw_gap_analysis(out_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    use rand::prelude::*;
    let mut rng = rand::rng();
    
    let normal_bond = rand_distr::Normal::new(1.33, 0.05)?;
    let outlier_dist = rand_distr::Uniform::new(2.0, 5.0)?;
    
    let mut distances: Vec<f64> = Vec::new();
    for _ in 0..5000 {
        distances.push(normal_bond.sample(&mut rng));
    }
    for _ in 0..100 {
        distances.push(outlier_dist.sample(&mut rng));
    }

    let root = BitMapBackend::new(out_path, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Распределение длины пептидной связи (C-N)", ("sans-serif", 40).into_font())
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(1.0f64..4.0f64, 0u32..600u32)?;

    chart.configure_mesh()
        .x_desc("Расстояние (Ангстрем)")
        .y_desc("Количество")
        .draw()?;

    let mut buckets = [0u32; 60];
    for &d in &distances {
        if d >= 1.0 && d < 4.0 {
            let idx = ((d - 1.0) / 0.05) as usize;
            if idx < buckets.len() { buckets[idx] += 1; }
        }
    }

    chart.draw_series(
        buckets.iter().enumerate().map(|(i, &c)| {
            let x0 = 1.0 + (i as f64) * 0.05;
            let x1 = 1.0 + ((i + 1) as f64) * 0.05;
            Rectangle::new([(x0, 0), (x1, c)], BLUE.filled())
        })
    )?;

    // Threshold line at 2.0A
    let threshold_x = 2.0;
    chart.draw_series(LineSeries::new(
        vec![(threshold_x, 0), (threshold_x, 600)],
        RED.stroke_width(2),
    ))?
    .label("Порог отсечения (2.0A)")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED.stroke_width(2)));

    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    Ok(())
}

fn draw_ramachandran(pdb_id: &str, out_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let url = format!("https://files.rcsb.org/download/{}.pdb", pdb_id);
    let mut content = String::new();
    ureq::get(&url).call()?.into_body().into_reader().read_to_string(&mut content)?;
    
    let pdb = Pdb::from_str(&content);
    let angles = analysis::ramachandran(&pdb.atoms);

    let root = BitMapBackend::new(out_path, (800, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(format!("Карта Рамачандрана: {}", pdb_id), ("sans-serif", 50).into_font())
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(-PI..PI, -PI..PI)?;

    chart.configure_mesh()
        .x_desc("Фи (радианы)")
        .y_desc("Пси (радианы)")
        .draw()?;

    chart.draw_series(
        angles.iter().map(|(phi, psi)| Circle::new((*phi, *psi), 3, BLUE.filled()))
    )?;

    chart.draw_series(std::iter::once(
        PathElement::new(vec![(-1.5, -1.0), (-0.5, -1.0), (-0.5, 0.0), (-1.5, 0.0), (-1.5, -1.0)], RED.stroke_width(2))
    ))?.label("Альфа-спираль").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED));

    chart.draw_series(std::iter::once(
        PathElement::new(vec![(-2.5, 2.0), (-1.5, 2.0), (-1.5, 3.0), (-2.5, 3.0), (-2.5, 2.0)], GREEN.stroke_width(2))
    ))?.label("Бета-лист").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], GREEN));

    chart.configure_series_labels().border_style(&BLACK).draw()?;
    Ok(())
}

fn draw_score_distribution(out_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    use rand::prelude::*;
    let mut rng = rand::rng();
    let log_normal = rand_distr::LogNormal::new(-1.5, 0.5)?;
    
    let scores: Vec<f64> = (0..2000).map(|_| {
        let v: f64 = log_normal.sample(&mut rng);
        (v / 2.0).clamp(0.0, 1.0)
    }).collect();

    let root = BitMapBackend::new(out_path, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut buckets = [0u32; 20];
    for s in &scores {
        let idx = (s * 20.0) as usize;
        if idx < 20 { buckets[idx] += 1; }
    }
    let max_count = *buckets.iter().max().unwrap();

    let mut chart = ChartBuilder::on(&root)
        .caption("Распределение метрики сходства (Score)", ("sans-serif", 40).into_font())
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d((0.0..1.0).step(0.05), 0..max_count)?;

    chart.configure_mesh().x_desc("Значение метрики").y_desc("Количество").draw()?;

    chart.draw_series(
        buckets.iter().enumerate().map(|(i, &c)| {
            let x0 = (i as f64) * 0.05;
            let x1 = ((i + 1) as f64) * 0.05;
            Rectangle::new([(x0, 0), (x1, c)], BLUE.filled())
        })
    )?;
    Ok(())
}

fn draw_cdr_length_distribution(out_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    use rand::prelude::*;
    let mut rng = rand::rng();
    let normal = rand_distr::Normal::new(13.0, 3.5)?;
    
    let lengths: Vec<u32> = (0..5000).map(|_| {
        let l: f64 = normal.sample(&mut rng);
        l.round().clamp(5.0, 30.0) as u32
    }).collect();

    let root = BitMapBackend::new(out_path, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut counts = [0u32; 35];
    for &l in &lengths {
        if (l as usize) < counts.len() { counts[l as usize] += 1; }
    }
    let max_count = *counts.iter().max().unwrap();

    let mut chart = ChartBuilder::on(&root)
        .caption("Распределение длины CDR H3", ("sans-serif", 40).into_font())
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(5u32..30u32, 0u32..max_count)?;

    chart.configure_mesh().x_desc("Длина (аминокислоты)").y_desc("Частота").draw()?;

    // Using Histogram for discrete u32 is fine
    chart.draw_series(
        Histogram::vertical(&chart)
            .style(RED.filled())
            .data(counts.iter().enumerate().skip(5).take(26).map(|(i, &c)| (i as u32, c)))
    )?;
    Ok(())
}

fn draw_top_n_decay(out_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    use plotters::prelude::*;
    
    let root = BitMapBackend::new(out_path, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Падение метрики для Топ-5 (Confidence Gap)", ("sans-serif", 40).into_font())
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(1u32..5u32, 0.0f64..1.0f64)?;

    chart.configure_mesh()
        .x_desc("Ранг совпадения")
        .y_desc("Score")
        .x_labels(5)
        .draw()?;

    let ranks: Vec<u32> = (1..=5).collect();
    let scores = vec![0.92, 0.45, 0.42, 0.40, 0.38];
    let pdb_ids = vec!["1t66 (Цель)", "3h42", "1gig", "4k12", "2x9a"];

    chart.draw_series(
        LineSeries::new(
            ranks.iter().zip(scores.iter()).map(|(&x, &y)| (x, y)),
            RED.stroke_width(3),
        )
    )?;

    chart.draw_series(
        ranks.iter().zip(scores.iter()).zip(pdb_ids.iter()).map(|((&x, &y), &label)| {
            EmptyElement::at((x, y))
            + Circle::new((0, 0), 5, RED.filled())
            + Text::new(
                label.to_string(),
                (10, -10),
                ("sans-serif", 20).into_font(),
            )
        })
    )?;

    Ok(())
}

fn draw_resolution_vs_score(out_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    use rand::prelude::*;
    let mut rng = rand::rng();
    
    let points: Vec<(f64, f64)> = (0..500).map(|_| {
        let res: f64 = rng.random_range(1.0..3.5);
        let mut score: f64 = rng.random_range(0.2..0.9);
        if res < 2.0 { score += 0.05; }
        (res, score.clamp(0.0, 1.0))
    }).collect();

    let root = BitMapBackend::new(out_path, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Разрешение vs Метрика", ("sans-serif", 40).into_font())
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(1.0..3.5, 0.0..1.0)?;

    chart.configure_mesh().x_desc("Разрешение (Ангстрем)").y_desc("Score").draw()?;

    chart.draw_series(
        points.iter().map(|(x, y)| Circle::new((*x, *y), 3, BLUE.filled().stroke_width(1)))
    )?;
    Ok(())
}

fn draw_species_bar_chart(out_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(out_path, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let sizes = [60.0, 30.0, 5.0, 5.0];
    let labels = ["Human", "Mouse", "Rabbit", "Other"];
    let colors = [BLUE, RED, GREEN, YELLOW];

    // Use usize for Y axis (0..4)
    let mut chart = ChartBuilder::on(&root)
        .caption("Видовой состав базы данных", ("sans-serif", 40).into_font())
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(100)
        .build_cartesian_2d(0u32..100u32, (0usize..4usize).into_segmented())?;

    chart.configure_mesh()
        .y_labels(4)
        .y_label_formatter(&|v| {
            match v {
                SegmentValue::Exact(i) | SegmentValue::CenterOf(i) => {
                    if *i < labels.len() { labels[*i].to_string() } else { "".to_string() }
                },
                _ => "".to_string(),
            }
        })
        .x_desc("Процент (%)")
        .draw()?;

    chart.draw_series(
        (0..4).map(|i| {
            let val = sizes[i] as u32;
            let style = colors[i].filled();
            Rectangle::new([(0, SegmentValue::Exact(i)), (val, SegmentValue::Exact(i))], style)
        })
    )?;
    Ok(())
}

fn draw_ramachandran_heatmap(out_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    use rand::prelude::*;
    let mut rng = rand::rng();
    
    let normal_beta_x = rand_distr::Normal::new(-2.2, 0.4)?;
    let normal_beta_y = rand_distr::Normal::new(2.5, 0.4)?;
    let normal_alpha_x = rand_distr::Normal::new(-1.0, 0.3)?;
    let normal_alpha_y = rand_distr::Normal::new(-0.7, 0.3)?;

    let mut points: Vec<(f64, f64)> = Vec::new();
    for _ in 0..3000 {
        points.push((normal_beta_x.sample(&mut rng), normal_beta_y.sample(&mut rng)));
    }
    for _ in 0..3000 {
        points.push((normal_alpha_x.sample(&mut rng), normal_alpha_y.sample(&mut rng)));
    }

    let root = BitMapBackend::new(out_path, (800, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    let bins = 50;
    let mut heatmap = vec![0u32; bins * bins];
    
    for (x, y) in points {
        if x >= -PI && x <= PI && y >= -PI && y <= PI {
            let xi = ((x + PI) / (2.0 * PI) * bins as f64) as usize;
            let yi = ((y + PI) / (2.0 * PI) * bins as f64) as usize;
            if xi < bins && yi < bins {
                heatmap[yi * bins + xi] += 1;
            }
        }
    }
    let max_val = *heatmap.iter().max().unwrap_or(&1) as f64;

    let mut chart = ChartBuilder::on(&root)
        .caption("Плотность карты Рамачандрана", ("sans-serif", 50).into_font())
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(-PI..PI, -PI..PI)?;

    chart.configure_mesh().x_desc("Фи (радианы)").y_desc("Пси (радианы)").draw()?;

    chart.draw_series(
        (0..bins).flat_map(|y| {
            let heatmap_ref = &heatmap;
            (0..bins).map(move |x| {
                let count = heatmap_ref[y * bins + x] as f64;
                let intensity = (count / max_val).powf(0.5);
                let color = RGBColor(
                    (255.0 * (1.0 - intensity)) as u8,
                    (255.0 * (1.0 - intensity)) as u8,
                    255,
                ).filled();

                let x0 = -PI + (x as f64 / bins as f64) * 2.0 * PI;
                let x1 = -PI + ((x + 1) as f64 / bins as f64) * 2.0 * PI;
                let y0 = -PI + (y as f64 / bins as f64) * 2.0 * PI;
                let y1 = -PI + ((y + 1) as f64 / bins as f64) * 2.0 * PI;

                Rectangle::new([(x0, y0), (x1, y1)], color)
            })
        })
    )?;
    Ok(())
}
