#![allow(non_snake_case)]

use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use nalgebra::{DMatrix, SymmetricEigen};
use serde::Serialize;
use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::{Arc, Mutex};
use std::thread;

// Formatting ----

fn magenta<T: std::fmt::Display>(value: T) -> String {
    format!("\x1b[35m{}\x1b[0m", value)
}

fn magenta_line<T: std::fmt::Display>(value: T) {
    println!("{}", magenta(value));
}

// CLI args ----

#[derive(Parser, Debug)]
#[command(
    name = "plasmid-model-cli",
    version,
    about = "Plasmid-host dynamics model CLI",
    long_about = None
)]
struct Args {
    /// Presence/absence matrix CSV. First column must contain tree tip names.
    #[arg(short, long)]
    matrix: PathBuf,

    /// Newick phylogenetic tree.
    #[arg(short, long)]
    tree: PathBuf,

    /// Number of phylogenetic eigenvectors to retain.
    #[arg(short, long, default_value_t = 10)]
    k_dims: usize,

    /// Number of warmup iterations per chain.
    #[arg(short = 'w', long, default_value_t = 2000)]
    warmup: usize,

    /// Number of post-warmup sampling iterations per chain.
    #[arg(short = 's', long, default_value_t = 2000)]
    samples: usize,

    /// Number of MCMC chains.
    #[arg(long, default_value_t = 4)]
    chains: usize,

    /// CmdStan progress refresh interval.
    #[arg(long, default_value_t = 50)]
    refresh: usize,

    /// Stan adapt_delta.
    #[arg(long, default_value_t = 0.95)]
    adapt_delta: f64,

    /// Stan random seed.
    #[arg(long, default_value_t = 123)]
    seed: u64,

    /// Output directory.
    #[arg(short, long, default_value = "./results")]
    out_dir: PathBuf,

    /// Optional CmdStan installation directory (expects bin/stansummary).
    #[arg(long)]
    cmdstan_dir: Option<PathBuf>,
}

// Stan data ----

#[derive(Serialize)]
struct StanData {
    N: usize,
    B: usize,
    K: usize,
    Y: Vec<Vec<u32>>,
    U: Vec<Vec<f64>>,
    lambda: Vec<f64>,
}

#[derive(Serialize)]
struct StanInit {
    mu_bar: f64,
    sigma_mu: Vec<f64>,
    L_block: Vec<Vec<f64>>,
    z_mu: Vec<f64>,
    sigma_phy: Vec<f64>,
    k: Vec<f64>,
    z_phy_std: Vec<Vec<f64>>,
}

// Newick tree ----

#[derive(Debug, Clone)]
struct Node {
    name: String,
    length: f64,
    parent: Option<usize>,
}

// Main ----

fn main() {
    let args = Args::parse();

    validate_args(&args);

    reset_output_dirs(&args.out_dir);

    let inputs_dir = args.out_dir.join("inputs");
    let raw_dir = args.out_dir.join("raw");
    let summary_dir = args.out_dir.join("summary");

    for dir in [&inputs_dir, &raw_dir, &summary_dir] {
        fs::create_dir_all(dir)
            .unwrap_or_else(|e| panic!("Failed to create {}: {}", dir.display(), e));
    }

    let stan_model = find_stan_model();

    let stansummary = find_stansummary(args.cmdstan_dir.as_deref());

    magenta_line("Preparing data...");
    println!();

    // Matrix input ----

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_path(&args.matrix)
        .unwrap_or_else(|e| {
            panic!(
                "Failed to open matrix '{}': {}",
                args.matrix.display(),
                e
            )
        });

    let headers = rdr
        .headers()
        .unwrap_or_else(|e| panic!("Failed to read matrix headers: {}", e))
        .clone();

    if headers.len() < 2 {
        panic!("Matrix must contain a tip column and at least one plasmid column");
    }

    let plasmid_names: Vec<String> = headers
        .iter()
        .skip(1)
        .map(|x| x.trim().to_string())
        .collect();

    let mut tips = Vec::new();
    let mut y_matrix = Vec::new();

    for (row_idx, result) in rdr.records().enumerate() {
        let record = result
            .unwrap_or_else(|e| panic!("Invalid CSV record on row {}: {}", row_idx + 2, e));

        if record.len() != headers.len() {
            panic!(
                "Matrix row {} has {} fields but header has {}",
                row_idx + 2,
                record.len(),
                headers.len()
            );
        }

        let tip = record[0].trim();

        if tip.is_empty() {
            panic!("Empty tip name on matrix row {}", row_idx + 2);
        }

        if tips.iter().any(|x: &String| x == tip) {
            panic!("Duplicate tip '{}' in matrix", tip);
        }

        tips.push(tip.to_string());

        let mut row = Vec::with_capacity(plasmid_names.len());

        for (col_idx, value) in record.iter().skip(1).enumerate() {
            let value = value.trim();

            let parsed: u32 = value.parse().unwrap_or_else(|_| {
                panic!(
                    "Invalid presence/absence value '{}' at row {}, column {}",
                    value,
                    row_idx + 2,
                    col_idx + 2
                )
            });

            if parsed > 1 {
                panic!(
                    "Presence/absence values must be 0 or 1; found {} at row {}, column {}",
                    parsed,
                    row_idx + 2,
                    col_idx + 2
                );
            }

            row.push(parsed);
        }

        y_matrix.push(row);
    }

    let N = tips.len();
    let B = plasmid_names.len();

    if N == 0 {
        panic!("Matrix contains no tips");
    }

    if B == 0 {
        panic!("Matrix contains no plasmids");
    }

    println!("{}", magenta(format!("  Tips: {}", N)));
    println!("{}", magenta(format!("  Plasmids: {}", B)));

    // Tree input ----

    let newick = fs::read_to_string(&args.tree)
        .unwrap_or_else(|e| panic!("Failed to read tree '{}': {}", args.tree.display(), e));

    let mut tree_nodes = parse_newick(&newick);

    add_terminal_epsilon(&mut tree_nodes);

    let mut tip_to_node = HashMap::new();

    for (i, node) in tree_nodes.iter().enumerate() {
        if !node.name.is_empty() {
            if tip_to_node.insert(node.name.clone(), i).is_some() {
                panic!("Duplicate named node '{}' in Newick tree", node.name);
            }
        }
    }

    for tip in &tips {
        if !tip_to_node.contains_key(tip) {
            panic!("Tip '{}' is present in the matrix but missing from the tree", tip);
        }
    }

    // Covariance ----

    let paths: Vec<Vec<usize>> = tips
        .iter()
        .map(|tip| {
            let node = *tip_to_node
                .get(tip)
                .unwrap_or_else(|| panic!("Tip '{}' missing from tree", tip));

            get_path_from_root(node, &tree_nodes)
        })
        .collect();

    let mut cov_mat = DMatrix::<f64>::zeros(N, N);

    for i in 0..N {
        for j in 0..=i {
            let covariance =
                calculate_covariance(&paths[i], &paths[j], &tree_nodes);

            cov_mat[(i, j)] = covariance;
            cov_mat[(j, i)] = covariance;
        }
    }

    cov_mat = (&cov_mat + cov_mat.transpose()) * 0.5;

    // Match the R preprocessing ----

    let initial_eigen = SymmetricEigen::new(cov_mat.clone());

    let mean_eval =
        initial_eigen.eigenvalues.iter().sum::<f64>() / N as f64;

    if !mean_eval.is_finite() || mean_eval <= 0.0 {
        panic!(
            "Invalid mean covariance eigenvalue: {}",
            mean_eval
        );
    }

    cov_mat /= mean_eval;

    // Eigen decomposition ----

    let eigen = SymmetricEigen::new(cov_mat);

    let mut eigen_pairs: Vec<(f64, Vec<f64>)> = (0..N)
        .map(|i| {
            let value = eigen.eigenvalues[i];

            let vector = eigen
                .eigenvectors
                .column(i)
                .iter()
                .copied()
                .collect();

            (value, vector)
        })
        .collect();

    eigen_pairs.sort_by(|a, b| {
        b.0.partial_cmp(&a.0)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let K = args.k_dims.min(N);

    if K == 0 {
        panic!("--k-dims must be greater than zero");
    }

    let total_variance: f64 =
        eigen_pairs.iter().map(|(value, _)| *value).sum();

    let retained_variance: f64 =
        eigen_pairs.iter().take(K).map(|(value, _)| *value).sum();

    let variance_explained = if total_variance > 0.0 {
        retained_variance / total_variance
    } else {
        0.0
    };

    let lambda: Vec<f64> = eigen_pairs
        .iter()
        .take(K)
        .map(|x| x.0)
        .collect();

    let mut U = vec![vec![0.0; K]; N];

    for (k, (_, eigenvector)) in eigen_pairs.iter().take(K).enumerate() {
        for n in 0..N {
            U[n][k] = eigenvector[n];
        }
    }

    println!(
        "{}",
        magenta(format!(
            "  Retained phylogenetic dimensions: {}",
            K
        ))
    );

    println!(
        "{}",
        magenta(format!(
            "  Variance explained: {:.2}%",
            variance_explained * 100.0
        ))
    );

    let target = 0.99;
    let mut cumulative = 0.0;
    let mut dims_needed = 0usize;

    if total_variance > 0.0 {
        for (i, (value, _)) in eigen_pairs.iter().enumerate() {
            cumulative += *value;

            if cumulative / total_variance >= target {
                dims_needed = i + 1;
                break;
            }
        }

        if dims_needed == 0 {
            dims_needed = N;
        }
    }

    println!(
        "{}",
        magenta(format!(
            "  Phylogenetic dimensions needed to explain at least 99% of variance: {}",
            dims_needed
        ))
    );

    println!();

    // Stan data ----

    let data_file = inputs_dir.join("data.json");

    let stan_data = StanData {
        N,
        B,
        K,
        Y: y_matrix,
        U,
        lambda,
    };

    write_json(&data_file, &stan_data);

    // Initial values ----

    let mut L_block = vec![vec![0.0; B]; B];

    for i in 0..B {
        L_block[i][i] = 1.0;
    }

    let init_data = StanInit {
        mu_bar: 0.0,
        sigma_mu: vec![0.5; B],
        L_block,
        z_mu: vec![0.0; B],
        sigma_phy: vec![0.1; B],
        k: vec![0.5; B],
        z_phy_std: vec![vec![0.0; K]; B],
    };

    let init_file = inputs_dir.join("init.json");

    write_json(&init_file, &init_data);

    // MCMC sampling ----

    let data_file = fs::canonicalize(&data_file)
        .unwrap_or_else(|e| panic!("Failed to resolve data file: {}", e));

    let init_file = fs::canonicalize(&init_file)
        .unwrap_or_else(|e| panic!("Failed to resolve init file: {}", e));

    let stan_model = fs::canonicalize(&stan_model)
        .unwrap_or_else(|e| panic!("Failed to resolve CmdStan model: {}", e));

    let samples_prefix = raw_dir.join("cmdstan_samples");
    let log_file = raw_dir.join("mcmc.log");

    let total_iterations = args.warmup + args.samples;
    let total_progress = total_iterations * args.chains;

    println!("{}", magenta("Running CmdStan sampling..."));
    println!();

    println!("{}", magenta(format!("  Model: {}", stan_model.display())));
    println!("{}", magenta(format!("  Chains: {}", args.chains)));
    println!("{}", magenta(format!("  Warmup: {}", args.warmup)));
    println!("{}", magenta(format!("  Sampling: {}", args.samples)));
    println!(
        "{}",
        magenta(format!("  Total iterations: {}", total_iterations))
    );
    println!("{}", magenta(format!("  Seed: {}", args.seed)));
    println!();

    let progress = ProgressBar::new(total_progress as u64);

    progress.set_style(
        ProgressStyle::with_template(
            "{spinner:.magenta} [{elapsed_precise}] {bar:40.magenta/white} {pos}/{len} ({percent}%) {msg:.magenta}",
        )
        .unwrap()
        .progress_chars("##-"),
    );

    progress.set_message(magenta("sampling"));

    let chain_iterations =
        Arc::new(Mutex::new(vec![0usize; args.chains]));

    // CmdStan command ----

    let mut child = Command::new(&stan_model)
        .arg("sample")
        .arg(format!("num_chains={}", args.chains))
        .arg(format!("num_warmup={}", args.warmup))
        .arg(format!("num_samples={}", args.samples))
        .arg("adapt")
        .arg(format!("delta={}", args.adapt_delta))
        .arg("random")
        .arg(format!("seed={}", args.seed))
        .arg("data")
        .arg(format!("file={}", data_file.display()))
        .arg(format!("init={}", init_file.display()))
        .arg("output")
        .arg(format!("file={}.csv", samples_prefix.display()))
        .arg(format!("refresh={}", args.refresh))
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap_or_else(|e| {
            panic!(
                "Failed to execute CmdStan model '{}': {}",
                stan_model.display(),
                e
            )
        });

    let stdout = child
        .stdout
        .take()
        .expect("Failed to capture CmdStan stdout");

    let stderr = child
        .stderr
        .take()
        .expect("Failed to capture CmdStan stderr");

    let stdout_progress = Arc::clone(&chain_iterations);
    let stdout_bar = progress.clone();

    let stdout_thread = thread::spawn(move || {
        read_cmdstan_stream(
            stdout,
            &stdout_bar,
            &stdout_progress,
            args.chains,
        )
    });

    let stderr_progress = Arc::clone(&chain_iterations);
    let stderr_bar = progress.clone();

    let stderr_thread = thread::spawn(move || {
        read_cmdstan_stream(
            stderr,
            &stderr_bar,
            &stderr_progress,
            args.chains,
        )
    });

    let status = child
        .wait()
        .unwrap_or_else(|e| panic!("Failed while waiting for CmdStan: {}", e));

    let stdout_log = stdout_thread
        .join()
        .unwrap_or_else(|_| panic!("CmdStan stdout reader thread panicked"));

    let stderr_log = stderr_thread
        .join()
        .unwrap_or_else(|_| panic!("CmdStan stderr reader thread panicked"));

    let log = format!(
        "===== STDOUT =====\n{}\n===== STDERR =====\n{}",
        stdout_log, stderr_log
    );

    fs::write(&log_file, log)
        .unwrap_or_else(|e| panic!("Failed to write {}: {}", log_file.display(), e));

    if !status.success() {
        progress.abandon_with_message(magenta("sampling failed"));

        panic!(
            "CmdStan model failed with status {}. See {}",
            status,
            log_file.display()
        );
    }

    progress.set_position(total_progress as u64);
    progress.finish_with_message(magenta("sampling complete"));

    progress.println("");
    println!();
    println!();

    // Raw CmdStan outputs ----

    println!("{}", magenta("Writing raw outputs..."));
    println!();

    let summary_csv = raw_dir.join("stansummary.csv");
    let summary_txt = raw_dir.join("stansummary.txt");

    let chain_files: Vec<PathBuf> = (1..=args.chains)
        .map(|chain| {
            let with_suffix = PathBuf::from(format!(
                "{}_{}.csv",
                samples_prefix.display(),
                chain
            ));

            if with_suffix.exists() {
                return with_suffix;
            }

            if args.chains == 1 {
                let single_chain = PathBuf::from(format!(
                    "{}.csv",
                    samples_prefix.display()
                ));

                if single_chain.exists() {
                    return single_chain;
                }

                panic!(
                    "Expected CmdStan output was not found. Checked '{}' and '{}'",
                    with_suffix.display(),
                    single_chain.display()
                );
            }

            panic!(
                "Expected CmdStan chain output was not found: {}",
                with_suffix.display()
            );
        })
        .collect();

    // stansummary ----

    let mut summary_command = Command::new(&stansummary);

    summary_command
        .arg("--percentiles=2.5,50,97.5")
        .arg(format!("--csv_filename={}", summary_csv.display()));

    for file in &chain_files {
        summary_command.arg(file);
    }

    let summary_output = summary_command
        .output()
        .unwrap_or_else(|e| {
            panic!(
                "Failed to execute stansummary '{}': {}",
                stansummary.display(),
                e
            )
        });

    let summary_stdout = String::from_utf8_lossy(&summary_output.stdout);
    let summary_stderr = String::from_utf8_lossy(&summary_output.stderr);

    fs::write(
        &summary_txt,
        format!(
            "{}{}",
            summary_stdout,
            if summary_stderr.is_empty() {
                String::new()
            } else {
                format!("\n===== STDERR =====\n{}", summary_stderr)
            }
        ),
    )
    .unwrap_or_else(|e| {
        panic!(
            "Failed to write {}: {}",
            summary_txt.display(),
            e
        )
    });

    if !summary_output.status.success() {
        panic!(
            "stansummary failed with status {}:\n{}",
            summary_output.status,
            summary_stderr
        );
    }

    // Posterior summaries ----

    println!("{}", magenta("Writing model summaries..."));
    println!();

    let output_csv = summary_dir.join("plasmid_summary.csv");

    extract_plasmid_summary(
        &summary_csv,
        &output_csv,
        &plasmid_names,
    );

    let phi_csv = summary_dir.join("phi.csv");
    let r_csv = summary_dir.join("r.csv");

    write_cooccurrence_scores(
        &chain_files,
        &phi_csv,
        &r_csv,
        &plasmid_names,
        N,
    );

    println!(
        "{}",
        magenta(format!(
            "Completed successfully. Results saved to {}",
            args.out_dir.display()
        ))
    );
}

// CmdStan paths ----

fn find_stan_model() -> PathBuf {
    let docker_path = PathBuf::from("/app/stan/model");

    if docker_path.exists() {
        return docker_path;
    }

    let local_path = PathBuf::from("./stan/model");

    if local_path.exists() {
        return local_path;
    }

    panic!(
        "Could not find compiled Stan model. Expected '{}' or '{}'.",
        docker_path.display(),
        local_path.display()
    );
}

fn find_stansummary(cmdstan_dir: Option<&Path>) -> PathBuf {
    if let Some(root) = cmdstan_dir {
        let preferred = root.join("bin/stansummary");
        if preferred.exists() {
            return preferred;
        }

        let windows = root.join("bin/stansummary.exe");
        if windows.exists() {
            return windows;
        }

        panic!(
            "--cmdstan-dir was provided, but stansummary was not found at '{}' or '{}'.",
            preferred.display(),
            windows.display()
        );
    }

    let candidates = [
        PathBuf::from("/opt/cmdstan/bin/stansummary"),
        PathBuf::from("/app/cmdstan/bin/stansummary"),
        PathBuf::from("./cmdstan/bin/stansummary"),
        PathBuf::from("./cmdstan/bin/stansummary.exe"),
    ];

    for path in candidates {
        if path.exists() {
            return path;
        }
    }

    panic!(
        "Could not find CmdStan stansummary. Use --cmdstan-dir <path> or provide one of: /opt/cmdstan/bin/stansummary, /app/cmdstan/bin/stansummary, ./cmdstan/bin/stansummary"
    );
}

// Argument validation ----

fn validate_args(args: &Args) {
    if !args.matrix.exists() {
        panic!(
            "Matrix file does not exist: {}",
            args.matrix.display()
        );
    }

    if !args.tree.exists() {
        panic!(
            "Tree file does not exist: {}",
            args.tree.display()
        );
    }

    if args.k_dims == 0 {
        panic!("--k-dims must be greater than zero");
    }

    if args.warmup == 0 {
        panic!("--warmup must be greater than zero");
    }

    if args.samples == 0 {
        panic!("--samples must be greater than zero");
    }

    if args.chains == 0 {
        panic!("--chains must be greater than zero");
    }

    if args.refresh == 0 {
        panic!("--refresh must be greater than zero");
    }

    if !(0.0..1.0).contains(&args.adapt_delta) {
        panic!("--adapt-delta must be > 0 and < 1");
    }

    if args.seed == 0 {
        panic!("--seed must be greater than zero");
    }

    if let Some(cmdstan_dir) = &args.cmdstan_dir {
        if !cmdstan_dir.exists() {
            panic!(
                "--cmdstan-dir does not exist: {}",
                cmdstan_dir.display()
            );
        }
    }
}

// JSON writing ----

fn write_json<T: Serialize>(path: &Path, value: &T) {
    let json = serde_json::to_string(value)
        .unwrap_or_else(|e| {
            panic!(
                "Failed to serialize {}: {}",
                path.display(),
                e
            )
        });

    fs::write(path, json)
        .unwrap_or_else(|e| {
            panic!(
                "Failed to write {}: {}",
                path.display(),
                e
            )
        });
}

// Output directory handling ----

fn reset_output_dirs(out_dir: &Path) {
    for subdir in ["inputs", "raw", "summary"] {
        let dir = out_dir.join(subdir);

        if dir.exists() {
            fs::remove_dir_all(&dir)
                .unwrap_or_else(|e| {
                    panic!(
                        "Failed to clear output directory '{}': {}",
                        dir.display(),
                        e
                    )
                });
        }
    }
}

// CmdStan output parsing ----

fn read_cmdstan_stream<R: Read>(
    reader: R,
    progress: &ProgressBar,
    chain_iterations: &Arc<Mutex<Vec<usize>>>,
    n_chains: usize,
) -> String {
    let reader = BufReader::new(reader);
    let mut output = String::new();

    for line in reader.lines() {
        let line = match line {
            Ok(line) => line,
            Err(_) => continue,
        };

        output.push_str(&line);
        output.push('\n');

        if let Some((chain, iteration, phase)) =
            parse_cmdstan_iteration(&line)
        {
            if chain == 0 || chain > n_chains {
                continue;
            }

            let mut state = chain_iterations
                .lock()
                .expect("Failed to lock chain progress");

            let index = chain - 1;

            if iteration > state[index] {
                state[index] = iteration;

                let total: usize = state.iter().sum();

                progress.set_position(total as u64);

                let active_chains =
                    state.iter().filter(|&&x| x > 0).count();

                progress.set_message(format!(
                    "{}: {}/{} chains",
                    phase,
                    active_chains,
                    n_chains
                ));
            }
        }
    }

    output
}

fn parse_cmdstan_iteration(
    line: &str,
) -> Option<(usize, usize, &'static str)> {
    // Example:
    // Chain [1] Iteration: 8300 / 10000 [ 83%] (Sampling)

    let chain_pos = line.find("Chain [")?;
    let after_chain = &line[chain_pos + 7..];

    let chain_end = after_chain.find(']')?;

    let chain: usize = after_chain[..chain_end]
        .trim()
        .parse()
        .ok()?;

    let iteration_pos = line.find("Iteration:")?;
    let after_iteration =
        line[iteration_pos + 10..].trim_start();

    let iteration_end = after_iteration
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(after_iteration.len());

    let iteration: usize = after_iteration[..iteration_end]
        .parse()
        .ok()?;

    let phase = if line.contains("(Warmup)") {
        "Warmup"
    } else if line.contains("(Sampling)") {
        "Sampling"
    } else {
        "Sampling"
    };

    Some((chain, iteration, phase))
}

// Summary extraction ----

fn extract_plasmid_summary(
    summary_csv: &Path,
    output_csv: &Path,
    plasmid_names: &[String],
) {
    let file = fs::File::open(summary_csv)
        .unwrap_or_else(|e| {
            panic!(
                "Failed to open {}: {}",
                summary_csv.display(),
                e
            )
        });

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .comment(Some(b'#'))
        .from_reader(file);

    let headers = rdr
        .headers()
        .unwrap_or_else(|e| {
            panic!(
                "Failed to read stansummary headers: {}",
                e
            )
        })
        .clone();

    let name_idx = find_column(&headers, &["name"])
        .unwrap_or_else(|| {
            panic!(
                "Could not find 'name' column in stansummary"
            )
        });

    let median_idx =
        find_column(&headers, &["50%", "Median"])
            .unwrap_or_else(|| {
                panic!(
                    "Could not find median column in stansummary. Headers: {:?}",
                    headers
                )
            });

    let p025_idx = find_column(&headers, &["2.5%"])
        .unwrap_or_else(|| {
            panic!(
                "Could not find 2.5% column in stansummary"
            )
        });

    let p975_idx = find_column(&headers, &["97.5%"])
        .unwrap_or_else(|| {
            panic!(
                "Could not find 97.5% column in stansummary"
            )
        });

    let mut summaries: HashMap<
        String,
        (String, String, String),
    > = HashMap::new();

    for result in rdr.records() {
        let record = result
            .unwrap_or_else(|e| {
                panic!(
                    "Failed to read stansummary: {}",
                    e
                )
            });

        let name = record[name_idx].to_string();

        summaries.insert(
            name,
            (
                record[median_idx].to_string(),
                record[p025_idx].to_string(),
                record[p975_idx].to_string(),
            ),
        );
    }

    let mut writer =
        csv::Writer::from_path(output_csv)
            .unwrap_or_else(|e| {
                panic!(
                    "Failed to create {}: {}",
                    output_csv.display(),
                    e
                )
            });

    writer
        .write_record([
            "Plasmid",
            "k_median",
            "k_2.5_CI",
            "k_97.5_CI",
            "sigma_median",
            "sigma_2.5_CI",
            "sigma_97.5_CI",
        ])
        .unwrap();

    for (i, plasmid) in plasmid_names.iter().enumerate() {
        let index = i + 1;

        let k_name = format!("k[{}]", index);
        let sigma_name = format!("sigma_phy[{}]", index);

        let k = summaries
            .get(&k_name)
            .unwrap_or_else(|| {
                panic!(
                    "Parameter '{}' was not found in stansummary",
                    k_name
                )
            });

        let sigma = summaries
            .get(&sigma_name)
            .unwrap_or_else(|| {
                panic!(
                    "Parameter '{}' was not found in stansummary",
                    sigma_name
                )
            });

        writer
            .write_record([
                plasmid,
                &k.0,
                &k.1,
                &k.2,
                &sigma.0,
                &sigma.1,
                &sigma.2,
            ])
            .unwrap();
    }

    writer
        .flush()
        .unwrap_or_else(|e| {
            panic!(
                "Failed to write {}: {}",
                output_csv.display(),
                e
            )
        });
}

fn find_column(
    headers: &csv::StringRecord,
    candidates: &[&str],
) -> Option<usize> {
    candidates.iter().find_map(|candidate| {
        headers
            .iter()
            .position(|x| x == *candidate)
    })
}

// Newick parser ----

fn parse_newick(newick: &str) -> Vec<Node> {
    let mut nodes = vec![Node {
        name: String::new(),
        length: 0.0,
        parent: None,
    }];

    let mut current = 0usize;
    let mut name = String::new();
    let mut length = String::new();
    let mut reading_length = false;
    let mut in_comment = false;

    for c in newick.chars() {
        if in_comment {
            if c == ']' {
                in_comment = false;
            }

            continue;
        }

        match c {
            '[' => {
                in_comment = true;
            }

            '(' => {
                let node = nodes.len();

                nodes.push(Node {
                    name: String::new(),
                    length: 0.0,
                    parent: Some(current),
                });

                current = node;
                name.clear();
                length.clear();
                reading_length = false;
            }

            ',' => {
                finish_node(
                    current,
                    &mut nodes,
                    &name,
                    &length,
                    reading_length,
                );

                let parent = nodes[current]
                    .parent
                    .unwrap_or_else(|| {
                        panic!("Malformed Newick tree")
                    });

                let node = nodes.len();

                nodes.push(Node {
                    name: String::new(),
                    length: 0.0,
                    parent: Some(parent),
                });

                current = node;
                name.clear();
                length.clear();
                reading_length = false;
            }

            ')' => {
                finish_node(
                    current,
                    &mut nodes,
                    &name,
                    &length,
                    reading_length,
                );

                current = nodes[current]
                    .parent
                    .unwrap_or_else(|| {
                        panic!("Malformed Newick tree")
                    });

                name.clear();
                length.clear();
                reading_length = false;
            }

            ':' => {
                nodes[current].name =
                    name.trim().to_string();

                name.clear();
                length.clear();
                reading_length = true;
            }

            ';' => {
                finish_node(
                    current,
                    &mut nodes,
                    &name,
                    &length,
                    reading_length,
                );

                break;
            }

            ' ' | '\n' | '\r' | '\t' => {}

            _ => {
                if reading_length {
                    length.push(c);
                } else {
                    name.push(c);
                }
            }
        }
    }

    nodes
}

fn finish_node(
    node: usize,
    nodes: &mut [Node],
    name: &str,
    length: &str,
    reading_length: bool,
) {
    if reading_length {
        if !length.trim().is_empty() {
            nodes[node].length =
                length.trim().parse::<f64>().unwrap_or_else(|_| {
                    panic!(
                        "Invalid branch length '{}' in Newick tree",
                        length
                    )
                });
        }
    } else if !name.trim().is_empty() {
        nodes[node].name =
            name.trim().to_string();
    }
}

// Add epsilon to terminal branches ----

fn add_terminal_epsilon(nodes: &mut Vec<Node>) {
    let mut has_child = vec![false; nodes.len()];

    for node in &mut *nodes {
        if let Some(parent) = node.parent {
            has_child[parent] = true;
        }
    }

    for i in 0..nodes.len() {
        if !has_child[i] && !nodes[i].name.is_empty() {
            nodes[i].length += 1e-6;
        }
    }
}

// Tree utilities ----

fn get_path_from_root(
    node: usize,
    nodes: &[Node],
) -> Vec<usize> {
    let mut path = Vec::new();
    let mut current = node;

    path.push(current);

    while let Some(parent) = nodes[current].parent {
        current = parent;
        path.push(current);
    }

    path.reverse();

    path
}

fn calculate_covariance(
    path_i: &[usize],
    path_j: &[usize],
    nodes: &[Node],
) -> f64 {
    let mut covariance = 0.0;

    for (a, b) in path_i.iter().zip(path_j.iter()) {
        if a == b {
            covariance += nodes[*a].length;
        } else {
            break;
        }
    }

    covariance
}

// Co-occurrence output ----

fn write_cooccurrence_scores(
    chain_files: &[PathBuf],
    phi_csv: &Path,
    r_csv: &Path,
    plasmid_names: &[String],
    n_tips: usize,
) {
    let n_plasmids = plasmid_names.len();

    let mut phi_draws: Vec<Vec<f64>> =
        vec![Vec::new(); n_plasmids * n_plasmids];

    let mut r_draws: Vec<Vec<f64>> =
        vec![Vec::new(); n_plasmids * n_plasmids];

    for chain_file in chain_files {
        let file = fs::File::open(chain_file)
            .unwrap_or_else(|e| {
                panic!(
                    "Failed to open {}: {}",
                    chain_file.display(),
                    e
                )
            });

        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(true)
            .comment(Some(b'#'))
            .from_reader(file);

        let headers = rdr
            .headers()
            .unwrap_or_else(|e| {
                panic!(
                    "Failed to read {}: {}",
                    chain_file.display(),
                    e
                )
            })
            .clone();

        let phy_cols: Vec<usize> = headers
            .iter()
            .enumerate()
            .filter_map(|(idx, name)| {
                if name.starts_with("phy[")
                    || name.starts_with("phy.")
                {
                    Some(idx)
                } else {
                    None
                }
            })
            .collect();

        if phy_cols.is_empty() {
            panic!(
                "No phy[...] columns were found in {}. The CmdStan model output does not include latent phylogenetic effects.",
                chain_file.display()
            );
        }

        let phy_indices = parse_phy_indices(
            &headers,
            &phy_cols,
            n_tips,
            n_plasmids,
        );

        for result in rdr.records() {
            let record = result
                .unwrap_or_else(|e| {
                    panic!(
                        "Failed to read {}: {}",
                        chain_file.display(),
                        e
                    )
                });

            let mut mat =
                vec![vec![0.0; n_plasmids]; n_tips];

            for (phy_index, (row_idx, col_idx)) in
                phy_indices.iter().enumerate()
            {
                let col = *phy_cols
                    .get(phy_index)
                    .unwrap_or_else(|| {
                        panic!(
                            "Phy column index out of range for {}",
                            chain_file.display()
                        )
                    });

                let raw_value = record[col].trim();

                let value =
                    raw_value.parse::<f64>().unwrap_or_else(|_| {
                        panic!(
                            "Invalid phy value '{}' in {}",
                            raw_value,
                            chain_file.display()
                        )
                    });

                mat[*row_idx][*col_idx] = value;
            }

            let (phi, r) =
                compute_phi_and_r(&mat);

            for i in 0..n_plasmids {
                for j in 0..n_plasmids {
                    let idx =
                        i * n_plasmids + j;

                    phi_draws[idx]
                        .push(phi[i][j]);

                    r_draws[idx]
                        .push(r[i][j]);
                }
            }
        }
    }

    if phi_draws.iter().all(Vec::is_empty) {
        panic!(
            "No posterior draws were found for phylogenetic co-occurrence analysis"
        );
    }

    let phi_matrix =
        signed_tail_probability_matrix(
            &phi_draws,
            n_plasmids,
        );

    let r_matrix =
        signed_tail_probability_matrix(
            &r_draws,
            n_plasmids,
        );

    write_named_matrix_csv(
        phi_csv,
        plasmid_names,
        &phi_matrix,
    );

    write_named_matrix_csv(
        r_csv,
        plasmid_names,
        &r_matrix,
    );
}

fn compute_phi_and_r(
    mat: &[Vec<f64>],
) -> (
    Vec<Vec<f64>>,
    Vec<Vec<f64>>,
) {
    let n_rows = mat.len();
    let n_cols = mat[0].len();

    if n_rows < 2 {
        panic!(
            "At least two tips are required to calculate co-occurrence scores"
        );
    }

    let mut phi =
        vec![vec![0.0; n_cols]; n_cols];

    let mut r =
        vec![vec![0.0; n_cols]; n_cols];

    for i in 0..n_cols {
        let x: Vec<f64> =
            (0..n_rows)
                .map(|row| mat[row][i])
                .collect();

        let mean_x =
            x.iter().sum::<f64>() /
            x.len() as f64;

        for j in 0..n_cols {
            let y: Vec<f64> =
                (0..n_rows)
                    .map(|row| mat[row][j])
                    .collect();

            let mean_y =
                y.iter().sum::<f64>() /
                y.len() as f64;

            let cov =
                (0..n_rows)
                    .map(|row| {
                        (x[row] - mean_x)
                            * (y[row] - mean_y)
                    })
                    .sum::<f64>()
                    / (n_rows - 1) as f64;

            let var_x =
                (0..n_rows)
                    .map(|row| {
                        (x[row] - mean_x).powi(2)
                    })
                    .sum::<f64>()
                    / (n_rows - 1) as f64;

            let var_y =
                (0..n_rows)
                    .map(|row| {
                        (y[row] - mean_y).powi(2)
                    })
                    .sum::<f64>()
                    / (n_rows - 1) as f64;

            let sd_x = var_x.sqrt();
            let sd_y = var_y.sqrt();

            phi[i][j] =
                if sd_x > 0.0 && sd_y > 0.0 {
                    cov / (sd_x * sd_y)
                } else {
                    0.0
                };
        }
    }

    let mut centered =
        vec![vec![0.0; n_cols]; n_rows];

    let row_means: Vec<f64> =
        mat.iter()
            .map(|row| {
                row.iter().sum::<f64>()
                    / row.len() as f64
            })
            .collect();

    for row in 0..n_rows {
        for col in 0..n_cols {
            centered[row][col] =
                mat[row][col] - row_means[row];
        }
    }

    let mut s =
        vec![vec![0.0; n_cols]; n_cols];

    for i in 0..n_cols {
        for j in 0..n_cols {
            s[i][j] =
                (0..n_rows)
                    .map(|row| {
                        centered[row][i]
                            * centered[row][j]
                    })
                    .sum::<f64>()
                    / (n_rows - 1) as f64;
        }
    }

    for i in 0..n_cols {
        for j in 0..n_cols {
            let sd_i =
                s[i][i].sqrt();

            let sd_j =
                s[j][j].sqrt();

            r[i][j] =
                if sd_i > 0.0 && sd_j > 0.0 {
                    s[i][j]
                        / (sd_i * sd_j)
                } else {
                    0.0
                };
        }
    }

    (phi, r)
}

fn signed_tail_probability_matrix(
    draws: &[Vec<f64>],
    n_plasmids: usize,
) -> Vec<Vec<f64>> {
    let mut matrix =
        vec![vec![0.0; n_plasmids]; n_plasmids];

    for i in 0..n_plasmids {
        for j in 0..n_plasmids {
            let idx =
                i * n_plasmids + j;

            let values =
                &draws[idx];

            if values.is_empty() {
                matrix[i][j] = 0.0;
                continue;
            }

            let median = {
                let mut sorted =
                    values.to_vec();

                sorted.sort_by(|a, b| {
                    a.partial_cmp(b)
                        .unwrap_or(
                            std::cmp::Ordering::Equal,
                        )
                });

                if sorted.len() % 2 == 0 {
                    (
                        sorted[sorted.len() / 2 - 1]
                            + sorted[sorted.len() / 2]
                    ) / 2.0
                } else {
                    sorted[sorted.len() / 2]
                }
            };

            let p_pos =
                values
                    .iter()
                    .filter(|&&x| x > 0.0)
                    .count() as f64
                    / values.len() as f64;

            let p_neg =
                values
                    .iter()
                    .filter(|&&x| x < 0.0)
                    .count() as f64
                    / values.len() as f64;

            let sign =
                if median > 0.0 {
                    1.0
                } else if median < 0.0 {
                    -1.0
                } else {
                    0.0
                };

            matrix[i][j] =
                sign * (1.0 - 2.0 * p_pos.min(p_neg));
        }
    }

    matrix
}

fn write_named_matrix_csv(
    path: &Path,
    plasmid_names: &[String],
    matrix: &[Vec<f64>],
) {
    let mut writer =
        csv::Writer::from_path(path)
            .unwrap_or_else(|e| {
                panic!(
                    "Failed to create {}: {}",
                    path.display(),
                    e
                )
            });

    let mut header =
        vec![String::new()];

    header.extend(
        plasmid_names.iter().cloned(),
    );

    writer
        .write_record(&header)
        .unwrap_or_else(|e| {
            panic!(
                "Failed to write header for {}: {}",
                path.display(),
                e
            )
        });

    for (i, row_name) in
        plasmid_names.iter().enumerate()
    {
        let mut record =
            vec![row_name.clone()];

        for j in 0..plasmid_names.len() {
            record.push(
                matrix[i][j].to_string(),
            );
        }

        writer
            .write_record(record)
            .unwrap_or_else(|e| {
                panic!(
                    "Failed to write row {} in {}: {}",
                    i,
                    path.display(),
                    e
                )
            });
    }

    writer
        .flush()
        .unwrap_or_else(|e| {
            panic!(
                "Failed to write {}: {}",
                path.display(),
                e
            )
        });
}

// CmdStan phy column parsing ----

fn parse_phy_indices(
    headers: &csv::StringRecord,
    phy_cols: &[usize],
    n_tips: usize,
    n_plasmids: usize,
) -> Vec<(usize, usize)> {
    let mut parsed =
        Vec::with_capacity(phy_cols.len());

    for idx in phy_cols {
        let name =
            &headers[*idx];

        let (i, j) =
            parse_phy_name(name);

        if i == 0 || i > n_tips {
            panic!(
                "Invalid phy row index {} in '{}'; expected 1..={}",
                i,
                name,
                n_tips
            );
        }

        if j == 0 || j > n_plasmids {
            panic!(
                "Invalid phy column index {} in '{}'; expected 1..={}",
                j,
                name,
                n_plasmids
            );
        }

        parsed.push((
            i - 1,
            j - 1,
        ));
    }

    parsed
}

fn parse_phy_name(
    name: &str,
) -> (usize, usize) {
    if let Some(body) =
        name.strip_prefix("phy[")
            .and_then(|s| s.strip_suffix("]"))
    {
        let mut parts =
            body.split(',');

        let i_str =
            parts.next().unwrap_or_else(|| {
                panic!(
                    "Malformed phy column name '{}': missing row index",
                    name
                )
            });

        let j_str =
            parts.next().unwrap_or_else(|| {
                panic!(
                    "Malformed phy column name '{}': missing column index",
                    name
                )
            });

        if parts.next().is_some() {
            panic!(
                "Malformed phy column name '{}': too many indices",
                name
            );
        }

        let i =
            i_str.trim().parse::<usize>()
                .unwrap_or_else(|_| {
                    panic!(
                        "Invalid phy row index in '{}'",
                        name
                    )
                });

        let j =
            j_str.trim().parse::<usize>()
                .unwrap_or_else(|_| {
                    panic!(
                        "Invalid phy column index in '{}'",
                        name
                    )
                });

        return (i, j);
    }

    let body =
        if let Some(s) =
            name.strip_prefix("phy.")
        {
            s
        } else {
            panic!(
                "Malformed phy column name '{}': expected phy[...] or phy.... format",
                name
            );
        };

    let mut parts =
        body.split('.');

    let i_str =
        parts.next().unwrap_or_else(|| {
            panic!(
                "Malformed phy column name '{}': missing row index",
                name
            )
        });

    let j_str =
        parts.next().unwrap_or_else(|| {
            panic!(
                "Malformed phy column name '{}': missing column index",
                name
            )
        });

    if parts.next().is_some() {
        panic!(
            "Malformed phy column name '{}': too many indices",
            name
        );
    }

    let i =
        i_str.trim().parse::<usize>()
            .unwrap_or_else(|_| {
                panic!(
                    "Invalid phy row index in '{}'",
                    name
                )
            });

    let j =
        j_str.trim().parse::<usize>()
            .unwrap_or_else(|_| {
                panic!(
                    "Invalid phy column index in '{}'",
                    name
                )
            });

    (i, j)
}
