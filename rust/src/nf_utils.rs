use std::{
    collections::HashSet,
    fs::{self, remove_dir, remove_dir_all, remove_file},
    io::{self, Write},
    path::{Path, PathBuf},
    process::{Command, Output},
};

use clap::{ArgAction, Parser};
use globset::{Glob, GlobSetBuilder};
use tempfile::NamedTempFile;
use walkdir::WalkDir;

/// Get name of last run
///
fn get_last_run() -> Result<String, std::io::Error> {
    let out: Output = Command::new("nextflow").arg("log").output()?;
    if out.status.success() {
        let last_run: String = String::from_utf8_lossy(&out.stdout)
            .split("\n")
            .filter(|x| *x != "")
            .into_iter()
            .last()
            .unwrap_or("")
            .to_string();
        let name = last_run
            .split_ascii_whitespace()
            .nth(3)
            .unwrap()
            .to_string();
        Ok(name)
    } else {
        let stderr: String = String::from_utf8_lossy(&out.stderr).to_string();
        panic!("{}", stderr);
    }
}

fn run_files(run_name: &str) -> Result<HashSet<PathBuf>, std::io::Error> {
    let out: Output = Command::new("nextflow").args(["log", run_name]).output()?;
    if out.status.success() {
        let last_run: HashSet<PathBuf> = String::from_utf8_lossy(&out.stdout)
            .split("\n")
            .filter(|x| *x != "")
            .map(|p| Path::new(p).to_path_buf())
            .collect();
        Ok(last_run)
    } else {
        let stderr: String = String::from_utf8_lossy(&out.stderr).to_string();
        panic!("{}", stderr);
    }
}

/// Remove empty directories under `dir`
///
fn remove_empty(dir: &str) -> () {
    let _: Vec<Result<(), std::io::Error>> = WalkDir::new(dir)
        .into_iter()
        .filter_map(|e| e.ok())
        .map(|p| p.into_path())
        .filter(|p| p.is_dir())
        .map(|d| remove_dir(d))
        .collect();
}

/// Get all work entries in work dir, including empty dirs
///
fn all_in_work() -> Result<HashSet<PathBuf>, std::io::Error> {
    let work_files: HashSet<PathBuf> = WalkDir::new("work")
        .max_depth(2)
        .min_depth(2)
        .into_iter()
        .filter_map(|e| e.ok())
        .map(|d| d.into_path())
        .collect();
    Ok(work_files)
}

/// Clean up work directory, optionally keeping the files from the last run
///     or a specified run
///
fn clean_work(keep_last: bool, to_keep: Vec<String>) -> Result<(), std::io::Error> {
    let mut untouched: HashSet<PathBuf> = HashSet::new();
    if keep_last {
        let last_run: String = get_last_run().unwrap();
        let last_files = run_files(&last_run).unwrap();
        untouched.extend(last_files);
    }
    for k in to_keep {
        let run: HashSet<PathBuf> = run_files(&k).unwrap();
        untouched.extend(run);
    }
    let all: HashSet<PathBuf> = all_in_work().unwrap();
    let diff: HashSet<&PathBuf> = all.difference(&untouched).collect();
    for d in diff.iter() {
        println!("{}", d.as_path().to_string_lossy());
    }
    loop {
        let mut input = String::new();
        println!("");
        println!("The above will be deleted. Continue? [Y/N] ");
        io::stdin()
            .read_line(&mut input)
            .expect("Failed to read line");
        input = input.trim().to_lowercase();
        if input == "y" {
            for d in diff.into_iter() {
                println!("Removing {}", d.as_path().to_string_lossy());
                let res = remove_file(d);
                if res.is_err() {
                    let res2 = remove_dir_all(d);
                    if res2.is_err() {
                        println!("Failed to remove {}", d.as_path().to_string_lossy());
                    }
                }
            }
            remove_empty("work");
            break;
        } else if input == "n" {
            println!("Delete files aborted");
            break;
        } else {
            println!("Command not recognized");
        }
    }
    Ok(())
}

fn clean_output(outdir: &str, prefixes: Vec<String>, editor: &str) -> Result<(), std::io::Error> {
    let mut builder: GlobSetBuilder = GlobSetBuilder::new();
    if prefixes.len() > 0 {
        for p in prefixes {
            let pattern = format!("{}*", p);
            builder.add(Glob::new(&pattern).unwrap());
        }
    } else {
        builder.add(Glob::new("*").unwrap());
    }
    let globber = builder.build().unwrap();
    let mut output_file_names: Vec<String> = WalkDir::new(outdir)
        .into_iter()
        .filter_map(|e| e.ok())
        .map(|d| d.into_path())
        .filter(move |p| {
            let fname = p.file_name().unwrap();
            globber.is_match(fname) && p.is_file()
        })
        .map(|p| p.to_str().unwrap_or("").to_string())
        .collect();
    if output_file_names.len() == 0 {
        println!("No matches");
        return Ok(());
    }
    output_file_names.insert(0, "--- FILES TO DELETE ---".to_string());
    let mut file = NamedTempFile::new()?;
    _ = file.write_all(output_file_names.join("\n").as_bytes())?;
    let path = &file.path();
    let name = path.to_string_lossy().to_string();
    let mut output = Command::new(editor).arg(name).spawn()?;
    if output.wait().unwrap().success() {
        let _to_delete: Vec<Result<(), io::Error>> = fs::read_to_string(path)?
            .split("\n")
            .map(|p| Path::new(p))
            .filter(|p| p.exists())
            .map(|p| remove_file(p))
            .collect();
    }
    _ = file.close()?;
    Ok(())
}

#[derive(Parser, Debug)]
struct Args {
    /// Editor to use for validating files to delete
    #[arg(short = 'e', long, default_value = "vim")]
    editor: String,

    /// Clean up work directory
    #[arg(short = 'w', long)]
    clean_work: bool,

    /// Clean up this output directory, optionally with prefixes to specify what files
    /// to remove
    #[arg(short = 'c', long, default_value = "")]
    clean_output: String,

    /// Show last run name
    #[arg(short = 's', long)]
    show_last: bool,

    /// List paths from last run
    #[arg(short = 'l', long)]
    list_last: bool,

    /// Prefix of output files to remove
    #[arg(short = 'p', long, action = ArgAction::Append)]
    prefix: Vec<String>,

    /// Keep output from last run when cleaning workdir
    #[arg(short = 'k', long)]
    keep_last: bool,

    /// Specific run names to keep when cleaning workdir
    #[arg(short = 'r', long, action = ArgAction::Append)]
    keep_run: Vec<String>,
}

pub fn main() {
    let args = Args::parse();
    if args.clean_work {
        _ = clean_work(args.keep_last, args.keep_run);
    } else if args.show_last {
        let last = get_last_run();
        println!("{}", last.unwrap());
    } else if args.list_last {
        let last_files = run_files(&get_last_run().unwrap());
        for l in last_files.unwrap() {
            println!("{}", l.file_name().unwrap().to_str().unwrap());
        }
    } else if args.clean_output != "" {
        _ = clean_output(&args.clean_output, args.prefix, &args.editor);
    } else {
        println!("No flag provided");
    }
}

#[ignore = "passed"]
#[test]
fn test_get_last_run() {
    use std::env;
    let p = Path::new("/home/shannc/Bio_SDD/tutorials/nextflow-test");
    assert!(env::set_current_dir(p).is_ok());
    println!("{:?}", get_last_run());
}

#[test]
fn test_get_others() {
    use std::env;
    let p = Path::new("/home/shannc/Bio_SDD/tutorials/nextflow-test");
    assert!(env::set_current_dir(p).is_ok());
    let files = all_in_work();
    println!("{:?}", files);
}
