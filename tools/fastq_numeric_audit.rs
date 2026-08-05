//! Exact, bounded-memory FASTQ validator for numeric ENA accession IDs.
//!
//! Build with:
//!   rustc -C opt-level=3 -C target-cpu=x86-64-v3 -o tools/bin/fastq-numeric-audit tools/fastq_numeric_audit.rs

use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

const BUFFER_BYTES: usize = 8 * 1024 * 1024;
const DOMAIN: &[u8] = b"fastq_numeric_accession_set_v1\0";

fn fail(message: impl AsRef<str>) -> ! {
    eprintln!("fastq-numeric-audit: {}", message.as_ref());
    std::process::exit(2);
}

fn read_line(reader: &mut impl BufRead, buffer: &mut Vec<u8>) -> io::Result<usize> {
    buffer.clear();
    reader.read_until(b'\n', buffer)
}

fn sequence_length(line: &[u8]) -> usize {
    let without_lf = line.strip_suffix(b"\n").unwrap_or(line);
    without_lf.strip_suffix(b"\r").unwrap_or(without_lf).len()
}

fn header_id(header: &[u8]) -> &[u8] {
    let body = &header[1..];
    let end = body
        .iter()
        .position(|byte| byte.is_ascii_whitespace())
        .unwrap_or(body.len());
    let value = &body[..end];
    if value.ends_with(b"/1") || value.ends_with(b"/2") {
        &value[..value.len() - 2]
    } else {
        value
    }
}

fn parse_numeric_id(value: &[u8]) -> Option<(&[u8], usize)> {
    let separator = value.iter().rposition(|byte| *byte == b'.')?;
    let prefix = &value[..separator];
    let digits = &value[separator + 1..];
    if prefix.is_empty() || digits.is_empty() {
        return None;
    }
    let mut number = 0usize;
    for byte in digits {
        if !byte.is_ascii_digit() {
            return None;
        }
        number = number
            .checked_mul(10)?
            .checked_add((byte - b'0') as usize)?;
    }
    Some((prefix, number))
}

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        fail("usage: fastq-numeric-audit FASTQ MATE NUMERIC_ID_MAX CANONICAL_OUTPUT");
    }
    let input = Path::new(&args[1]);
    let mate: usize = args[2]
        .parse()
        .unwrap_or_else(|_| fail("MATE must be a non-negative integer"));
    let numeric_id_max: usize = args[3]
        .parse()
        .ok()
        .filter(|value| *value > 0)
        .unwrap_or_else(|| fail("NUMERIC_ID_MAX must be a positive integer"));
    let canonical_output = Path::new(&args[4]);

    let input_handle = File::open(input)?;
    let mut reader = BufReader::with_capacity(BUFFER_BYTES, input_handle);
    let mut header = Vec::with_capacity(512);
    let mut sequence = Vec::with_capacity(512);
    let mut separator = Vec::with_capacity(512);
    let mut quality = Vec::with_capacity(512);
    let mut bitmap = vec![0u8; (numeric_id_max + 7) / 8];
    let mut accession_prefix: Option<Vec<u8>> = None;
    let mut records = 0usize;

    loop {
        if read_line(&mut reader, &mut header)? == 0 {
            break;
        }
        records += 1;
        if read_line(&mut reader, &mut sequence)? == 0
            || read_line(&mut reader, &mut separator)? == 0
            || read_line(&mut reader, &mut quality)? == 0
        {
            fail(format!(
                "truncated FASTQ record {records} in {}",
                input.display()
            ));
        }
        if !header.starts_with(b"@") || !separator.starts_with(b"+") {
            fail(format!(
                "malformed FASTQ record {records} in {}",
                input.display()
            ));
        }
        if sequence_length(&sequence) != sequence_length(&quality) {
            fail(format!(
                "sequence/quality length mismatch in FASTQ record {records} in {}",
                input.display()
            ));
        }

        let id = header_id(&header);
        let (prefix, numeric_id) = parse_numeric_id(id).unwrap_or_else(|| {
            fail(format!(
                "FASTQ record {records} does not have a numeric accession ID"
            ))
        });
        match &accession_prefix {
            None => accession_prefix = Some(prefix.to_vec()),
            Some(expected) if expected.as_slice() == prefix => {}
            Some(expected) => fail(format!(
                "FASTQ record {records} changes accession prefix from {:?} to {:?}",
                String::from_utf8_lossy(expected),
                String::from_utf8_lossy(prefix)
            )),
        }
        if numeric_id == 0 || numeric_id > numeric_id_max {
            fail(format!(
                "FASTQ record {records} numeric ID {numeric_id} is outside 1..{numeric_id_max}"
            ));
        }
        let bit = numeric_id - 1;
        let mask = 1u8 << (bit % 8);
        if bitmap[bit / 8] & mask != 0 {
            fail(format!(
                "duplicate numeric accession ID {} in {}",
                String::from_utf8_lossy(id),
                input.display()
            ));
        }
        bitmap[bit / 8] |= mask;
    }

    let output_handle = File::create(canonical_output)?;
    let mut output = BufWriter::with_capacity(BUFFER_BYTES, output_handle);
    output.write_all(DOMAIN)?;
    write!(output, "{mate}\0{numeric_id_max}\0")?;
    if let Some(prefix) = accession_prefix {
        output.write_all(&prefix)?;
    }
    output.write_all(b"\0")?;
    output.write_all(&bitmap)?;
    output.flush()?;
    println!("{{\"records\":{records}}}");
    Ok(())
}
