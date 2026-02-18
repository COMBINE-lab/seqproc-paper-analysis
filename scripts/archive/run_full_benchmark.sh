#!/bin/bash
# Benchmark script for Full Dataset (SRR13948564) - Fwd and RC passes

DATA_DIR="/home/ubuntu/combine-lab/seqproc-paper-analysis-clean/data"
CONFIG_DIR="/home/ubuntu/combine-lab/seqproc-paper-analysis-clean/configs"
SEQPROC_BIN="/home/ubuntu/combine-lab/seqproc/target/release/seqproc"
MATCHBOX_BIN="/home/ubuntu/combine-lab/matchbox/target/release/matchbox"

FULL_FQ="${DATA_DIR}/SRR13948564_full.fastq"

echo "Benchmarking Seqproc Forward..."
/usr/bin/time -v $SEQPROC_BIN --geom ${CONFIG_DIR}/seqproc/splitseq_singleend_primer.geom --file1 $FULL_FQ --out1 /dev/null --threads 8 2> seqproc_fwd_bench.log

echo "Benchmarking Seqproc RC..."
/usr/bin/time -v $SEQPROC_BIN --geom ${CONFIG_DIR}/seqproc/splitseq_singleend_rc.geom --file1 $FULL_FQ --out1 /dev/null --threads 8 2> seqproc_rc_bench.log

echo "Benchmarking Matchbox Forward..."
/usr/bin/time -v $MATCHBOX_BIN -e 0.2 -t 8 -s ${CONFIG_DIR}/matchbox/splitseq_singleend.mb $FULL_FQ > /dev/null 2> matchbox_fwd_bench.log

echo "Benchmarking Matchbox RC..."
/usr/bin/time -v $MATCHBOX_BIN -e 0.2 -t 8 -s ${CONFIG_DIR}/matchbox/splitseq_singleend_rc.mb $FULL_FQ > /dev/null 2> matchbox_rc_bench.log

echo "Done."
