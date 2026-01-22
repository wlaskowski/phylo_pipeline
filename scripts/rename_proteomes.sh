#!/usr/bin/env bash

INPUT_DIR="data/proteomes/ncbi_dataset/data"
OUT_DIR="data/proteomes/renamed"

mkdir -p "$OUT_DIR"

for f in "$INPUT_DIR"/*/protein.faa; do
	genome=$(basename "$(dirname "$f")")
	out="${OUT_DIR}/${genome}.faa"

	awk -v g="${genome}|" '
		/^>/ {
			sub(/^>/, ">" g);
			print;
			next;
		}
		{ print }
	' "$f" > "$out"
done
