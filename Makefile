run:
	# cargo run --release tests/test.csv tests/induced_barcode.txt tests/sample_barcode.txt tests/seq_result.txt
	cargo run --release -- --seq-pool tests/test.csv --out-induced tests/induced_barcode.txt --out-bio tests/sample_barcode.txt --out-sequen tests/seq_result.txt

test:
	cargo run --release -- --seq-pool tests/test.csv \
    --out-induced tests/induced_barcode.txt \
    --out-bio tests/sample_barcode.txt \
    --out-sequen tests/seq_result.txt \
    --cell-cycle 10 \
    --pcr-cycle 15 \
    --init-cell 10:10:10

clean:
	rm -f tests/*.txt
