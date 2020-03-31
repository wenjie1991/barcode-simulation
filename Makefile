run:
	mkdir -p tests
	cargo run --release tests/test.csv tests/induced_barcode.txt tests/sample_barcode.txt tests/seq_result.txt
