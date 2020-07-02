run:
	# cargo run --release tests/test.csv tests/induced_barcode.txt tests/sample_barcode.txt tests/seq_result.txt
	cargo run --release -- --seq-pool tests/test.csv --out-induced tests/induced_barcode.txt --out-bio tests/sample_barcode.txt --out-sequen tests/seq_result.txt

test:
	# cargo run --release -- \
	# --seq-pool tests/test.csv \
    # --out-barcode tests/output \
    # --cell-cycle 10 \
    # --pcr-cycle 30 \
    # --pcr-effi 0.7 \
    # --init-cell 10:10:10
	cargo run --release -- --n-cycle 20 \
    --efficiency-mean 0.5 \
    --efficiency-sd 0.2 \
    --mutation-rate 0.00001 \
    --out-prefix ./tests/pcr \
    --template ./tests/dna_template.txt

clean:
	rm -f tests/*.txt
