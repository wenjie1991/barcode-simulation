# pcr-simulation

The tiny naughty program is used to generate null hypothesis distribution of barcode frequency for lineage tricing. 
And adding the noise into the barocde.

# Usage

```bash
cargo run --release -- --n-cycle 25 \
    --efficiency-mean 0.5 \
    --efficiency-sd 0.2 \
    --mutation-rate 0.00001 \
    --out-prefix ./tests/pcr \
    --template ./tests/dna_template.txt
```
