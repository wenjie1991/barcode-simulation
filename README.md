# naughty

The tiny naughty program is used to generate null hypothesis distribution of barcode frequency for lineage tricing. 
And adding the noise into the barocde.

# Usage

Exp:
```bash
naughty --seq-pool tests/test.csv \
    --out-induced tests/induced_barcode.txt \
    --out-bio tests/sample_barcode.txt \
    --out-sequen tests/seq_result.txt \
    --cell-cycle 10 \
    --pcr-cycle 15 \
    --init-cell 10:10:10
    --
```
