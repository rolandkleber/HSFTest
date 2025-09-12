REM Run HSF P0 analysis on a single yearly RSTN file,
REM always writing to fixed filenames.

echo "[INFO] Processing $INFILE"

REM Step 1: parse → fixed file
python S10_rstn_parse.py "C:/HSFTESTS/hsf_lab/HSF000_Template/projects/HSF001_Sun/data/raw/pale_noontime-flux_2018.txt" -o ../data/10_2018parsed.csv

REM Step 2: fit spectrum → fixed file
python S20_fit_spectrum.py "C:/HSFTESTS/hsf_lab/HSF000_Template/projects/HSF001_Sun/data/raw/pale_noontime-flux_2018.txt" -o ../data/20_2018parsedspec.csv
