@echo off
python 10_linewidth_TP_sweep.py ^
  --hitran ..\data\raw\H2O_lines_7150-7220.csv ^
  --auto-strongest ^
  --T 250,300,350,400,500,700 ^
  --P 0.2,0.5,1.0 ^
  --instr-fwhm 0.10 ^
  --win 1.0 ^
  --out ..\results\MOL-01_H2O_TP_Voigt
