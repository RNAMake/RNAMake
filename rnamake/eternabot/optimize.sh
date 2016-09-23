rm -f strategies/weights.txt strategies/mean.txt strategies/stdev.txt
rm -f strategies/normalized/*
python ensemble_design.py sparse optimization_run "((((....))))" NNNNNNNNNNNN 1
