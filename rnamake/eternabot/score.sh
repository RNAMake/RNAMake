rm -f strategies/weights.txt strategies/mean.txt strategies/stdev.txt
rm -f strategies/normalized/*
python score_designs.py
cp data.csv ../data.csv
