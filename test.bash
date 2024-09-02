echo "RUNNING TEST MAIN.PY ... "
rm -rf ../test_dir/output/
python3 commec/cli.py screen -d ../../common-mechanism-dbs/ -o ../test_dir/output/ ../test_dir/igem_test_queries.fasta -t 1