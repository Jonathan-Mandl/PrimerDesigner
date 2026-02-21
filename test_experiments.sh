echo "running mul_greedy_100_cds.py"
python -m Experiments.mul_greedy_100_cds

echo "running mul_greedy.py"
python -m Experiments.mul_greedy

echo "running mul_ILP.py"
python -m Experiments.mul_ILP

echo "running var_ILP_length.py"
python -m Experiments.var_ILP_length

echo "running var_ILP_num.py"
python -m Experiments.var_ILP_num

echo "running var_ILP_proteins.py"
python -m Experiments.var_ILP_proteins