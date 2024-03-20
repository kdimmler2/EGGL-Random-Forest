import math

with open('results/current_iteration/training/training.table', 'rt') as infile:
    features = len(infile.readlines())

mtry_param = math.ceil(math.sqrt(features))

forest_size = 5000
mtry = mtry_param
n_features = features
max_branch = 6000
loss_function = 1

# Writing to the file using f-strings
for i in range(1,11):
    with open('results/current_iteration/training/subsets/subset' + str(i) + '/params.txt', 'wt') as file:
        file.write(f"#Input files\n"
                   f"training=data/training.txt\n"
                   f"testing=data/testing.txt\n\n"
                   f"#Forest parameters\n"
                   f"ForestSize={forest_size}\n"
                   f"mtry={mtry}\n"
                   f"N_features={n_features}\n"
                   f"max_branch={max_branch}\n\n"
                   f"LossFunction={loss_function}")
