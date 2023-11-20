forest_size = 500
mtry = 100
n_features = 24861
max_branch = 2000
loss_function = 4
false_positive_cost = 1
false_negative_cost = 1

# Writing to the file using f-strings
for i in range(1,11):
    with open('input_files/GWAS/subsets/subset' + str(i) + '/params.txt', 'wt') as file:
        file.write(f"#Input files\n"
                   f"training=data/training.txt\n"
                   f"testing=data/testing.txt\n\n"
                   f"#Forest parameters\n"
                   f"ForestSize={forest_size}\n"
                   f"mtry={mtry}\n"
                   f"N_features={n_features}\n"
                   f"max_branch={max_branch}\n\n"
                   f"LossFunction={loss_function}\n"
                   f"#Personalized cost function\n"
                   f"false_positive_cost={false_positive_cost}\n"
                   f"false_negative_cost={false_negative_cost}\n")
