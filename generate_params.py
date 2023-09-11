infile = open('params.txt', 'rt')
lines = infile.readlines()

for i in range(1,11):
    outfile = open('input_files/GWAS/parameters/params_sub' + str(i) + '.txt', 'wt')
    for line in lines:
        line = line.rstrip()
        if 'training' in line:
            print('training=input_files/GWAS/subsets/training_sub' + str(i) + '.txt', file=outfile)
        elif 'testing' in line:
            print('training=input_files/GWAS/subsets/testing_sub' + str(i) + '.txt', file=outfile)
        else:
            print(line, file=outfile)

