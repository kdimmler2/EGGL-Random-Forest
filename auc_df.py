import os

outfile = open('results/auc_df_tuning5.txt', 'wt')

print('iteration', 'variant_number', 'subset1', 'subset2', 'subset3', 'subset4', 'subset5', 'subset6', 'subset7', 'subset8', 'subset9', 'subset10',
        sep='\t',
        file=outfile)

itr = 0

for root,dirs,files in os.walk('results'):
    for dir in dirs:
        if dir[0:9] == 'iteration' and 'bak' not in dir:
               itr += 1

print(itr)

for i in range(1,itr + 1):
    auc_path = 'results/iteration' + str(i) + '/training/testing_auc.list'
    infile1 = open(auc_path, 'rt')

    feature_path = 'results/iteration' + str(i) + '/training/training.table'
    with open(feature_path, 'rt') as infile2:
        features = len(infile2.readlines())

    aucs = []

    for line in infile1:
        line = line.rstrip()
        aucs.append(line)

    print('iteration' + str(i), str(features), aucs[0], aucs[1], aucs[2], aucs[3], aucs[4], aucs[5], aucs[6], aucs[7], aucs[8], aucs[9],
            sep='\t',
            file=outfile)

#for i in range(1,
#
#files_and_dirs = os.listdirs('results')
#
#for item in files_and_dirs:
#    if 


#for root,dirs,files in os.walk('results'):
#    for dir in dirs:
#        if dir[0:9] ==  'iteration' and 'bak' not in dir:
#            for file in files:
#                print(file)
#                if file == 'testing_auc.list':
#                    print(file)
