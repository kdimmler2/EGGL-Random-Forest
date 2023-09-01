import random
import os

infile = open('input_files/GWAS/GWAS_RanFoG_Input.txt', 'rt')

lines = infile.readlines()

random.shuffle(lines)

tenth_of_lines = round(len(lines) * 0.1)
start_10 = 0
end_10 = tenth_of_lines
start_90 = end_10
end_90 = len(lines)

for i in range(1,11):
    if os.path.exists('input_files/GWAS/subsets/training_sub1.txt'):
        start_10 = end_10
        end_10 = start_10 + tenth_of_lines
        start_90 = 0 
        end_90 = start_10
        start_90_2 = end_10
        end_90_2 = len(lines)

        print(start_10)
        print(end_10)
        print(start_90)
        print(end_90)
        print(start_90_2)
        print(str(end_90_2) + '\n')

        with open('input_files/GWAS/subsets/training_sub' + str(i) + '.txt', 'wt') as f:
            f.writelines(lines[start_90:end_90])
            f.writelines(lines[start_90_2:end_90_2])
        with open('input_files/GWAS/subsets/testing_sub' + str(i) + '.txt', 'wt') as f:
            f.writelines(lines[start_10:end_10])

    else:
        with open('input_files/GWAS/subsets/training_sub1.txt', 'wt') as f:
            f.writelines(lines[start_90:end_90])
        with open('input_files/GWAS/subsets/testing_sub1.txt', 'wt') as f:
            f.writelines(lines[start_10:end_10])

