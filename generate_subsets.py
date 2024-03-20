import random
import os

infile = open('results/iteration1/training/training_RanFoG_Input.txt', 'rt')

lines = infile.readlines()

random.shuffle(lines)

tenth_of_lines = round(len(lines) * 0.1)
start_10 = 0
end_10 = tenth_of_lines
start_90 = end_10
end_90 = len(lines)

if not os.path.exists('results/iteration1/training/subsets'):
    os.mkdir('results/iteration1/training/subsets')

for i in range(1,11):
    if not os.path.exists('results/iteration1/training/subsets/subset' + str(i)):
        os.mkdir('results/iteration1/training/subsets/subset' + str(i))
    if not os.path.exists('results/iteration1/training/subsets/subset' + str(i) + '/data'):
        os.mkdir('results/iteration1/training/subsets/subset' + str(i) + '/data')

for i in range(1,11):
    if os.path.exists('results/iteration1/training/subsets/subset1/data/training.txt'):
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

        with open('results/iteration1/training/subsets/subset' + str(i) + '/data/training.txt', 'wt') as f:
            f.writelines(lines[start_90:end_90])
            f.writelines(lines[start_90_2:end_90_2])
        with open('results/iteration1/training/subsets/subset' + str(i) + '/data/testing.txt', 'wt') as f:
            f.writelines(lines[start_10:end_10])

    else:
        with open('results/iteration1/training/subsets/subset1/data/training.txt', 'wt') as f:
            f.writelines(lines[start_90:end_90])
        with open('results/iteration1/training/subsets/subset1/data/testing.txt', 'wt') as f:
            f.writelines(lines[start_10:end_10])

