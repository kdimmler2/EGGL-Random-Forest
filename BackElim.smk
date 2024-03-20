import pandas as pd
import random
import os

include: 'src/utils.py'

rule all:
    input:
        'results/current_iteration/training/training.table',
        'results/current_iteration/training/training_RanFoG_Input.txt',
        expand('results/current_iteration/training/subsets/subset{num}/data/testing.txt', num = [str(i) for i in range(1,11)]),
        expand('results/current_iteration/training/subsets/subset{num}/data/training.txt', num = [str(i) for i in range(1,11)]),
        expand('results/current_iteration/training/subsets/subset{num}/params.txt', num = [str(i) for i in range(1,11)]),
        expand('results/current_iteration/training/subsets/subset{num}/TimesSelected.txt', num = [str(i) for i in range(1,11)]),
        expand('results/current_iteration/training/subsets/subset{num}/Variable_Importance.txt', num = [str(i) for i in range(1,11)]),
        expand('results/current_iteration/training/subsets/subset{num}/Trees.txt', num = [str(i) for i in range(1,11)]),
        expand('results/current_iteration/training/subsets/subset{num}/EGBV.txt', num = [str(i) for i in range(1,11)]),
        expand('results/current_iteration/training/subsets/subset{num}/Trees.test', num = [str(i) for i in range(1,11)]),
        expand('results/current_iteration/training/subsets/subset{num}/Predictions.txt', num = [str(i) for i in range(1,11)]),
        expand('results/current_iteration/training/subsets/subset{num}/auc/training_auc_values.txt', num = [str(i) for i in range(1,11)]),
        expand('results/current_iteration/training/subsets/subset{num}/auc/testing_auc_values.txt', num = [str(i) for i in range(1,11)]),
        'results/current_iteration/training/training_auc.list',
        'results/current_iteration/training/testing_auc.list',
        'results/current_iteration/training/feature_importance_scores.txt',
        'results/current_iteration/training/top_80_features.txt',
        'results/current_iteration/training/BackElim.vcf.gz',

rule gatk_table_BE:
    input:
        training_vcf = 'results/previous_iteration/training/BackElim.vcf.gz', 
    output:
        table = 'results/current_iteration/training/training.table', 
    resources:
        time    = 60,
        mem_mb  = 24000,
        cpus    = 4,
    shell:
        '''
            gatk VariantsToTable \
                -V {input.training_vcf} \
                -F CHROM -F POS -F REF -F ALT -GF GT \
                --split-multi-allelic true \
                -O {output.table}
        '''

rule ranfog_input_BE:
    input:
        table = 'results/current_iteration/training/training.table',
        phenos = config['training_pheno_file'],
    output:
        ranfog_input = 'results/current_iteration/training/training_RanFoG_Input.txt',
    resources:
        time    = 120,
        mem_mb  = 24000,
        cpus    = 4,
    run:
        df = pd.read_csv(input.table, delimiter='\t')

        df2 = df.transpose()

        for c in range(len(df2.columns)):
            for r in range(4, len(df2)):
                ref = df2.iat[2, c]
                alt = df2.iat[3, c]
                gt = df2.iat[r,c]
                if df2.iat[r,c] == ref + '/' + ref or df2.iat[r,c] == ref + '|' + ref:
                    df2.iat[r, c] = 0
                elif df2.iat[r,c] == ref + '/' + alt or df2.iat[r,c] == ref + '|' + alt:
                    df2.iat[r, c] = 1
                elif df2.iat[r,c] == alt + '/' + alt or df2.iat[r,c] == alt + '|' + alt:
                    df2.iat[r, c] = 2
                elif alt in str(df2.iat[r, c]):
                    df2.iat[r, c] = 1
                else:
                    df2.iat[r, c] = 0

        df3 = df2[4:]
        df4 = df3.reset_index()
        for r in range(len(df4)):
            split = df4.at[r, 'index'].split('.')
            df4.at[r, 'index'] = split[0]
        df4.insert(0, 'phenos', 5)

        infile = open(input.phenos, 'rt')
        phenos = {}
        sexes = {}
        for line in infile:
            line = line.rstrip()
            split = line.split(' ')
            phenos[split[1]] = split[0]
            sexes[split[1]] = split[2]

        df4.insert(2, 'sex', '')

        for r in range(len(df4)):
            if df4.at[r, 'index'] in phenos:
                df4.at[r, 'phenos'] = phenos[df4.at[r, 'index']]
            if df4.at[r, 'index'] in sexes:
                df4.at[r, 'sex'] = sexes[df4.at[r, 'index']]
        df4.to_csv(output.ranfog_input, header=False, index=False, sep=' ')

rule generate_subsets_BE:
    input:
        ranfog_input = 'results/current_iteration/training/training_RanFoG_Input.txt',
    output:
        testing_subs = expand('results/current_iteration/training/subsets/subset{num}/data/testing.txt', num = [str(i) for i in range(1,11)]),
        training_subs = expand('results/current_iteration/training/subsets/subset{num}/data/training.txt', num = [str(i) for i in range(1,11)]),
    resources:
        time    = 20,
        mem_mb  = 12000,
        cpus    = 4,
    script:
        'generate_subsets_BE.py'

rule generate_params_BE:
    input:
        all_data = expand('results/current_iteration/training/subsets/subset{num}/data/testing.txt', num = [str(i) for i in range(1,11)]),
    output:
        params_file = expand('results/current_iteration/training/subsets/subset{num}/params.txt', num = [str(i) for i in range(1,11)]),
    resources:
        time    = 20,
        mem_mb  = 12000,
        cpus    = 4,
    script:
        'generate_params_BE.py'

rule random_forest_BE:
    input:
        params_file = expand('results/current_iteration/training/subsets/subset{num}/params.txt', num = [str(i) for i in range(1,11)]),
    output:
        times_selected = 'results/current_iteration/training/subsets/subset{num}/TimesSelected.txt',
        variable_importance = 'results/current_iteration/training/subsets/subset{num}/Variable_Importance.txt',
        trees_train = 'results/current_iteration/training/subsets/subset{num}/Trees.txt',
        egbv = 'results/current_iteration/training/subsets/subset{num}/EGBV.txt',
        trees_test = 'results/current_iteration/training/subsets/subset{num}/Trees.test',
        predictions = 'results/current_iteration/training/subsets/subset{num}/Predictions.txt',
    resources:
        time    = 120,
        mem_mb  = 60000,
        cpus    = 4,
    shell:
        '''
            cp RanFoG.jar results/current_iteration/training/subsets/subset{wildcards.num}

            cd results/current_iteration/training/subsets/subset{wildcards.num}

            java -jar RanFoG.jar
        '''

rule auc_roc_BE:
    input:
        egbv = 'results/current_iteration/training/subsets/subset{num}/EGBV.txt'
    output:
        training_file = 'results/current_iteration/training/subsets/subset{num}/auc/training_auc_values.txt',
        testing_file = 'results/current_iteration/training/subsets/subset{num}/auc/testing_auc_values.txt',
    resources:
        time    = 10,
        mem_mb  = 20000,
    shell:
        '''
            cp AUC.R results/current_iteration/training/subsets/subset{wildcards.num}/auc/

            cd results/current_iteration/training/subsets/subset{wildcards.num}/auc/

            pwd

            Rscript AUC.R
        '''

rule gather_aucs_BE:
    input:
        training_aucs = expand('results/current_iteration/training/subsets/subset{num}/auc/training_auc_values.txt', num = [str(i) for i in range(1,11)]),
        testing_aucs = expand('results/current_iteration/training/subsets/subset{num}/auc/testing_auc_values.txt', num = [str(i) for i in range(1,11)]),
    output:
        training_aucs = 'results/current_iteration/training/training_auc.list',
        testing_aucs = 'results/current_iteration/training/testing_auc.list',
    resources:
        time    = 10,
        mem_mb  = 20000,
    run:
        outfile = open(output.training_aucs, 'wt')

        for file in input.training_aucs:
            infile = open(file, 'rt')
            line = infile.readline().rstrip()
            print(line, file=outfile)

        outfile = open(output.testing_aucs, 'wt')

        for file in input.testing_aucs:
            infile = open(file, 'rt')
            line = infile.readline().rstrip()
            print(line, file=outfile)

rule get_final_VI_BE:
    input:
        testing_acus = 'results/current_iteration/training/testing_auc.list',
    output:
        feature_importance = 'results/current_iteration/training/feature_importance_scores.txt',
    resources:
        time    = 10,
        mem_mb  = 20000,
    shell:
        '''
            Rscript get_bottom_features_BE_MEDIAN.R
        '''

rule top_80_features:
    input:
        table = 'results/current_iteration/training/training.table',
        feature_importance = 'results/current_iteration/training/feature_importance_scores.txt',
    output:
        top80 = 'results/current_iteration/training/top_80_features.txt',
    resources:
        time    = 10,
        mem_mb  = 20000,
    run:
        infile1 = open(input.feature_importance, 'rt')
        infile2 = open(input.table, 'rt')
        outfile = open(output.top80, 'wt')

        lines = infile1.readlines()
        lines = lines[1:]

        num_features = len(lines)

        twenty = round(0.2 * num_features)

        print(twenty)

        features = []

        top80 = lines[twenty:num_features+1]

        print(len(top80))

        for line in top80:
            line = line.rstrip()
            split = line.split(',')
            if split[0] != 1:
                features.append(int(split[0]))

        lines = infile2.readlines()

        #lines = lines[1:]

        print(len(lines))
        print(max(features))

        important_lines = []

        for feature in features:
            ind = feature-1
            important_lines.append(lines[ind])

        for line in important_lines:
            line = line.rstrip()
            split = line.split('\t')
            print(split[0] + '_' + split[1] + '_' + split[2] + '_' + split[3], file=outfile)

rule top_80_variants_BE:
    input:
        training_vcf = config['training_vcf'],
        top80 = 'results/current_iteration/training/top_80_features.txt',
    output:
        first_iteration_vcf = 'results/current_iteration/training/BackElim.vcf.gz',
    resources:
        time    = 10,
        mem_mb  = 2000,
    shell:
        '''
            bcftools view -Oz -o {output.first_iteration_vcf} --include ID=@{input.top80} {input.training_vcf}

            gatk IndexFeatureFile -I {output.first_iteration_vcf}

        '''

























