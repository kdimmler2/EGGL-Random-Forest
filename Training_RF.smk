import pandas as pd

include: 'src/utils.py'

rule all:
    input:
        'input_files/training/generation/training.table',
        'input_files/training/generation/training_RanFogInput.txt',

rule gatk_table:
    input:
        training_vcf = config['training_vcf'], 
    output:
        table = 'input_files/training/generation/training.table', 
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

rule ranfog_input:
    input:
        table = 'input_files/training/generation/training.table',
    output:
        ranfog_input = 'input_files/training/generation/training_RanFogInput.txt',
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

        infile = open(config['pheno_file'], 'rt')
        phenos = {}
        for line in infile:
            line = line.rstrip()
            split = line.split(' ')
            phenos[split[1]] = split[0]

        for r in range(len(df4)):
            if df4.at[r, 'index'] in phenos:
                df4.at[r, 'phenos'] = phenos[df4.at[r, 'index']]

        df4.to_csv(output.ranfog_input, header=False, index=False, sep=' ')

