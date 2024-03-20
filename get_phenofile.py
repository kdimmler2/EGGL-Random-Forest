infile1 = open('/panfs/jay/groups/27/mccuem/dimml002/NuGEN/Sample_Organization/Spring2024/tuning_set10.tsv', 'rt')
infile2 = open('/panfs/jay/groups/27/mccuem/dimml002/NuGEN/Sample_Organization/Spring2024/phenotyped.tsv', 'rt')
outfile = open('tuning10.phenos', 'wt')

line = infile2.readline()

d = {}

for line in infile2:
	line = line.rstrip()
	split = line.split('\t')
	if split[0] == 'NA' and split[4] != 'unknown':
		if split[4] == 'mare':
			d[split[1]] = 2
		elif split[4] == 'gelding':
			d[split[1]] = 1
		elif split[4] == 'stallion':
			d[split[1]] = 0
		elif split[4] == 'male':
			d[split[1]] = 1
		elif split[4] == 'colt':
			d[split[1]] = 0
		elif split[4] == 'Ridgling':
			d[split[1]] = 1
		else:
			print(split[4])
	else:
		if split[0] != 'NA' and split[4] != 'unknown':
			if split[4] == 'mare':
				d[split[0]] = 2
			elif split[4] == 'gelding':
				d[split[0]] = 1
			elif split[4] == 'stallion':
				d[split[0]] = 0
			elif split[4] == 'male':
				d[split[0]] = 1
			elif split[4] == 'colt':
				d[split[0]] = 0
			elif split[4] == 'Ridgling':
				d[split[0]] = 1
			else:
				print(split[4])

print(len(d))

d2 = {}

infile2 = open('/panfs/jay/groups/27/mccuem/dimml002/NuGEN/Sample_Organization/Spring2024/phenotyped.tsv', 'rt')

line = infile2.readline()

for line in infile2:
	line = line.rstrip()
	split = line.split('\t')
	if split[0] == 'NA':
		if split[5] == 'control':
			d2[split[1]] = 0
		elif split[5] == 'case':
			d2[split[1]] = 1
		else:
			print(split[0])
	else:
		if split[5] == 'control':
			d2[split[0]] = 0
		elif split[5] == 'case':
			d2[split[0]] = 1
		else:
			print(split[0])

print(len(d2))
print(d['M6116'])

for line in infile1:
	line = line.rstrip()
	split = line.split('\t')
	if split[0] in d:
		print(str(d2[split[0]]) + ' ' + split[0] + ' ' + str(d[split[0]]), file=outfile)
