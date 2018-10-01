import os

#Get current working directory the files are in
cwd = os.getcwd()
#Get actual folder names
x = next(os.walk(cwd))[1]

file = cwd + '/' + x[0] + '/output.txt'

with open(file) as f:
    bmean = []
    bbest = []
    gens = []

    line = f.readline()
    while line:
        row = line.strip().split()
        i = 0

        while i < len(row):
            if row[i] == "Gen:":
                gens.append(row[i+1])

            if row[i] == "Best_Mean:":
                bmean.append(row[i+1])
        
            sif row[i] == "Best_Best:":
                bmean.append(row[i+1])


print()        