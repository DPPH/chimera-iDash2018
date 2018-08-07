import csv
import numpy

covariate_file = open('covariates.csv', 'r')
sheet = csv.reader(covariate_file)

rows = []
for row in sheet:
    rows.append(row)

rows.remove(rows[0])

# numbers
n = len(rows)
k = 3

# variables
S = []  # SNP matrix
X = []  # co-variate matrix
Y = []  # output vector

sum_up = [0, 0, 0]
cnt = [0, 0, 0]

for row in rows:
    Y += [int(row[1])]
    temp = []
    for i in [2, 3, 4]:
        if row[i].isdigit():
            temp.append(int(row[i]))
            sum_up[i - 2] += int(row[i])
            cnt[i - 2] += 1
        else:
            temp.append(-1)
    X.append(temp)

# mean filling the covariate matrix
for row in X:
    for i in range(3):
        if row[i] < 0:
            row[i] = sum_up[i] / cnt[i]

X_T = [[X[j][i] for j in range(n)] for i in range(k)]

k_plus = k + 1
X_0 = []
for row in X:
    X_0.append([1] + row)

X_0_T = [[X_0[j][i] for j in range(n)] for i in range(k_plus)]

# reading the snp matrix
snp_file = open("snpMat.txt")
snp_file.readline()
for _ in range(n):
    line = list(map(int, snp_file.readline().split()))
    S.append(line)

m = len(S[0])

Q, R = numpy.linalg.qr(X_0)
Q_T = [[Q[j][i] for j in range(n)] for i in range(k_plus)]

