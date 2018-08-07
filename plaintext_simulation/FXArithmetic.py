import numpy
from FixedPoint import FXnum


def matmul(A, B):
    if len(A[0]) != len(B):
        Exception("Length Mismatch")

    if isinstance(B[0], list):
        C = [[0 for _ in range(len(B[0]))] for _ in range(len(A))]
        for i in range(len(A)):
            for j in range(len(B[0])):
                for k in range(len(A[0])):
                    C[i][j] += A[i][k] * B[k][j]
        return C

    else:
        C = [0 for _ in range(len(A))]
        for i in range(len(A)):
            for k in range(len(A[0])):
                C[i] += A[i][k] * B[k]
        return C


def matrixInverse(A, FXFam):
    if len(A) != len(A[0]):
        Exception("Non-square Matrix!")
    n = len(A)
    B = [[float(A[i][j]) for j in range(n)] for i in range(n)]
    A_inv = numpy.linalg.inv(B)
    A_inv_fixed = [[FXnum(A_inv[i][j], FXFam) for j in range(n)] for i in range(n)]

    return A_inv_fixed


def cholesky_resolution(A, FXFam):
    k = len(A)
    AA = [[A[i][j] for i in range(k)] for j in range(k)]
    for i in range(k):
        temp = FXnum(AA[i][i], FXFam)
        alpha = FXnum(1, FXFam) / temp.sqrt()
        AA[i] = [alpha * AA[i][j] for j in range(k)]
        for j in range(k):
            AA[j][i] = AA[j][i] * alpha
        for j in range(i+1, k):
            gamma = FXnum(-AA[i][i], FXFam)
            AA[j] = [AA[j][l] + gamma * AA[i][l] for l in range(k)]
            for l in range(k):
                AA[l][j] += gamma * AA[l][i]
    return AA

