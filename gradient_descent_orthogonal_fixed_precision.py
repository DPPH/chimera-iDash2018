import math
from FixedPoint import FXnum, FXfamily
import FXArithmetic

from data_reader import Q, Q_T, Y, n, k_plus
import plotter

PRECISION_NUM = 12
PRECISION = FXfamily(PRECISION_NUM)


####################################################################


def sigmoid(x):
    if x > 30:
        return FXnum(1, PRECISION)
    elif x < -30:
        return FXnum(0, PRECISION)
    return FXnum(1 / (1 + math.e ** (-x)), PRECISION)


Q_fx = [[FXnum(Q[i][j], PRECISION) for j in range(len(Q[0]))] for i in range(len(Q))]
Q_T_fx = [[FXnum(Q_T[i][j], PRECISION) for j in range(len(Q_T[0]))] for i in range(len(Q_T))]

Y_fx = [FXnum(Y[i], PRECISION) for i in range(len(Y))]
min_log_likelihood = 10000
min_iter = 0
min_step = 0

range_min = 1000000
range_max = -1000000

for i in range(n):
    for j in range(k_plus):
        if Q[i][j] > range_max:
            range_max = Q[i][j]
        if Q[i][j] < range_min:
            range_min = Q[i][j]

for STEP in [0.5, 1, 2, 4, 8]:
    for ITERATION_NUM in [5, 10, 50, 100]:

        iteration_cnt = 0

        beta = [FXnum(0, PRECISION) for _ in range(k_plus)]  # initial values
        iteration_range = []
        beta_history = [[] for _ in range(k_plus)]
        error_history = [[] for _ in range(k_plus)]
        real_beta = [3.832768685338099, -3.8888042017301028, 5.128408314113946, 1.6671674898802589]

        for _ in range(ITERATION_NUM):
            iteration_cnt += 1
            temp_vec_1 = FXArithmetic.matmul(Q_fx, beta)
            temp_vec_2 = [FXnum(Y_fx[i] - sigmoid(temp_vec_1[i]), PRECISION) for i in range(n)]
            temp_vec_3 = FXArithmetic.matmul(Q_T_fx, temp_vec_2)

            beta = [FXnum(beta[i] + STEP * temp_vec_3[i], PRECISION) for i in range(k_plus)]

            for i in range(k_plus):
                if float(beta[i]) > range_max:
                    range_max = float(beta[i])
                if float(beta[i]) < range_min:
                    range_min = float(beta[i])

            # print("iteration # %d: %s" % (iteration_cnt, beta))

            for index in range(k_plus):
                beta_history[index].append(float(beta[index]))
            for index in range(k_plus):
                error_history[index].append(
                    math.fabs(float(beta[index]) - real_beta[index]) / math.fabs(real_beta[index]))
            iteration_range.append(iteration_cnt)

        print("------------------   Step=%.2f Iter=%d   ---------------------" % (STEP, iteration_cnt))

        print("beta: %s" % [float(beta[i]) for i in range(len(beta))])
        X_beta_vec = FXArithmetic.matmul(Q_fx, beta)
        p_vec = [sigmoid(X_beta_vec[i]) for i in range(n)]
        # print("probability vector: %s" % p_vec)

        # error = 0
        # for i in range(n):
        #     error += (Y[i] - float(p_vec[i])) ** 2
        # error /= n

        # print("Coefficient Error : %f" % error)

        # likelihood = 0
        # for i in range(n):
        #     likelihood += math.log(1 + math.e ** X_beta_vec[i], math.e) - Y[i] * X_beta_vec[i]
        # print("-log(Likelihood): %f" % likelihood)

        # if likelihood < min_log_likelihood:
        #     min_log_likelihood = likelihood
        #     min_step = min_step
        #     min_iter = iteration_cnt

        W = [[FXnum(0, PRECISION) for _ in range(n)] for _ in range(n)]
        for i in range(n):
            W[i][i] = p_vec[i] * (FXnum(1, PRECISION) - p_vec[i])
        # U = [W[i][i] * X_beta_vec[i] + Y[i] - p_vec[i] for i in range(n)]

        G = FXArithmetic.matmul(
            FXArithmetic.matmul(
                Q_T_fx
                , W
            ), Q_fx
        )

        for i in range(k_plus):
            for j in range(k_plus):
                if G[i][j] > range_max:
                    range_max = G[i][j]
                if G[i][j] < range_min:
                    range_min = G[i][j]

        ch = FXArithmetic.cholesky_resolution(G, PRECISION)

        print([list(map(float, ch[i])) for i in range(k_plus)])

        for i in range(k_plus):
            for j in range(k_plus):
                if ch[i][j] > range_max:
                    range_max = ch[i][j]
                if ch[i][j] < range_min:
                    range_min = ch[i][j]

        print("Ranges: [%f, %f]" % (range_min, range_max))
        print()

        # try:
        #     G_inv = FXArithmetic.matrixInverse(G, PRECISION)
        #     # print("G: %s" % [[float(G[i][j]) for j in range(k_plus)] for i in range(k_plus)])
        #     # print("G_Inverse: %s" % [[float(G_inv[i][j]) for j in range(k_plus)] for i in range(k_plus)])
        #     print()
        #     ID = FXArithmetic.matmul(G, G_inv)
        #     # print("ID: %s" % [[float(ID[i][j]) for j in range(k_plus)] for i in range(k_plus)])
        #     print()
        # except:
        #     print("Singular Matrix")

        ######################################## PLOTTING ######################################

        plotter.value_iteration_plot(
            iteration_range,
            beta_history,
            "Finding coefficients with gradient descent\n" +
            format("step=%.2f, precision=%d bits" % (STEP, PRECISION_NUM)),
            "/results/gradient_descent_orthogonal_fixed_precision/"
            + format("precision_%d/beta_graph/step_%.2f" % (PRECISION_NUM, STEP)),
            format("/iter_%d" % iteration_cnt)
        )

        plotter.value_iteration_plot(
            iteration_range,
            error_history,
            "Finding coefficients with gradient descent\n error calculation,"
            + format("step=%.2f, precision=%d bits" % (STEP, PRECISION_NUM)),
            "/results/gradient_descent_orthogonal_fixed_precision/"
            + format("precision_%d/error_graph/step_%.2f" % (PRECISION_NUM, STEP)),
            format("/iter_%d" % iteration_cnt)
        )

# print("____________RESULTS___________")
# print("min_log_likelihood=%f" % min_log_likelihood)
# print("min_iter=%f" % min_iter)
# print("min_step=%f" % min_step)
