import math

import numpy

import plotter
from data_reader import X_0, X_0_T, Y, n, k_plus


def sigmoid(x):
    if x > 30:
        return 1
    elif x < -30:
        return 0
    return 1 / (1 + math.e ** (-x))


min_log_likelihood = 10000
min_iter = 0
min_step = 0

for STEP in [0.0000001, 0.00000005, 0.00000001]:
    for ITERATION_NUM in [100, 1000, 2000, 4000, 5000, 10000]:

        iteration_cnt = 0
        beta = [0 for _ in range(k_plus)]  # initial values

        iteration_range = []
        beta_history = [[] for _ in range(k_plus)]
        error_history = [[] for _ in range(k_plus)]

        real_beta = [-0.0002439672516133533, 0.023700806140536313, 0.032954116581954275, -0.02316932312584364]

        for _ in range(ITERATION_NUM):
            iteration_cnt += 1
            temp_vec_1 = numpy.matmul(X_0, beta)
            temp_vec_2 = [Y[i] - sigmoid(temp_vec_1[i]) for i in range(n)]
            temp_vec_3 = numpy.matmul(X_0_T, temp_vec_2)

            beta = [beta[i] + STEP * temp_vec_3[i] for i in range(k_plus)]

            # print("iteration # %d: %s" % (iteration_cnt, beta))

            for index in range(k_plus):
                beta_history[index].append(beta[index])
            for index in range(k_plus):
                error_history[index].append(math.fabs(beta[index] - real_beta[index]) / math.fabs(real_beta[index]))
            iteration_range.append(iteration_cnt)

        print("---------------   Step=%f, Iter=%d   ---------------------" % (STEP, iteration_cnt))

        print(beta)
        X_beta_vec = numpy.matmul(X_0, beta)
        p_vec = [sigmoid(X_beta_vec[i]) for i in range(n)]
        print(p_vec)

        error = 0
        for i in range(n):
            error += (Y[i] - p_vec[i]) ** 2
        error /= n

        print("Coefficient Error = %f" % error)

        likelihood = 0
        for i in range(n):
            likelihood += math.log(1 + math.e ** X_beta_vec[i], math.e) - Y[i] * X_beta_vec[i]
        print("-log(Likelihood): %f" % likelihood)

        if likelihood < min_log_likelihood:
            min_log_likelihood = likelihood
            min_step = min_step
            min_iter = iteration_cnt

        # W = [[0 for _ in range(n)] for _ in range(n)]
        # for i in range(n):
        #     W[i][i] = p_vec[i] * (1 - p_vec[i])
        # U = [W[i][i] * X_beta_vec[i] + Y[i] - p_vec[i] for i in range(n)]
        #
        # AAA = numpy.matmul(
        #     numpy.matmul(
        #         X_0_T
        #         , W
        #     ), X_0
        # )
        # temp_matrix = numpy.linalg.inv(AAA)

        print()

        # ################################# PLOTTING ######################################

        plotter.value_iteration_plot(
            iteration_range,
            beta_history,
            "Finding coefficients with gradient descent\n" +
            format("step=%f" % STEP),
            "/results/gradient_descent/beta_graph/"
            + format("step_%f" % STEP),
            format("/iter_%d" % iteration_cnt)
        )

        plotter.value_iteration_plot(
            iteration_range,
            error_history,
            "Finding coefficients with gradient descent\n" +
            format("error calculation, step=%.8f" % STEP),
            "/results/gradient_descent/error_graph/"
            + format("step_%.8f" % STEP),
            format("/iter_%d" % iteration_cnt)
        )

print("____________RESULTS___________")
print("min_log_likelihood=%f" % min_log_likelihood)
