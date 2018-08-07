import math
import numpy
import plotter
from data_reader import Q, Q_T, Y, n, k_plus


def sigmoid(x):
    if x > 30:
        return 1
    elif x < -30:
        return 0
    return 1 / (1 + math.e ** (-x))


min_log_likelihood = 10000
min_iter = 0
min_step = 0

for STEP in {1, 2, 3, 4, 5, 6, 7, 8}:  # [0.001, 0.0005, 0.0002, 0.0001, 0.00005, 0.00001]:
    for ITERATION_NUM in [10, 100, 1000]:  # [5, 10, 100, 1000, 2000, 4000, 5000, 10000]:

        iteration_cnt = 0

        beta = [0 for _ in range(k_plus)]  # initial values

        iteration_range = []
        beta_history = [[] for _ in range(k_plus)]
        error_history = [[] for _ in range(k_plus)]

        real_beta = [3.832768685338099, -3.8888042017301028, 5.128408314113946, 1.6671674898802589]

        for _ in range(ITERATION_NUM):
            iteration_cnt += 1
            temp_vec_1 = numpy.matmul(Q, beta)
            temp_vec_2 = [Y[i] - sigmoid(temp_vec_1[i]) for i in range(n)]
            temp_vec_3 = numpy.matmul(Q_T, temp_vec_2)

            beta = [beta[i] + STEP * temp_vec_3[i] for i in range(k_plus)]

            # print("iteration # %d: %s" % (iteration_cnt, beta))

            for index in range(k_plus):
                beta_history[index].append(beta[index])
            for index in range(k_plus):
                error_history[index].append(math.fabs(beta[index] - real_beta[index]) / math.fabs(real_beta[index]))
            iteration_range.append(iteration_cnt)

        print("---------------     Step=%f, Iter=%d     ---------------------" % (STEP, iteration_cnt))

        print(beta)
        X_beta_vec = numpy.matmul(Q, beta)
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

        W = [[0 for _ in range(n)] for _ in range(n)]
        for i in range(n):
            W[i][i] = p_vec[i] * (1 - p_vec[i])
        U = [W[i][i] * X_beta_vec[i] + Y[i] - p_vec[i] for i in range(n)]

        G = numpy.matmul(
            numpy.matmul(
                Q_T
                , W
            ), Q
        )
        temp_matrix = numpy.linalg.inv(G)
        print("G: %s" % G)
        print()

        ################################# PLOTTING ######################################

        plotter.value_iteration_plot(
            iteration_range,
            beta_history,
            "Finding coefficients with gradient descent, Orthogonal\n" +
            format("step=%.2f" % STEP),
            "/results/gradient_descent_orthogonal/beta_graph/"
            + format("step_%.2f" % STEP),
            format("/iter_%d" % iteration_cnt)
        )

        plotter.value_iteration_plot(
            iteration_range,
            error_history,
            "Finding coefficients with gradient descent, Orthogonal\n error calculation,"
            + format("step=%.2f" % STEP),
            "/results/gradient_descent_orthogonal/error_graph/"
            + format("step_%.2f" % STEP),
            format("/iter_%d" % iteration_cnt)
        )


print("____________RESULTS___________")
print("min_log_likelihood=%f" % min_log_likelihood)
print("min_iter=%f" % min_iter)
print("min_step=%f" % min_step)
