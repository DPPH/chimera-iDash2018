import os
import matplotlib.pyplot as plt
import pandas as pd


def value_iteration_plot(x_values, y_values, title, figure_dir, figure_name):
    PATH = figure_dir
    if not os.path.exists(os.getcwd() + PATH):
        os.makedirs(os.getcwd() + PATH)
    NAME = format(PATH + figure_name)

    # data
    df = pd.DataFrame({'x': x_values, 'beta[0]': y_values[0], 'beta[1]': y_values[1],
                       'beta[2]': y_values[2], 'beta[3]': y_values[3]})
    num = 0
    for column in df.drop('x', axis=1):
        num += 1
        plt.plot(df['x'], df[column], marker='', linewidth=1, alpha=0.9, label=column)

    # Add legend
    plt.legend(loc=2, ncol=2)
    plt.xlabel("# of iterations")
    plt.ylabel("value")
    plt.title(title)

    # plot
    # plt.plot('x', 'y1', data=df, linestyle='-', marker='o')
    fig1 = plt.gcf()
    plt.show()
    plt.draw()
    fig1.savefig(os.getcwd() + NAME + ".png")
