import numpy as np
import matplotlib.pyplot as plt


def plot(x: list, y: list, label: str = ""):
    """
    Plot funtion y with given x.
    """

    plt.plot(x, y, label=label)

    plt.title("line chart")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()

    plt.savefig(label + ".png")
    plt.show()
