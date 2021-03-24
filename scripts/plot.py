import matplotlib.pyplot as plt

plt.grid()
plt.title(r"LDPC NMS performance (10GBPS ETHERNET, 1723/2048)")
plt.xlabel(r"$E_b/N_0$")
plt.ylabel(r"$BER$")
plt.yscale("log")

SNR = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
BER = [1.15e-01, 1.07e-01, 9.58e-02, 8.35e-02, 7.04e-02, 6.02e-02, 4.18e-02, 1.05e-02, 2.50e-04]
BER2 = [1.04e-01, 9.05e-02, 7.90e-02, 6.67e-02, 5.53e-02, 4.37e-02, 2.79e-02, 4.49e-03, 3.31e-05]
plt.plot(SNR, BER, color="blue", marker="o", markerfacecolor="none", label="mine")
plt.plot(SNR, BER2, color="red", marker="o", markerfacecolor="none", label="aff3ct")
plt.legend()

plt.show()
