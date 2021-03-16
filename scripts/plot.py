import matplotlib.pyplot as plt

plt.grid()
plt.title(r"LDPC NMS performance (10GBPS ETHERNET, 1723/2048)")
plt.xlabel(r"$E_b/N_0$")
plt.ylabel(r"$BER$")
plt.yscale("log")

SNR = [0.0, 0.75, 1.50, 2.5, 2.75, 3.25, 3.75]
BER = [1.16e-01, 1.00e-01, 8.35e-02, 5.96e-02, 5.08e-02, 3.00e-02, 2.93e-03]
BER2 = [9.97e-02, 8.12e-02, 6.35e-02, 4.06e-02, 3.29e-02, 8.12e-03, 9.57e-05]
plt.plot(SNR, BER, color="blue", marker="o", markerfacecolor="none", label="mine")
plt.plot(SNR, BER2, color="red", marker="o", markerfacecolor="none", label="aff3ct")
plt.legend()

plt.show()
