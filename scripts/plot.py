import matplotlib.pyplot as plt

plt.grid()
plt.title(r"LDPC NMS performance (10GBPS ETHERNET, 1723/2048)")
plt.xlabel(r"$E_b/N_0$")
plt.ylabel(r"$BER$")
plt.yscale("log")

SNR = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
BER = [9.75e-02, 8.43e-02, 7.17e-02, 6.00e-02, 4.99e-02, 4.18e-02, 3.36e-02]
BER2 = [1.01e-01, 8.81e-02, 7.63e-02, 6.42e-02, 5.24e-02, 4.15e-02, 2.67e-02]
plt.plot(SNR, BER, color="blue", marker="o", markerfacecolor="none", label="mine")
plt.plot(SNR, BER2, color="red", marker="o", markerfacecolor="none", label="aff3ct")
plt.legend()

plt.show()
