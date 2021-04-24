import matplotlib.pyplot as plt

plt.grid()
plt.title(r"NBLDPC EMS performance (LDPC_N96_K48_GF256_d1_exp)")
plt.xlabel(r"$E_b/N_0$")
plt.ylabel(r"$SER$")
plt.yscale("log")

SNR = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0,
       3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
BER = [9.82e-01, 9.64e-01, 9.77e-01, 8.72e-01, 8.85e-01, 8.48e-01, 7.08e-01,
       3.10e-01, 2.98e-01, 1.38e-01, 7.70e-02, 2.15e-02, 9.96e-03, 4.86e-03, 4.97e-04]
BER2 = [9.52e-01, 9.42e-01, 9.25e-01, 8.96e-01, 8.29e-01, 7.33e-01, 6.02e-01	,
        4.25e-01, 2.84e-01, 1.45e-01, 6.85e-02, 2.17e-02, 6.43e-03, 1.76e-03, 3.02e-04]
plt.plot(SNR, BER, color="blue", marker="o",
         markerfacecolor="none", label="mine")
plt.plot(SNR, BER2, color="red", marker="o",
         markerfacecolor="none", label="NBLDPC")
plt.legend()

plt.show()
