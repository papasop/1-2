import mpmath
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import IPython.display as display

mpmath.mp.dps = 100  # 提高精度
t_values = np.arange(1000, 100001, 5000)  # 增大间隔减少噪声
sigma_values = [0.5, 0.6]
theta_deriv_3 = {sigma: [] for sigma in sigma_values}

# 计算 theta(t) 和三阶导数
for sigma in sigma_values:
    theta_values = []
    for t in t_values:
        t_mp = float(t)
        s = mpmath.mpc(sigma, t_mp)
        # Hardy Z 调整
        gamma_arg = mpmath.gamma(0.25 + t_mp * mpmath.mpc(0, 0.5)) / mpmath.sqrt(mpmath.pi)
        chi_t = mpmath.power(mpmath.pi, -t_mp * mpmath.mpc(0, 0.5)) * gamma_arg
        theta_adjusted = mpmath.arg(chi_t)
        zeta_val = mpmath.zeta(s)
        z_t = zeta_val * mpmath.exp(-mpmath.mpc(0, 1) * theta_adjusted)
        if abs(mpmath.im(z_t)) < 1e-10:  # 确保 Z(t) 接近实值
            theta_t = 0
        else:
            theta_t = mpmath.arg(z_t)
        theta_deg = float(theta_t * 180 / mpmath.pi)
        theta_values.append(theta_deg)

    # 解缠绕相位
    theta_unwrapped = np.unwrap(np.array(theta_values) * np.pi / 180) * 180 / np.pi

    # 数值微分（五点公式提高精度）
    h = 5000  # 匹配 t 间隔
    for i in range(len(t_values)):
        if i > 2 and i < len(t_values) - 3:
            theta_i = theta_unwrapped[i]
            theta_ip2 = theta_unwrapped[i + 2]
            theta_ip1 = theta_unwrapped[i + 1]
            theta_im1 = theta_unwrapped[i - 1]
            theta_im2 = theta_unwrapped[i - 2]
            theta_triple = (theta_ip2 - 2 * theta_ip1 + 2 * theta_im1 - theta_im2) / (2 * h ** 3)
        else:
            theta_triple = 0
        theta_deriv_3[sigma].append(theta_triple)

# 表格
df = pd.DataFrame({
    "t": t_values,
    "θ'''(t) σ=0.5 [deg/unit³]": np.round(theta_deriv_3[0.5], 5),
    "θ'''(t) σ=0.6 [deg/unit³]": np.round(theta_deriv_3[0.6], 5)
})
display.display(df)

# 绘制
plt.plot(t_values, theta_deriv_3[0.5], 'b-', label="σ = 0.5")
plt.plot(t_values, theta_deriv_3[0.6], 'r-', label="σ = 0.6")
plt.title("Third Derivative θ'''(t) vs t with Hardy Z Adjustment")
plt.xlabel("t")
plt.ylabel("θ'''(t) [deg/unit³]")
plt.legend()
plt.grid()
plt.show()

# 统计分析
std_0_5 = np.std(theta_deriv_3[0.5])
std_0_6 = np.std(theta_deriv_3[0.6])
print(f"Standard deviation σ=0.5: {std_0_5:.5e}")
print(f"Standard deviation σ=0.6: {std_0_6:.5e}")

