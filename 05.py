# 安装依赖（如未安装）
!pip install mpmath scipy tqdm --quiet

# 导入必要模块
import numpy as np
import matplotlib.pyplot as plt
from mpmath import zeta, arg, mpc
from scipy.stats import entropy
from scipy.signal import savgol_filter
from tqdm import tqdm
import pandas as pd
from IPython.display import display  # ✅ 修正处

# 设置 σ 范围为 0.01 到 0.99（步长 0.01），共 99 个点
sigma_scan = np.arange(0.01, 1.00, 0.01)
t_values = np.linspace(10, 50, 1000)  # t 范围

# 初始化熵结果
entropy_values = []

# 主循环：对每个 σ 计算相位三阶导数的熵
for sigma in tqdm(sigma_scan, desc="扫描 σ ∈ [0.01, 0.99]"):
    theta_vals = []
    for t in t_values:
        s = mpc(sigma, t)
        val = arg(zeta(s))
        theta_vals.append(float(val))

    # 解包相位角，平滑去噪
    theta_vals = np.unwrap(theta_vals)
    theta_vals = savgol_filter(theta_vals, 51, 3)

    # 计算三阶导数
    theta_prime = np.gradient(theta_vals, t_values)
    theta_double_prime = np.gradient(theta_prime, t_values)
    theta_triple_prime = np.gradient(theta_double_prime, t_values)

    # 剔除边缘效应
    triple_clean = theta_triple_prime[100:-100]

    # 计算 Shannon 熵
    hist, _ = np.histogram(triple_clean, bins=100, density=True)
    hist = hist[hist > 0]
    shannon_entropy = entropy(hist, base=2)
    entropy_values.append(shannon_entropy)

# 构建 DataFrame
df = pd.DataFrame({
    "σ": sigma_scan,
    "Shannon Entropy": entropy_values
})
display(df)  # ✅ 显示表格

# 绘制熵 vs σ 图
plt.figure(figsize=(10, 5))
plt.plot(df["σ"], df["Shannon Entropy"], marker='o', color='darkblue')
plt.title("Shannon Entropy of ζ Phase Derivative vs σ ∈ [0.01, 0.99]", fontsize=14)
plt.xlabel("σ", fontsize=12)
plt.ylabel("Shannon Entropy", fontsize=12)
plt.grid(True)
plt.show()
