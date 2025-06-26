# 安装必要库（如未安装）
!pip install mpmath scipy tqdm

import numpy as np
import matplotlib.pyplot as plt
from mpmath import zeta, arg, mpc
from scipy.stats import entropy, skew, kurtosis
from tqdm import tqdm

# 样本参数
t_values = np.linspace(10, 50, 2000)
sigma_list = [0.4, 0.5, 0.6]

# 存储结果
results = {}

for sigma in sigma_list:
    print(f"σ={sigma}:")
    
    # 计算相位角 θ(t)
    theta_vals = []
    for t in tqdm(t_values):
        s = mpc(sigma, t)
        val = arg(zeta(s))
        theta_vals.append(float(val))

    # 导数计算
    theta_vals = np.unwrap(theta_vals)
    theta_prime = np.gradient(theta_vals, t_values)
    theta_double_prime = np.gradient(theta_prime, t_values)
    theta_triple_prime = np.gradient(theta_double_prime, t_values)

    # 三阶导数标准差与信息熵
    std_dev = np.std(theta_triple_prime)
    hist, _ = np.histogram(theta_triple_prime, bins=100, density=True)
    hist = hist[hist > 0]
    shannon_entropy = entropy(hist, base=2)

    # 偏度与峰度
    skewness = skew(theta_triple_prime)
    kurt = kurtosis(theta_triple_prime)

    # 存入结果
    results[sigma] = {
        "三阶导数标准差 (Stability)": std_dev,
        "导数分布熵 (Shannon Entropy)": shannon_entropy,
        "偏度 (Skewness)": skewness,
        "峭度 (Kurtosis)": kurt
    }

# 转为 DataFrame
import pandas as pd
df = pd.DataFrame(results).T
df.index.name = "σ"
df
