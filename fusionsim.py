import pandas as pd
import numpy as np
from sklearn.decomposition import PCA


def concatenate_and_reduce(files, n_components):
    """
    将相似性矩阵拼接并降维
    :param files: CSV 文件列表，每个文件包含一个药物-药物相似性矩阵
    :param n_components: 降维后的成分数
    :return: 降维后的相似性矩阵 DataFrame
    """
    matrices = []

    # 逐个读取文件并添加到矩阵列表中
    for file in files:
        matrix = pd.read_csv(file, index_col=0)
        matrices.append(matrix)

    # 确保所有矩阵具有相同的行和列顺序
    drug_ids = matrices[0].index
    concatenated_matrix = np.hstack([matrix.values for matrix in matrices])

    # 执行PCA降维
    pca = PCA(n_components=n_components)
    reduced_matrix = pca.fit_transform(concatenated_matrix)

    # 将降维后的矩阵转换为DataFrame
    reduced_df = pd.DataFrame(reduced_matrix, index=drug_ids)

    return reduced_df


# 文件路径列表
files = [
    'drug_similarity_disease.csv',
    'drug_similarity_side_effects.csv',
    'drug_similarity_interactions.csv'
]

# 设置降维后的成分数量
n_components = 100  # 根据需求调整降维后的维度

# 组合并降维相似性矩阵
reduced_similarity_matrix = concatenate_and_reduce(files, n_components)

# 将降维后的矩阵保存为 CSV 文件
output_file = 'reduced_drug_similarity_matrix.csv'
reduced_similarity_matrix.to_csv(output_file, index=True, header=False)

print(f"Reduced similarity matrix saved to {output_file}")
