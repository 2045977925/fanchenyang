import pandas as pd
import numpy as np

# 为Luo数据集中的药物编号匹配对应的SMILES
def drug_smiles():
    # 读取drug.txt中的DB编号到列表
    with open('mat_data/ID_drug.txt', 'r') as drug_file:
        drug_ids = [line.strip() for line in drug_file]

    # 读取DB_Drug.csv为DataFrame
    db_drug_df = pd.read_csv('DrugBank/DB_Drug.csv')

    # 筛选出DB编号对应的SMILES
    selected_data = db_drug_df[db_drug_df['DRUGBANK_ID'].isin(drug_ids)][['DRUGBANK_ID', 'SMILES']]

    # 确保输出目录存在
    output_dir = 'mat_data'
    import os
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 将结果保存到drug_smiles.csv
    selected_data.to_csv('mat_data/drug_smiles.csv', index=False, header=False)

    print("SMILES信息已保存到mat_data/drug_smiles.csv")

def matrix_hang_lie_name():
    # 读取疾病编号
    with open('mat_data/ID_disease.txt', 'r', encoding='utf-8') as file:
        disease_ids = [line.strip().split()[0] for line in file]  # 假设每行是 "编号 名称"，只读取编号

    # 读取蛋白质编号
    with open('mat_data/ID_drug.txt', 'r', encoding='utf-8') as file:
        protein_ids = [line.strip() for line in file]

    # 读取相互作用矩阵
    # 假设矩阵是用空格或逗号分隔的，需要根据实际情况调整 sep 参数
    interaction_matrix = pd.read_csv('mat_data/mat_drug_disease.txt', sep=' ', header=None)

    # 验证尺寸匹配
    if len(protein_ids) != interaction_matrix.shape[0]:
        raise ValueError("The number of protein IDs does not match the number of rows in the interaction matrix.")

    if len(disease_ids) != interaction_matrix.shape[1]:
        raise ValueError("The number of disease IDs does not match the number of columns in the interaction matrix.")

    # 为矩阵添加行和列名
    interaction_matrix.index = protein_ids
    interaction_matrix.columns = disease_ids

    # 保存带有行列名的矩阵
    interaction_matrix.to_csv('mat_data/drug_disease.csv', index=True, header=True)

    print("已保存 protein_disease.csv")

# 根据药物-药物相互作用矩阵计算药物-药物相似性
def calculate_jaccard_similarity_for_drug_interactions():
    def jaccard_similarity(vector_a, vector_b):
        intersection = np.logical_and(vector_a, vector_b).sum()
        union = np.logical_or(vector_a, vector_b).sum()
        return intersection / union if union != 0 else 0

    df = pd.read_csv('mat_data/protein_association/protein_protein.csv', index_col=0)
    drugs = df.index
    n = len(drugs)
    similarity_matrix = pd.DataFrame(np.zeros((n, n)), index=drugs, columns=drugs)

    for i, drug_a in enumerate(drugs):
        for j, drug_b in enumerate(drugs):
            if i <= j:
                sim = jaccard_similarity(df.loc[drug_a].values, df.loc[drug_b].values)
                similarity_matrix.at[drug_a, drug_b] = sim
                similarity_matrix.at[drug_b, drug_a] = sim

    similarity_matrix.to_csv('sim_network/protein_sim/protein_protein_sim.csv')

def calculate_jaccard_similarity_for_protein_interactions():
    def jaccard_similarity(vector_a, vector_b):
        intersection = np.logical_and(vector_a, vector_b).sum()
        union = np.logical_or(vector_a, vector_b).sum()
        return intersection / union if union != 0 else 0

    # 读取蛋白质-蛋白质相互作用矩阵
    df = pd.read_csv('mat_data/protein_association/protein_protein.csv', index_col=0)

    # 验证矩阵是方阵
    if df.shape[0] != df.shape[1]:
        raise ValueError("输入矩阵不是方阵，行数与列数不匹配。")

    proteins = df.index
    n = len(proteins)
    similarity_matrix = pd.DataFrame(np.zeros((n, n)), index=proteins, columns=proteins)

    for i, protein_a in enumerate(proteins):
        for j, protein_b in enumerate(proteins):
            if i <= j:  # 只计算上三角矩阵部分
                vector_a = df.loc[protein_a].values
                vector_b = df.loc[protein_b].values
                sim = jaccard_similarity(vector_a, vector_b)
                similarity_matrix.at[protein_a, protein_b] = sim
                similarity_matrix.at[protein_b, protein_a] = sim

    # 将相似性矩阵保存到文件
    similarity_matrix.to_csv('sim_network/protein_sim/protein_protein_sim.csv')

# 药物-疾病关联、药物-副作用关联、蛋白质-疾病关联
def compute_jaccard_similarity_for_drug_disease():
    # 读取CSV文件，第一列是药物编号，第一行是疾病编号
    df = pd.read_csv('mat_data/protein_association/protein_disease.csv', index_col=0)

    # 获取药物编号列表
    drug_ids = df.index.tolist()

    # 初始化相似度矩阵
    similarity_matrix = pd.DataFrame(index=drug_ids, columns=drug_ids, dtype=float)

    # 构建每个药物对应的疾病集合
    disease_sets = {drug_id: set(df.columns[df.loc[drug_id] == 1].tolist()) for drug_id in drug_ids}

    # 计算Jaccard相似度
    for i, drug_i in enumerate(drug_ids):
        for j, drug_j in enumerate(drug_ids):
            if i <= j:  # 只计算上三角矩阵，避免重复计算
                set_i = disease_sets[drug_i]
                set_j = disease_sets[drug_j]

                # 计算 Jaccard 系数
                intersection = len(set_i & set_j)
                union = len(set_i | set_j)

                # 处理空集的情况，避免除以零
                jaccard_index = intersection / union if union != 0 else 0

                # 更新相似度矩阵
                similarity_matrix.loc[drug_i, drug_j] = jaccard_index
                similarity_matrix.loc[drug_j, drug_i] = jaccard_index  # 利用对称性

    # 将相似度矩阵保存为 CSV 文件，包含行列名
    similarity_matrix.to_csv('sim_network/protein_sim/protein_disease_sim.csv', index=True, header=True)


if __name__ == '__main__':
    # os.remove('mat_data/mat_drug_protein_remove_homo.txt')
    # drug_smiles()
    # matrix_hang_lie_name()
    """
    input_output_pairs = [
        ('mat_data/drug_association/drug_disease.csv', 'sim_network/drug_sim/drug_disease_sim.csv'),
        ('mat_data/drug_association/drug_drug.csv', 'sim_network/drug_sim/drug_drug_sim.csv'),
        ('mat_data/drug_association/drug_protein.csv', 'sim_network/drug_sim/drug_protein_sim.csv'),
        ('mat_data/drug_association/drug_se.csv', 'sim_network/drug_sim/drug_se_sim.csv')
    ]
    """
    # calculate_jaccard_similarity_for_drug_interactions()
    # calculate_jaccard_similarity_for_protein_interactions()
    compute_jaccard_similarity_for_drug_disease()

