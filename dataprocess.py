import pandas as pd
import os
from rdkit import Chem

def drug_in_hetero():
    cid_db_df = pd.read_csv('davis/CIDDB.txt', sep='\s+', header=None, names=['CID', 'DB'])
    drug_df = pd.read_csv('mat_data/ID_drug.txt', header=None, names=['DB'])
    unique_db_in_cid_db = cid_db_df['DB'].dropna().unique()
    all_exist = all(db in drug_df['DB'].values for db in unique_db_in_cid_db)
    if all_exist:
        print('CIDDB.txt中所有的DB编号都存在于drug.txt中')
    else:
        print('CIDDB.txt中有一些DB编号不在drug.txt中')
    missing_db = [db for db in unique_db_in_cid_db if db not in drug_df['DB'].values]
    if missing_db:
        print('以下DB编号在drug.txt中缺失：')
        print(missing_db)

# 读取TXT文件中的药物编号
def read_txt_file(file_path):
    with open(file_path, 'r') as file:
        drug_ids = {line.strip() for line in file}  # 使用集合来存储药物编号，去除换行符
    return drug_ids

# 读取CSV文件中的药物编号
def read_csv_file(file_path):
    df = pd.read_csv(file_path)
    if 'DrugID' in df.columns:
        drug_ids = set(df['DrugID'])
    else:
        raise ValueError("CSV文件中没有找到名为'DrugID'的列")
    return drug_ids

# 判断CSV文件中的药物ID是否都在TXT文件中
def compare_drug_ids(txt_file_path, csv_file_path):
    txt_drug_ids = read_txt_file(txt_file_path)
    csv_drug_ids = read_csv_file(csv_file_path)

    # 检查CSV中的药物ID是否都在TXT文件中
    missing_ids = csv_drug_ids - txt_drug_ids
    if not missing_ids:
        print("CSV文件中的所有药物ID都存在于TXT文件中。")
    else:
        print("以下药物ID在TXT文件中不存在：")
        for drug_id in missing_ids:
            print(drug_id)

# 为SMILES匹配:
def davis_SMILES_in_DB():
    # 读取CSV文件
    df1 = pd.read_csv('davis/davis_drug_smiles.csv')
    df2 = pd.read_csv('DrugBank/DB_Drug.csv')

    # 提取SMILES列
    smiles_1 = set(df1['Drug SMILES'])
    smiles_2 = set(df2['SMILES'])

    # 查找在df1中但不在df2中的SMILES
    missing_smiles = smiles_1 - smiles_2

    # 输出结果
    if missing_smiles:
        print("The following SMILES are in 1.csv but not in 2.csv:")
        for smiles in missing_smiles:
            print(smiles)
    else:
        print("All SMILES in 1.csv are present in 2.csv.")

# 标准化SMILES
def standardize_smiles(smiles):
    try:
        # 从SMILES生成分子对象
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Invalid SMILES: {smiles}")
            return None

        # 对分子进行规范化处理
        standardized_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        return standardized_smiles
    except Exception as e:
        print(f"Error in processing SMILES {smiles}: {e}")
        return None


def process_SMILES_csv():
    # 读取CSV文件
    df = pd.read_csv('davis/davis_drug_smiles.csv')

    # 确保“Drug SMILES”列存在
    if 'Drug SMILES' not in df.columns:
        raise ValueError("The file does not contain a 'Drug SMILES' column.")

    # 对“Drug SMILES”列中的每个SMILES进行标准化
    df['Standardized SMILES'] = df['Drug SMILES'].apply(standardize_smiles)

    # 显示或保存结果
    print(df.head())  # 显示处理后的数据框的前几行
    # 将结果保存到新的CSV文件中
    df.to_csv('davis/standardized_drug_smiles.csv', index=False)

# 查看DB编号是否存在于异构信息中
def read_file_to_set(file_path, column_index, expected_num_columns):
    """读取文件并提取指定列的值到集合中"""
    values = set()
    with open(file_path, 'r') as file:
        for line_number, line in enumerate(file, 1):
            parts = line.split()
            if len(parts) < expected_num_columns:
                print(f"Skipping line {line_number}: expected {expected_num_columns} columns, found {len(parts)}")
                continue
            values.add(parts[column_index].strip())
    return values


def check_db_existence(file1_path, file2_path):
    # 从第一个文件中提取DB编号
    db_numbers_file1 = read_file_to_set(file1_path, 1, 2)  # 第二列DB编号，预期两列

    # 从第二个文件中提取DB编号
    db_numbers_file2 = read_file_to_set(file2_path, 0, 1)  # 唯一列DB编号，预期一列

    # 找出在两个文件中都存在的DB编号
    common_db_numbers = db_numbers_file1.intersection(db_numbers_file2)

    # 输出结果
    print("Common DB Numbers:")
    for db_number in sorted(common_db_numbers):
        print(db_number)
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

def filter_DrugID_without_smiles():
    # 定义药物编号
    drug_ids = ['DB00375', 'DB00994', 'DB01258']
    # 要处理的文件列表
    files = ['mat_data/drug_association/drug_disease.csv',
             'mat_data/drug_association/drug_drug.csv',
             'mat_data/drug_association/drug_protein.csv',
             'mat_data/drug_association/drug_se.csv',
             'mat_data/protein_association/protein_drug.csv']
    for file in files:
        # 读取 CSV 文件为 DataFrame
        df = pd.read_csv(file, index_col=0)

        # 行过滤：删除行索引与药物编号匹配的行
        df = df[~df.index.isin(drug_ids)]

        # 列过滤：删除列名称与药物编号匹配的列
        df = df.loc[:, ~df.columns.isin(drug_ids)]

        # 将过滤后的 DataFrame 写回文件
        df.to_csv(file)

    print("处理完成，相关行和列已删除。")


if __name__ == '__main__':
    # drug_in_hetero()
    # os.remove('mat_data/mat_drug_protein_remove_homo.txt')
    # davis_SMILES_in_DB()
    # process_SMILES_csv()
    # check_db_existence('davis/CIDDB.txt', 'mat_data/ID_drug.txt')
    # drug_smiles()
    # matrix_hang_lie_name()
    # disease_ID()
    # without_DBID()
    filter_DrugID_without_smiles()