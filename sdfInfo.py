import pandas as pd
from rdkit import Chem


def extract_info_from_sdf_to_csv():
    # 读取SDF文件，不进行消毒
    supplier = Chem.SDMolSupplier('DrugBank/DB_Drug_structures.sdf', sanitize=False)

    # 存储提取的数据
    data = []

    # 遍历分子
    for mol in supplier:
        if mol is None:  # 跳过无法解析的分子
            continue

        try:
            # 提取所需的信息，去掉尖括号和空格
            smiles = mol.GetProp("SMILES") if mol.HasProp("SMILES") else "N/A"
            formula = mol.GetProp("FORMULA") if mol.HasProp("FORMULA") else "N/A"
            drugbank_id = mol.GetProp("DRUGBANK_ID") if mol.HasProp("DRUGBANK_ID") else "N/A"
            generic_name = mol.GetProp("GENERIC_NAME") if mol.HasProp("GENERIC_NAME") else "N/A"

            # 调试输出
            print(f"SMILES: {smiles}, FORMULA: {formula}, DRUGBANK_ID: {drugbank_id}, GENERIC_NAME: {generic_name}")

            # 将提取的信息添加到数据列表中
            data.append({
                "SMILES": smiles,
                "FORMULA": formula,
                "DRUGBANK_ID": drugbank_id,
                "GENERIC_NAME": generic_name
            })
        except Exception as e:
            # 如果在提取过程中出现错误，则打印错误并继续
            print(f"Error processing molecule: {e}")

    # 将数据转换为DataFrame
    df = pd.DataFrame(data)

    # 保存为CSV文件
    df.to_csv('DrugBank/DB_Drug.csv', index=False)


if __name__ == '__main__':
    # 调用函数进行提取和保存
    extract_info_from_sdf_to_csv()
