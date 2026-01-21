import os
import numpy as np
from Bio import PDB

def calculate_rg(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    
    # 提取所有 CA 原子的坐标
    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    coords.append(residue['CA'].get_coord())
    
    coords = np.array(coords)
    
    if len(coords) == 0:
        return None

    # 计算几何中心 (Center of Geometry)
    center = np.mean(coords, axis=0)
    
    # 计算到中心的平方距离之和
    diff = coords - center
    squared_distances = np.sum(diff**2, axis=1)
    
    # 计算 Rg
    rg = np.sqrt(np.mean(squared_distances))
    return rg

def process_folder(folder_path):
    print(f"{'PDB_File':<20} | {'Rg (Å)':<10}")
    print("-" * 35)
    
    # 遍历文件夹中的 pdb 文件
    for filename in sorted(os.listdir(folder_path)):
        if filename.endswith(".pdb"):
            file_path = os.path.join(folder_path, filename)
            try:
                rg = calculate_rg(file_path)
                if rg is not None:
                    print(f"{filename:<20} | {rg:.4f}")
                else:
                    print(f"{filename:<20} | No CA atoms found")
            except Exception as e:
                print(f"{filename:<20} | Error: {e}")

# 使用示例：将 '.' 替换为你的文件夹路径
process_folder('.')
