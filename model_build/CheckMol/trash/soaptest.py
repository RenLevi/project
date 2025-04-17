from ase.io import read
from dscribe.descriptors import SOAP
from dscribe.kernels import AverageKernel
import numpy as np

# 创建示例原子结构
mol1 = read('cal\\nequipOpt.traj')
mol2 = read('cal\\POSCAR')
# 设置SOAP描述符参数
soap = SOAP(
    species=[ 'Ru','H',"C", "O"],
    periodic=True,
    r_cut=5.0,
    n_max=8,
    l_max=6,
)

# 计算SOAP描述符
soap1 = soap.create(mol1)
soap2 = soap.create(mol2)


from sklearn.metrics.pairwise import cosine_similarity
similarity_matrix = cosine_similarity(soap1, soap2)

print("原子环境相似度矩阵:")
print(similarity_matrix)

# 可以计算平均相似度作为整体分子相似度
mean_similarity = np.mean(similarity_matrix)
print(f"平均原子环境相似度: {mean_similarity:.4f}")