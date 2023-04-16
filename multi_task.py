import warnings
import anndata
import scanpy as sc
import Cell_BLAST as cb
import os

warnings.simplefilter("ignore")
cb.config.RANDOM_SEED = 0

# data_dir = "/home/cnic02/wangyn/data/figshare"
# fat = anndata.read_h5ad(os.path.join(data_dir, 'Fat_droplet.h5ad'))
data_dir = "/home/cnic02/wangyn/data/figshare"
fat = anndata.read_h5ad(os.path.join(data_dir, 'Fat_droplet.h5ad'))

axes = cb.data.find_variable_genes(fat)
# model = cb.directi.fit_DIRECTi(
#     fat, genes=fat.var.query("variable_genes").index, supervision='cell_ontology_class',
#     latent_dim=100, cat_dim=7
# )
# model.save("/home/cnic02/wangyn/projects/Cell_BLAST/models/directi_fat")
# del model
model = cb.directi.DIRECTi.load("/home/cnic02/wangyn/projects/Cell_BLAST/models/directi_fat")
# fat.obsm["X_gau"],fat.obsm["X_cat"], fat.obsm["X_latent"]
fatty= model.inference(fat)
print(len(fatty))
sc.pp.neighbors(fat, use_rep="X_latent")
sc.tl.umap(fat)
sc.pl.umap(fat, color="cell_ontology_class", palette="tab20")