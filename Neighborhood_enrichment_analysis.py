import scanpy as sc
import squidpy as sq
import spatialdata as sd
import matplotlib.pyplot as plt

zarr_path = "/xenium.zarr"
sdata = sd.read_zarr(zarr_path)
adata = sdata.tables["table"]

print(adata)
marker_genes = [
  "KLRC1", "NCAM1", "KLRK1",
  "CD14", "CD68", "CD163",
  "SNHG15", "SOST", "KIT",
  "MZB1", "XBP1", "FCRL5",
  "CD3E", "CD8A", "CD2",
  "COL5A1", "ANXA2", "MMP14",
  "PECAM1", "POSTN", "MYH9",
  "MS4A1", "CD19", "CD79A",
  "COL4A1", "COL4A2", "COL18A1",
  "TUBB", "STMN1", "GZMB"
]

# Dot plot 
sc.pl.dotplot(
  adata,
  var_names=marker_genes,
  groupby="cell_type",
  standard_scale="var",
  dot_max=0.7,
  cmap="RdBu_r",
  title="Marker gene expression across annotated cell types"
)

# Neighborhood enrichment analysis
sq.gr.spatial_neighbors(
  adata,
  spatial_key="spatial"
)

sq.gr.nhood_enrichment(
  adata,
  cluster_key="leiden"
)

sq.pl.nhood_enrichment(
  adata,
  cluster_key="leiden",
  cmap="Blues",
  title="Neighborhood enrichment of Leiden clusters"
)
