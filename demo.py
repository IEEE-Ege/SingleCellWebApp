
import matplotlib 
import seaborn
import scanpy as sc
import numpy as np
import pandas as pd
from pca_analysis import PCAAnalyzer  # Yazdığın kodu içe aktar

# Örnek veri seti oluştur
adata = sc.datasets.pbmc3k()  # 3k PBMC hücre datası

# PCA işlemi için sınıfı başlat
pca_analyzer = PCAAnalyzer(adata)

# PCA Çalıştır
pca_analyzer.run_pca()

# PCA'yı Görselleştir
pca_analyzer.plot_pca(color="CST3")  # CST3, gen ifadesi için örnek bir sütun

# Varyans Oranını Çizdir
pca_analyzer.plot_variance()

# Sonuçları Kaydet
pca_analyzer.save_results()

# İşlenmiş veriye erişim
adata_processed = pca_analyzer.get_adata()
print("Processed AnnData object:", adata_processed)
