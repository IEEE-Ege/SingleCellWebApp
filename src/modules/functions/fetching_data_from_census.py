!pip install cellxgene_census
import cellxgene_census 
from typing import List, Optional
import pandas as pd
import anndata

def census_function(
    organism: str = "Homo sapiens",
    cell_filter: Optional[str] = None,
    gene_filter: Optional[str] = None,
    metadata_columns: Optional[List[str]] = None,
    return_anndata: bool = False,
    census_version: str = "stable"
) -> pd.DataFrame | anndata.AnnData:

    if metadata_columns is None:        # if there are no parameters set from the user just get the parameters below.
        metadata_columns = [
            "dataset_id","assay","suspension_type",
            "sex","tissue_general","tissue","cell_type"
        ]

    with cellxgene_census.open_soma(census_version=census_version) as census: #connecting to the census dataset and census version is to decide which census we are fetching the data
        if return_anndata:   # if true the get_anndata function is called
            adata = cellxgene_census.get_anndata(
                census=census,
                organism=organism,
                var_value_filter=gene_filter,
                obs_value_filter=cell_filter,
                column_names={"obs": metadata_columns},
            )
            return adata
        else:
            obs = census["census_data"][organism.lower().replace(" ", "_")].obs # we are getting the cell informations about the given organism like "homo sapiens"
            table = obs.read(
                value_filter=cell_filter,
                column_names=metadata_columns
            )
            data = table.concat().to_pandas()
            return data # the data is like an excel table

#testing the function
if __name__ == "__main__":
    # example usage
    print("Test is running...")
    
    # test with assumed parameters
    default_result = census_function()
    print("\nResult with default parameters (DataFrame):")
    print(default_result.head())
    
    adata = census_function(
        organism="Homo sapiens",
        cell_filter="tissue == 'brain' and sex == 'male'",
        gene_filter="feature_id in ['ENSG00000161798','ENSG00000188229']",
        metadata_columns=["cell_type","tissue"],
        return_anndata=True
    )
    print("\nAnndata object:")
    print(adata)
