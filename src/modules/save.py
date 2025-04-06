def save_result(adata, result_file):
    adata.write(result_file)
    return result_file