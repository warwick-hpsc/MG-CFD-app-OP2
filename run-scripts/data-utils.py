import pandas as pd
from sets import Set

def grep(text, filepath):
    found = False
    f = open(filepath, "r")
    for line in f:
        if text in line:
            found = True
            break
    return found
    
def clean_pd_read_csv(filepath):
    df = pd.read_csv(filepath, keep_default_na=False)
    return df

def get_data_colnames(df):
    mg_cfd_data_colnames = ["iters", "computeTime", "syncTime", "count"]
    op2_data_colnames = ["count", "total time", "plan time", "mpi time", "GB used", "GB total"]
    data_colnames = list(Set(mg_cfd_data_colnames+op2_data_colnames).intersection(Set(df.columns.values)))
    return data_colnames

def get_job_id_colnames(df):
    job_id_colnames = list(Set(df.columns.values).difference(Set(get_data_colnames(df))))

    if len(job_id_colnames) == 0:
        print("ERROR: get_job_id_colnames() failed to find any.")
        sys.exit(-1)
    return job_id_colnames