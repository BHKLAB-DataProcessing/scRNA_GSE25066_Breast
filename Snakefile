from os import path
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"],
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)

prefix = config["prefix"]
filename = config["filename"]

rule get_object:
    output:
        S3.remote(prefix + filename)
    resources:
        mem_mb = 8000,
        disk_mb = 10000
    shell:
        """
        Rscript scripts/wrapper_GSE25066-2.R {prefix} {filename}
        """
