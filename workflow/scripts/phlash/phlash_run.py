import sys
import os
import argparse

# Get the number of threads from environment variable or default to 1
threads = os.environ.get("OMP_NUM_THREADS", "1")

# JaX configuration
os.environ["XLA_FLAGS"] = (
    f"--xla_cpu_multi_thread_eigen=false "
    f"--xla_force_host_platform_device_count={threads}"
)

# Do not use async dispatch to avoid issues with thread management in some environments
os.environ["JAX_FLAGS"] = (
    f"--jax_cpu_enable_async_dispatch=false "
    f"--jax_backend_target=cpu"
)


# Math libraries
os.environ["OMP_NUM_THREADS"] = threads
os.environ["MKL_NUM_THREADS"] = threads
os.environ["OPENBLAS_NUM_THREADS"] = threads
os.environ["VECLIB_MAXIMUM_THREADS"] = threads
os.environ["NUMEXPR_NUM_THREADS"] = threads
os.environ["BLIS_NUM_THREADS"] = threads
os.environ["NUMBA_NUM_THREADS"] = threads

# Memory management
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"

# Rest of them

import phlash
import pysam
import numpy as np
import pandas as pd


def get_chrom_length(vcf_obj, chrom):

    if chrom in vcf_obj.header.contigs:
        length = vcf_obj.header.contigs[chrom].length
        if length:
            return length

    last_pos = 1
    try:
        for record in vcf_obj.fetch(chrom):
            last_pos = record.pos
        return last_pos
    except Exception as e:
        raise ValueError(f"Could not determine length for chromosome {chrom}: {e}")



def run_phlash_analysis(args):
    
    with open(args.chroms, 'r') as f:
        chrom_names = [line.strip() for line in f if line.strip()]
    
    with open(args.samples_file, 'r') as f:
        sample_ids = [line.strip() for line in f if line.strip()]
    
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"El archivo BCF no existe en la ruta: {args.input}")

    vcf_obj = pysam.VariantFile(args.input)
    
    chroms_to_fit = []

    for chrom in chrom_names:
        if chrom not in vcf_obj.header.contigs:
            print(f"Warning: Cromosome {chrom} not found in BCF header. Skipping.")
            continue

        chrom_len = get_chrom_length(vcf_obj, chrom)
        print(f"Preparando {chrom}: longitud 1-{chrom_len} con {len(sample_ids)} muestras.")

        chroms_to_fit.append(
            phlash.contig(
                args.input,
                samples=sample_ids,
                region=f"{chrom}:1-{chrom_len}"
            )
        )
    
    vcf_obj.close()

    
    print(f"Running phlash.fit (mu={args.mu})...")
    results = phlash.fit(chroms_to_fit, mutation_rate=args.mu)

    rescaled_results = [dm.rescale(args.mu) for dm in results]
    
    times = np.array([dm.eta.t[1:] for dm in rescaled_results])
    
    T = np.geomspace(times.min(), times.max(), 1000)
    
    data_to_save = {"time_generations": T}
    
    for i, dm in enumerate(rescaled_results):
        data_to_save[f"n_{i+1}"] = dm.eta(T, Ne=True)

    df_output = pd.DataFrame(data_to_save)
    df_output.to_csv(args.output, sep='\t', index=False)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Análisis demográfico con phlash usando un único BCF")
    parser.add_argument("-i", "--input", required=True, help="Ruta al archivo BCF único")
    parser.add_argument("-c", "--chroms", required=True, help="Archivo .txt con los cromosomas a incluir")
    parser.add_argument("-s", "--samples_file", required=True, help="Archivo .txt con IDs de individuos")
    parser.add_argument("-m", "--mu", type=float, required=True, help="Mutation rate (ej: 1.25e-8)")
    parser.add_argument("-o", "--output", required=True, help="Archivo TSV de salida")

    args = parser.parse_args()
    run_phlash_analysis(args)
