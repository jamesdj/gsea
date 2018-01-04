import os
import subprocess

import pandas as pd

gsea_dir = os.path.abspath(os.path.dirname(__file__))


def write_gmt(sets_to_genes, gmt_file, append=False):
    strings = []
    for s, genes in sets_to_genes.items():
        try:
            string = '{}\t\t{}'.format(s, '\t'.join(map(str, genes)))
        except Exception as e:
            print(s)
            print(genes)
            raise(e)
        strings.append(string)
    with open(gmt_file, 'a' if append else 'w') as f:
        f.write('\n'.join(strings))
        f.write('\n')


def read_gct(gct_path, use_description=False, description_func=lambda x: x, transpose=True, dropna_axes=[]):
    df = pd.read_table(gct_path, skiprows=2, index_col=0)
    for axis in dropna_axes:
        df = df.dropna(axis=axis)
    if use_description:
        df.index = [description_func(d) for d in df.Description]
    df = df.drop("Description", axis=1)
    df = df.dropna(axis=0)
    if transpose:
        df = df.T
    return df


def write_gct(tdf, path):
    df = tdf.T
    n_rows, n_cols = df.shape
    with open(path, 'w') as f:
        f.write("#1.2\n")
        f.write("{}\t{}\n".format(n_rows, n_cols))
        header_row = 'Name\tDescription\t' + '\t'.join(df.columns) + '\n'
        f.write(header_row)
        for i in range(n_rows):
            name = df.index[i]
            values = df.iloc[i]
            row = '{}\t{}\t{}\n'.format(name, name, '\t'.join(map(str, values)))
            f.write(row)


def r_ssgsea(exp_df, genesets, wd='/tmp', weight=0.75, return_path=False, verbose=False):
    if not os.path.exists(wd):
        os.mkdir(wd)
    exp_file = os.path.join(wd, 'exp.gct')
    write_gct(exp_df, exp_file)
    genesets_file = os.path.join(wd, 'genesets.gmt')
    if len(genesets) == 1:
        only_gs_name, only_gs = list(genesets.items())[0]
        genesets[only_gs_name + ".1"] = only_gs
    write_gmt(genesets, genesets_file)
    out_file = os.path.join(wd, 'ssgsea_out.gct')

    r_ssgsea_script_path = os.path.join(gsea_dir, 'run_ssgsea.R')
    ccba_source_path = os.path.join(gsea_dir, 'CCBA.ssgsea.R')
    command = 'Rscript {} {} {} {} {} {}'.format(r_ssgsea_script_path, exp_file, genesets_file, weight, out_file, ccba_source_path)
    print(command)
    #os.system(command)
    try:
        cmnd_output = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True, universal_newlines=True);
        if verbose:
            print(cmnd_output)
    except subprocess.CalledProcessError as exc:
        print("Status : FAIL", exc.returncode, exc.output)
    else:
        print("Output: \n{}\n".format(cmnd_output))


    if return_path:
        return out_file
    result_df = read_gct(out_file)
    #print(result_df.head())
    return result_df