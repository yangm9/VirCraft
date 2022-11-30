from os import path

def description(name,version):
    basename=path.basename(name)
    readme=f'''Name: {basename}
Discription:
    A flexible pipeline for metaviromic data analysis.
Usage:
    {name} <module> [opts] -o <outdir>
        module: a optional functional module, including assembly, identify, votus, classify, quantify, func_annot and host_prid.
        opts: options described below in the section of Options.
        outdir: output directory.
Author:
    yangming, yangm@idsse.ac.cn
Version:
    {version}, 2022-10-24 11:07'''
    return readme
