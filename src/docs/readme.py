from os import path

def description(name,version):
    basename=path.basename(name)
    readme=f'''Name: {basename}
Discription:
    A flexible pipeline for metaviromic data analysis.
Usage:
    {name} <module> [opts] -o <output_file>
Author:
    yangming, yangm@idsse.ac.cn
Version:
    {version}, 2022-10-24 11:07'''
    return readme
