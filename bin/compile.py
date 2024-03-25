import py_compile
import os

# 获取当前目录下的所有Python文件
python_files = [f for f in os.listdir('.') if f.endswith('.py')]

# 遍历所有Python文件并编译成.pyc文件
for python_file in python_files:
    try:
        py_compile.compile(python_file)
        print(f"{python_file} 编译成功！")
    except py_compile.PyCompileError as e:
        print(f"{python_file} 编译失败：{e}")
