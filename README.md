# VirCraft
VirCraft is an automatic viromic analysis pipeline.

## 1 软件结构
![Overall workflow of VirCraft](docs/Overall_workflow_of_VirCraft.png)

## 2 安装和使用教程

```
sh install.sh
```

## 3 结果文件说明

## 5 注意事项
当前版本只能生成脚本并不能直接运行，请生成脚本后自行运行。
## 6 版本更新日志


**VirCraft-v0.0.1版**
```
初始版本。
```

**VirCraft-v0.0.2版**
```
初始版本。
```

**VirCraft-v0.0.3版**
```
尝试多路分析失败的版本。
```

**VirCraft-v0.0.4版**
```
1.暂时不支持多线程的版本。
2.quantify模块完成，包括统计多样性分析、散点图、柱状图等。
3.大势所趋，用import argparse代替from optparse import OptionParser。
```

**VirCraft-v0.0.5版**
```
添加AMGs分析，暂定用dramv和vibrant分析，尚未整合。
```

**VirCraft-v0.0.6版**
```
1.手动curation部分的自动化。
2.整合compare结果和classify结果。
3.完善host_pred模块。
```
