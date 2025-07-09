# SJET: Solar Jet Extraction Tool

An interactive multi-algorithm solar jet extraction tool for quantitative analysis of solar jet phenomena using high-resolution solar observations.

## Overview

SJET (Solar Jet Extraction Tool) is a Python-based tool designed for extracting and analyzing solar jets from FITS format observational data. The tool integrates multiple thresholding algorithms with morphological operations through an interactive graphical user interface.

### Key Features

- **Five Thresholding Methods**: Manual, Otsu, Adaptive, Percentile, and Log-enhanced
- **Interactive Interface**: Real-time parameter adjustment with immediate visual feedback
- **ROI Selection**: Flexible region-of-interest selection for targeted analysis
- **Morphological Operations**: Opening and closing operations with region merging
- **Geometric Parameter Extraction**: Automated calculation of jet length, width, curvature, and deflection angles
- **Circular Region Analysis**: Novel start-end point identification method
- **Bézier Curve Modeling**: Quadratic curve fitting for jet axis representation
- **Four-Panel Visualization**: Original data, binary mask, extracted jet, and edge detection
- **Data Export**: FITS and PNG format output with metadata preservation

## Quick Start

1. **Launch SJET**:
   ```bash
   python SJET.py
   ```

2. **Load FITS File**:
   - Click "Open FITS File" and select your data
   - Tool automatically handles preprocessing

3. **Configure Analysis**:
   - Choose thresholding method
   - Adjust parameters using sliders
   - Optionally define ROI
   - Apply morphological operations

4. **Extract and Save**:
   - Monitor visualization updates
   - Click "Save Results" to export

## Geometric Parameter Analysis

Use the standalone analysis function:

```python
from Geometric_Parameter_Extraction import analyze_jet_circular_regions

results = analyze_jet_circular_regions(
    file_path='your_jet_mask.fits',
    visualize=True,
    save_results=True
)

print(f"Length: {results['length']:.2f} pixels")
print(f"Width: {results['average_width']:.2f} pixels")
print(f"Rotation: {results['rotation_angle_deg']:.2f} degrees")
```

## Supported Data

- **Solar Orbiter/EUI**: HRI_EUV 174 Å
- **SDO/AIA**: Multi-wavelength (304 Å tested)
- **Format**: sunpy compatible imaging files

## Thresholding Methods

1. **Manual**: User-controlled percentage thresholding
2. **Otsu**: Automatic histogram-based threshold
3. **Adaptive**: Local neighborhood thresholding
4. **Percentile**: Statistical distribution-based
5. **Log-enhanced**: Logarithmic transformation for weak signals

## Data Access

- **Solar Orbiter**: https://soar.esac.esa.int/soar/
- **SDO/AIA**: http://jsoc.stanford.edu/

## Citation

```bibtex
@article{tan2025sjet,
    title={SJET: An Interactive Multi-Algorithm Solar Jet Extraction Tool},
    author={Tan, Song and Warmuth, Alexander and Schuller, Frédéric and Shen, Yuandeng and Mitchell, Jake A. J. and Liu, Zedong},
    journal={xxx},
    year={2025},
    note={Submitted}
}
```

## Contact

- **Issues**: Submit via GitHub Issues
- **Contact**: awarmuth@aip.de
- **Institution**: Leibniz Institute for Astrophysics Potsdam (AIP)

## License

MIT License - see LICENSE file for details.

---

# SJET: 太阳喷流提取工具

## 概述

SJET（太阳喷流提取工具）是一个基于Python开发的交互式多算法太阳喷流特征提取工具，专门用于从FITS格式的太阳观测数据中提取和分析太阳喷流现象。

### 主要功能

- **五种阈值算法**: 手动阈值、Otsu方法、自适应阈值、百分位数阈值和对数增强阈值
- **交互式界面**: 实时参数调整，即时可视化反馈
- **ROI区域选择**: 灵活的感兴趣区域选择功能
- **形态学操作**: 开运算、闭运算和区域合并功能  
- **几何参数提取**: 自动计算喷流长度、宽度、曲率和偏转角
- **圆形区域分析**: 创新的起始-终点识别方法
- **贝塞尔曲线建模**: 二次曲线拟合喷流轴线
- **四面板可视化**: 同时显示原始数据、二值掩膜、提取的喷流和边缘检测结果
- **数据导出**: 支持FITS和PNG格式输出，保留完整元数据


## 快速开始

1. **启动SJET**:
   ```bash
   python SJET.py
   ```

2. **加载FITS文件**:
   - 点击"Open FITS File"按钮选择数据文件
   - 工具自动处理数据预处理

3. **配置分析参数**:
   - 选择阈值方法
   - 使用滑块调整参数
   - 可选择定义感兴趣区域
   - 应用形态学操作

4. **提取和保存结果**:
   - 监控可视化更新
   - 点击"Save Results"导出结果

## 几何参数分析

使用独立的分析函数：

```python
from Geometric_Parameter_Extraction import analyze_jet_circular_regions

results = analyze_jet_circular_regions(
    file_path='your_jet_mask.fits',
    visualize=True,
    save_results=True
)

print(f"长度: {results['length']:.2f} 像素")
print(f"宽度: {results['average_width']:.2f} 像素") 
print(f"旋转角: {results['rotation_angle_deg']:.2f} 度")
```

## 支持的数据类型

- **Solar Orbiter/EUI**: HRI_EUV 174 Å观测数据
- **SDO/AIA**: 多波长数据（已测试304 Å）
- **格式**: sunpy兼容的标准FITS文件

## 阈值方法说明

1. **手动阈值**: 用户控制的百分比阈值
2. **Otsu方法**: 基于直方图的自动阈值
3. **自适应阈值**: 局部邻域阈值处理
4. **百分位数**: 基于统计分布的阈值
5. **对数增强**: 对数变换增强弱信号

## 数据获取

- **Solar Orbiter数据**: https://soar.esac.esa.int/soar/
- **SDO/AIA数据**: http://jsoc.stanford.edu/

## 引用格式

```bibtex
@article{tan2025sjet,
    title={SJET: An Interactive Multi-Algorithm Solar Jet Extraction Tool},
    author={Tan, Song and Warmuth, Alexander and Schuller, Frédéric and Shen, Yuandeng and Mitchell, Jake A. J. and Liu, Zedong},
    journal={xxx},
    year={2025},
    note={Submitted}
}
```

## 联系方式

- **问题反馈**: 通过GitHub Issues提交
- **联系邮箱**: awarmuth@aip.de  
- **机构**: 莱布尼茨天体物理研究所 (AIP)

## 许可证

MIT许可证 - 详见LICENSE文件
