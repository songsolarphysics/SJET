# SJET: Solar Jet Extraction Tool

An interactive multi-algorithm solar jet extraction tool for quantitative analysis of solar jet phenomena.

## Overview

SJET (Solar Jet Extraction Tool) is a Python-based interactive tool designed for extracting and analyzing solar jets from FITS format observational data. The tool integrates five thresholding algorithms with morphological operations through a graphical user interface, and provides automated geometric parameter extraction including a novel Gaussian FWHM width measurement.

### Key Features

- **Five Thresholding Methods**: Manual, Otsu, Adaptive, Percentile, and Log-enhanced
- **Interactive Interface**: Real-time parameter adjustment with immediate visual feedback
- **ROI Selection**: Flexible region-of-interest selection (rectangle or polygon) for targeted analysis
- **Morphological Operations**: Opening and closing operations with distance-based and force-merge region merging
- **Geometric Parameter Extraction**: Automated calculation of jet length, width (boundary-based, area-based, and Gaussian FWHM), curvature, and deflection angles
- **Circular Region Analysis**: Novel start/end point identification method exploiting jet morphological asymmetry
- **Bézier Curve Modeling**: Quadratic curve fitting for jet axis representation, with manual control point override for complex morphologies
- **Gaussian FWHM Width**: Built-in cross-sectional Gaussian fitting on the original image for a threshold-independent width measure
- **Four-Panel Visualization**: Original data, binary mask, extracted jet, and edge detection
- **Data Export**: FITS, PNG, and ASCII metadata output with full traceability

## Quick Start

1. **Launch SJET**:
   ```bash
   python SJET.py
   ```

2. **Load FITS File**:
   - Click "Open FITS File" and select your data
   - Tool automatically handles NaN replacement and preprocessing

3. **Configure Analysis**:
   - Choose thresholding method and adjust parameters using sliders
   - Optionally define ROI (rectangle or polygon)
   - Apply morphological operations (opening/closing, region merging)

4. **Extract and Save**:
   - Monitor the four-panel visualization
   - Click "Save Results" to export FITS masks, PNG visualizations, and ASCII parameter records

## Geometric Parameter Analysis

Use the standalone analysis function:

```python
from Geometric_Parameter_Extraction import analyze_jet_circular_regions

# Basic usage (automatic control point)
results = analyze_jet_circular_regions(
    file_path='your_jet_mask.fits',
    data=original_image_array,   # numpy 2D array from sunpy_map.data
    visualize=True
)

print(f"Length:        {results['curve_length']:.2f} px")
print(f"Width (boundary): {results['average_width']:.2f} px")
print(f"Width (FWHM):  {results['fwhm_mean']:.2f} ± {results['fwhm_std']:.2f} px")
print(f"Rotation:      {results['rotation_angle_deg']:.2f} degrees")

# Manual control point override (for C-shaped or force-merged structures)
results = analyze_jet_circular_regions(
    file_path='your_jet_mask.fits',
    data=original_image_array,
    control_point_override=(row, col),  # specify in pixel coordinates
    visualize=True
)
```

> **Note**: The `data` parameter should be the raw intensity array (`sunpy.map.Map(...).data`) 
> from the **same observation frame** as the mask. It is used for Gaussian FWHM computation 
> and does not need to be normalised.

## Width Measurements

SJET provides three complementary width parameters:

| Parameter | Method | Description |
|---|---|---|
| `average_width` | Boundary-based | Mean perpendicular distance between mask edges at 10 locations |
| `average_width_by_area` | Area-based | Total mask area divided by curve length |
| `fwhm_mean` ± `fwhm_std` | Gaussian FWHM | Fitted on original image intensity; threshold-independent |

The Gaussian FWHM is recommended as a reproducible measure comparable to traditional slice-fitting methods. A large discrepancy between the boundary width and the FWHM may indicate over- or under-extraction and should prompt threshold re-evaluation.

## Supported Data

- **Solar Orbiter/EUI**: HRI_EUV 174 Å (tested with data release 6.0)
- **SDO/AIA**: Multi-wavelength (304 Å tested)
- **Format**: Any sunpy-compatible FITS file

## Thresholding Methods

| Method | Best suited for |
|---|---|
| Manual | High-contrast jets; direct control |
| Otsu | Images with bimodal intensity histograms |
| Adaptive | Non-uniform backgrounds; on-disk jets |
| Percentile | Long-tail intensity distributions |
| Log-enhanced | Low-contrast or faint jet structures |

> The block size for Adaptive thresholding is user-adjustable via an interactive slider. 
> Smaller values improve sensitivity to fine structures; larger values give smoother 
> segmentation. Adjust when switching between instruments with different spatial resolutions 
> (e.g., HRI at 195 km px⁻¹ vs. AIA at 435 km px⁻¹).

## Sensitivity and Limitations

- Endpoint identification is robust to moderate threshold variations (~80–90% percentile range).
- Area-based width and rotation angle are more sensitive to threshold choice; use Gaussian FWHM as a consistency check.
- For strongly curved, C-shaped, or force-merged structures, the automatic control point may deviate from the true jet spine. Use `control_point_override` to correct the Bézier axis.
- The large standard deviation of Gaussian FWHM across the ten cross-sections is a useful diagnostic: it may indicate that the fitted axis does not follow the jet spine consistently.
- Feature selection (which structure to analyse, at which time step) remains a manual decision that SJET does not automate.

## Data Access

- **Solar Orbiter/EUI**: https://soar.esac.esa.int/soar/
- **SDO/AIA**: http://jsoc.stanford.edu/

## Citation

If you use SJET in your research, please cite:

```bibtex
@article{tan2026sjet,
    title   = {{SJET}: An Interactive Solar Jet Extraction Tool},
    author  = {Tan, Song and Warmuth, Alexander and Schuller, Fr{\'e}d{\'e}ric 
               and Shen, Yuandeng and Fang, Yue and Mitchell, Jake A. J. 
               and Liu, Zedong},
    journal = {RAS Techniques and Instruments},
    year    = {2026},
    note    = {Accepted 2026 April 09}
}
```

## Contact

- **Issues**: Submit via GitHub Issues
- **Email**: stan@aip.de / tansong430@gmail.com
- **Institution**: Leibniz Institute for Astrophysics Potsdam (AIP)

## License

MIT License — see LICENSE file for details.

---

# SJET: 太阳喷流提取工具

## 概述

SJET（太阳喷流提取工具）是一个基于 Python 开发的交互式多算法太阳喷流特征提取工具，专门用于从 FITS 格式的太阳观测数据中提取和分析太阳喷流现象。工具集成了五种阈值算法与形态学操作，并新增了基于高斯拟合的 FWHM 宽度测量功能。

### 主要功能

- **五种阈值算法**: 手动阈值、Otsu 方法、自适应阈值、百分位数阈值和对数增强阈值
- **交互式界面**: 实时参数调整，即时可视化反馈
- **ROI 区域选择**: 矩形或多边形感兴趣区域选择
- **形态学操作**: 开运算、闭运算和区域合并功能（距离合并与强制合并）
- **几何参数提取**: 自动计算喷流长度、宽度（边界法、面积法、高斯 FWHM）、曲率和偏转角
- **圆形区域分析**: 基于形态不对称性的创新起始终点识别方法
- **贝塞尔曲线建模**: 二次曲线拟合喷流轴线，支持手动控制点覆盖
- **高斯 FWHM 宽度**: 基于原始图像强度的截面高斯拟合，与阈值选择无关
- **四面板可视化**: 同时显示原始数据、二值掩膜、提取的喷流和边缘检测结果
- **数据导出**: 支持 FITS、PNG 和 ASCII 格式输出，保留完整元数据

## 快速开始

1. **启动 SJET**:
   ```bash
   python SJET.py
   ```

2. **加载 FITS 文件**:
   - 点击 "Open FITS File" 按钮选择数据文件
   - 工具自动处理 NaN 值替换和数据预处理

3. **配置分析参数**:
   - 选择阈值方法并使用滑块调整参数
   - 可选择定义感兴趣区域（矩形或多边形）
   - 应用形态学操作

4. **提取和保存结果**:
   - 监控四面板可视化更新
   - 点击 "Save Results" 导出结果

## 几何参数分析

```python
from Geometric_Parameter_Extraction import analyze_jet_circular_regions

# 基本用法（自动控制点）
results = analyze_jet_circular_regions(
    file_path='your_jet_mask.fits',
    data=original_image_array,   # sunpy_map.data 的 numpy 二维数组
    visualize=True
)

print(f"曲线长度: {results['curve_length']:.2f} px")
print(f"边界宽度: {results['average_width']:.2f} px")
print(f"FWHM 宽度: {results['fwhm_mean']:.2f} ± {results['fwhm_std']:.2f} px")
print(f"旋转角:   {results['rotation_angle_deg']:.2f} 度")

# 手动控制点覆盖（适用于 C 形或强制合并后的复杂形态）
results = analyze_jet_circular_regions(
    file_path='your_jet_mask.fits',
    data=original_image_array,
    control_point_override=(row, col),
    visualize=True
)
```

## 宽度参数说明

| 参数 | 方法 | 说明 |
|---|---|---|
| `average_width` | 边界法 | 10 个截面处掩膜边界垂直距离的均值 |
| `average_width_by_area` | 面积法 | 掩膜总面积除以曲线长度 |
| `fwhm_mean` ± `fwhm_std` | 高斯 FWHM | 基于原始图像强度拟合，与阈值无关 |

推荐同时报告边界宽度和 FWHM，两者物理含义互补。FWHM 与传统切片拟合方法直接可比。

## 敏感性与局限性

- 端点识别对中等阈值变化（80–90% 百分位范围）较为鲁棒。
- 面积宽度和旋转角对阈值更敏感，建议用 FWHM 作为一致性检验。
- 对于强弯曲、C 形或强制合并后的结构，自动控制点可能偏离真实轴线，此时建议使用 `control_point_override` 手动修正。
- 目标结构的选择（分析哪个特征、哪个时刻）仍是用户的主观判断，SJET 不自动化这一决策。

## 数据获取

- **Solar Orbiter/EUI**: https://soar.esac.esa.int/soar/
- **SDO/AIA**: http://jsoc.stanford.edu/

## 引用格式

```bibtex
@article{tan2026sjet,
    title   = {{SJET}: An Interactive Solar Jet Extraction Tool},
    author  = {Tan, Song and Warmuth, Alexander and Schuller, Fr{\'e}d{\'e}ric 
               and Shen, Yuandeng and Fang, Yue and Mitchell, Jake A. J. 
               and Liu, Zedong},
    journal = {RAS Techniques and Instruments},
    year    = {2026},
    note    = {Accepted 2026 April 09}
}
```

## 联系方式

- **问题反馈**: 通过 GitHub Issues 提交
- **联系邮箱**: stan@aip.de / tansong430@gmail.com
- **机构**: 莱布尼茨天体物理研究所 (AIP)

## 许可证

MIT 许可证 — 详见 LICENSE 文件
