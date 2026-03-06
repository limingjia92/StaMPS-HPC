# StaMPS-HPC: GMT Color Palette Tables (`cptfiles`)

**Version:** 1.0.0  
**Maintainer:** Mingjia Li  
**Context:** StaMPS-HPC Optimization Project  
**Date:** March 2026

---

## Module Overview

This directory contains approximately 30 **Generic Mapping Tools (GMT)** Color Palette Table (`.cpt`) files. 

These files are included to significantly expand and enhance the default MATLAB colormap options for visualizing InSAR processing results (e.g., deformation rates, topography, coherence, and phase). By leveraging standard GMT color scales, StaMPS-HPC can produce publication-quality, visually intuitive maps that are standard in the geophysics community.

Within the StaMPS-HPC framework, these files are seamlessly read and parsed by the custom MATLAB function `cptcmap.m` and are directly applied during the execution of plotting commands like `ps_plot.m`.

---

## Usage & Configuration

You do not need to manually edit or compile any files in this directory. 

To change the colormap used for your StaMPS plots, you can globally set the color scheme parameter in your MATLAB command window before calling `ps_plot.m`. 

### Example Command
To set the color scheme to the `GMT_relief.cpt` palette, use the following command (omitting the `.cpt` extension):

```matlab
% Set the plotting colormap to GMT_relief
setparm('plot_color_scheme', 'GMT_relief')

% Generate the plot (e.g., mean velocity)
ps_plot('v')
```

### Common Colormap Recommendations

While you can explore all 30+ files, here are a few commonly used options for InSAR data:

 * **GMT_red2green / GMT_seis** : Excellent for diverging data like deformation velocity (where 0 is a neutral color, positive is blue/cool, and negative is red/warm).

 * **GMT_relief / GMT_globe** : Ideal for topographic (DEM) mapping.

 * **GMT_wysiwyg** : Useful for wrapping phase visualization.

### Adding Custom Colormaps

You are welcome to supplement this directory with your own custom GMT .cpt files to further expand your visualization options. If you choose to add new palettes, please ensure that their internal formatting is standard and can be correctly parsed by the cptcmap.m function.

---

## Important Note

Do not delete this directory. If the specified .cpt file is missing, the ps_plot.m function will fail to locate the requested color scheme via cptcmap.m and may result in an error or force a fallback to suboptimal default MATLAB colormaps (like jet).