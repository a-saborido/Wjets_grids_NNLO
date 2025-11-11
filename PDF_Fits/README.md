# README

This folder contains tools to produce PDF plots from **xFitter** outputs.

---

## Requirements

- **xFitter** 
- **ROOT**  
- **LHAPDF** (only required for LHAPDF overlays)

To set up the environment:

```bash
source ~/xfitter/tools/setup.sh
```

---

## Generate xFitter output and ROOT files

Run the xFitter fit and create ROOT-compatible PDF output:
```bash
source ~/xfitter/tools/setup.sh
nohup xfitter &
xfitter-draw output --bands --root --q2all
# or restrict to a single scale:
xfitter-draw output --bands --root --q2 10
```

ROOT files are saved as:
```
output/plots.root
```

To list available graphs:
```bash
rootls output/plots.root:Graphs
```

---

## Quick plotting with `plot_pdf.C`

The script **`plot_pdf.C`** provides a simple way to compare two PDF fits.

Example:
```bash
root -l -b -q 'plot_pdf.C("uv",1.9,"output_HERA+DY/plots.root","uv","output_HERA+DY+absyj1/plots.root")'
```

**Explanation of arguments:**
1. `"uv"` — the PDF flavour to plot (e.g. `uv`, `dv`, `g`, `sbar`, etc.)  
2. `1.9` — the Q2 scale in GeV2  
3. `"output_HERA+DY/plots.root"` — reference ROOT file (baseline fit)  
4. `"uv"` — the PDF flavour in the comparison file (choose same as first)  
5. `"output_HERA+DY+absyj1/plots.root"` — comparison ROOT file  

The first output argument is the reference PDF. The second is divided by the central value of the reference and displayed in a ratio plot.

---

## Using the more complete script: `plot_pdf_thesis.C`

The script **`plot_pdf_thesis.C`** extends `plot_pdf.C`. It adds support for LHAPDF comparisons and uncertainty breakdowns.

### Key Features
- Supports both **replica** and **(a)symmetric Hessian** error sets
- Adds parameterisation uncertainty bands for total uncertainty visualization
- Provides **three panels** in one canvas:
  1. **Main panel:** absolute PDF distributions with uncertainty bands
  2. **Ratio panel:** comparison between PDFs or to a reference
  3. **Relative-uncertainty panel:** fractional uncertainties for each dataset

---

### How to use it

- Run with a command like:
   ```bash
   root -l -b -q 'plot_pdf_thesis.C("uv",10,
       "output_HERA_wzProduction/plots.root","uv",
       "output_HERA_wzProduction_WpJ+WmJ_ptw/plots.root",
       "ATLASepWZVjet20-EIG")'
   ```

   This command:
   - Plots **uv** distributions at **Q² = 10 GeV²**
   - Compares two ROOT files (baseline and comparison)
   - Overlays the **ATLASepWZVjet20-EIG** LHAPDF set for reference

- Optionally include parameterisation uncertainties to show augmented total bands:
   ```bash
   root -l -b -q 'plot_pdf_thesis.C("uv",10,
       "output_HERA_wzProduction/plots.root","uv",
       "output_HERA_wzProduction_WpJ+WmJ_ptw/plots.root",
       "ATLASepWZVjet20-EIG",
       "output_HERA_wzProduction_NegativeGluon/plots.root",
       "output_HERA_wzProduction_WpJ+WmJ_ptw_NegativeGluon/plots.root")'
   ```

   This option incorporates the difference between PDFs obtained using the two parameterisations as an additional source of uncertainty in the total error band.

---

### Example for testing a reference

```bash
root -l -b -q 'plot_pdf_thesis.C("sbar",10,
    "output_HERA_wzProduction_NegativeGluon/plots.root","sbar",
    " ","ATLAS-epWZ16-EIG",
    "output_HERA_WZ/plots.root"," ")'
```

---

### Output

The script automatically produces a PDF file named like:
```
pdf_<flavour>_Q2_<value>_withLHAPDF_relunc.pdf
```
For example:
```
pdf_uv_Q2_10_withLHAPDF_relunc.pdf
```

