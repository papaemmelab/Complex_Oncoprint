# generate_complex_oncoprint

**generate_complex_oncoprint** is an R function that utilizes the **jokergoo/ComplexHeatmap** package to create **oncoprints** for visualizing variants, including **mutations** (required), **copy number variants (CNVs)** (optional), and **structural variants (SVs)** (optional). The function provides publication-ready visuals with easily adjustable aesthetics.

The colors follow mutation effects based on the **Isabl pipeline default formats** (Elli Papaemanuil’s Lab).

---

## 🚀 Installation

To install **ComplexHeatmap**, ensure you have **Bioconductor** installed:

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
```

---

## 📌 Usage

### Simple Example

```r
source("./generate_complex_oncoprint.R")

ht <- generate_complex_oncoprint(
    muts = MUT, 
    show.sample.names = TRUE, 
    show.border = FALSE,
    min.freq = 1, 
    show.title = TRUE, 
    title.str = "Example 1 - YOHOOOO oncoprint with NO ERROR!", 
    save.name = "Example_1.A",
    save.path = "./example_oncoprints/"
)
```

---

## 📊 Features

- **Supports**:
  - Mutations (required)
  - Copy Number Variants (**CNVs**, optional)
  - Structural Variants (**SVs**, optional)
- **Customizable aesthetics** for publication-quality figures.
- **Annotation banners**: Add categorical variables from a lookup table as **bottom annotations**:
  ```r
  show.another.banner = TRUE
  banner.name = c("FEATURE.COL.1", "FEATURE.COL.2", ...)
  ```
- **Paired sample oncoprint** for pre/post visualization (*experimental, may throw errors*).
- **Timepoint visualization**: Option to show a **heatbar** for samples from the same patient.
- **Annotation barplots**: Examples available in `BOTTOM_annot` script for specific projects (**not yet integrated** into parameter-based runs).

---

## 📁 File Structure

- **Example script & documentation**: Refer to example files for more details on function usage.
- **Recognized terms** for mutation, CNV, and SV types (column: `EFFECT`):  
  - Refer to: `~/Complex_Oncoprint/sub_function/test_required_fields.R`.

---

## ⚠️ Notes

- **Paired sample oncoprint** is available but **may need updates**.
- If you need a **specific mutation representation**, please reach out.

📩 Contact: **rahnaman@mskcc.org**

---

## 🙌 Acknowledgments

Thanks to **Zuguang Gu** for the **ComplexHeatmap** package!  
🔗 [https://github.com/jokergoo/ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap)

---

**Last updated:** *February 2025*
