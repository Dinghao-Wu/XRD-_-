[xrd_to_umap.md](https://github.com/user-attachments/files/22290300/xrd_to_umap.md)

# xrd_to_dataset.py — Build an ML‑ready XRD dataset (numeric‑only CSV)

This utility turns a folder of local XRD profiles into a **machine‑learning ready matrix** by
reading two‑column files (2θ, intensity), **resampling to a common 2θ grid**, and **L2 normalizing** each pattern.
It writes a **numeric‑only** CSV for features and a separate CSV that maps each row back to the original file.

---

## Highlights

- **No external APIs** — works entirely on your local files.
- **Robust text parsing** — accepts `.txt/.csv/.xy/.xye/.dat` with space/tab/comma/semicolon separators and skips header/comment lines (`# ; ! * //`).
- **Vendor formats (optional)** — supports `.ras/.uxd/.raw/.cpi/...` via [`xylib-py`](https://pypi.org/project/xylib-py/) if installed.
- **XRDML** support (`.xrdml`) via a lightweight XML reader.
- **Common 2θ grid** — default **10–80°** with **0.02°** step (configurable).
- **Per‑sample L2 normalization** — emphasizes **shape** over absolute intensity.
- **Clean output** — feature matrix is **numeric‑only** (no headers, no file column), with a separate `targets.csv` mapping.

---

## Installation

```bash
pip install numpy pandas
# optional (to read .ras/.uxd/.raw/.cpi/...)  
pip install xylib-py
```

> `xylib-py` is a binary package. If you don’t install it, the script will still work for text/XRDML files and will **skip** vendor‑specific binaries.

---

## Command‑line usage

```bash
python xrd_to_dataset.py ./mydata --xmin 10 --xmax 80 --step 0.02
```

**Arguments**

- `input_dir` (positional): folder containing your XRD files.
- `--xmin` (float, default **10.0**): lower bound of the common 2θ grid (degrees).  
- `--xmax` (float, default **80.0**): upper bound of the common 2θ grid (degrees).  
- `--step` (float, default **0.02**): grid step (degrees).  

The number of columns in the output matrix equals:
\[
N_\text{points} = \left\lfloor \frac{\text{xmax}-\text{xmin}}{\text{step}} \right\rfloor + 1.
\]

---

## Input expectations

Each file should contain **two numeric columns** per line:

1. **2θ (deg)**
2. **Intensity**

The reader tolerates:
- Blank lines
- Comment/header lines beginning with `#`, `;`, `!`, `*`, or `//`
- Mixed separators: spaces, tabs, commas, or semicolons

Files are **sorted by x** and duplicate x values are removed before interpolation.

Supported types:
- **ASCII two‑column**: `.txt`, `.csv`, `.xy`, `.xye`, `.dat`
- **XRDML**: `.xrdml`
- **Vendor formats** (with `xylib-py`): `.ras`, `.uxd`, `.raw`, `.cpi`, `.udf`, `.xdd` (first data block is read)

---

## What the script does

1. **Discover files** in `input_dir` (common extensions; if none match, it tries all files in the folder).  
2. **Read & parse** each file  
   - ASCII → two‑column parser (auto separator, skip headers)  
   - XRDML → minimal XML reader (`<dataPoints>`)  
   - Vendor → `xylib` if available
3. **Sort & deduplicate** by 2θ.  
4. **Resample to a common 2θ grid** with linear interpolation (`out‑of‑range → 0`).  
5. **L2 normalize** each resampled vector:
   \[
   \mathbf{y}' = \frac{\mathbf{y}}{\max(\|\mathbf{y}\|_2,\, \varepsilon)}\quad (\varepsilon=10^{-12})
   \]
   Samples with (near)‑zero norm are **skipped**.
6. **Write outputs** (see below).

---

## Outputs

All files are written to the **current working directory**.

- **`calcd_patterns.csv`** — **numeric‑only** feature matrix (no header row, no file column, no index).  
  - Shape: `n_samples × n_points`  
  - Each row corresponds to one input file (after filtering/normalization).

- **`targets.csv`** — row mapping and labels:  
  ```text
  file,label
  sample_A.txt,0
  sample_B.csv,1
  ...
  ```
  - `label` equals the **row index** in `calcd_patterns.csv` (0‑based).  
  - Use this to map rows back to filenames.

**Console log** prints how many files were processed or skipped and the grid summary.

---

## Loading the outputs

### NumPy
```python
import numpy as np, pandas as pd
X = np.loadtxt("calcd_patterns.csv", delimiter=",")  # numeric-only matrix
targets = pd.read_csv("targets.csv")                 # file ↔ label mapping

# sanity check: L2 norms should be ~1
import numpy.linalg as LA
print(np.round(LA.norm(X, axis=1), 6)[:5])
```

### pandas
```python
import pandas as pd
X_df = pd.read_csv("calcd_patterns.csv", header=None)  # no header in file
targets = pd.read_csv("targets.csv")
```

---

## Tips & good practices

- **Choose a grid that covers all files.** If a file’s x‑range falls outside `[xmin, xmax]`, most of its data becomes zero and the sample may be skipped (`norm==0`).  
- **L2 normalization** makes samples comparable and is compatible with cosine/Euclidean distances used downstream.  
- This script **does not** perform baseline removal, smoothing, Kα2 stripping, or 2θ zero‑shift correction. If you need those, apply them **before** running this tool.  
- Vendor formats: install `xylib-py`. If not installed, such files will be skipped with a console note.

---

## Troubleshooting

- `No candidate XRD files.` → Check the folder path and file extensions.  
- `Grid error: --xmax must be > --xmin` → Fix your bounds; typical 2θ ranges are `5–80`, `10–90`, etc.  
- `norm==0` skip → The resampled vector is (near) all zeros, usually because the grid doesn’t overlap the file’s 2θ range.  
- `xylib not installed` → Install `xylib-py` to read `.ras/.uxd/.raw/.cpi/...`.

---

## How it works (internals)

- **File discovery**: looks for common XRD extensions; if none found, tries every file in the folder.  
- **Parsing order**: ASCII two‑column → XRDML → xylib (vendor).  
- **Interpolation**: `numpy.interp` with `left=0`, `right=0`.  
- **Grid size**: `n = floor((xmax - xmin)/step) + 1` points.  
- **Numeric output**: written via `numpy.savetxt` (comma‑separated), intentionally **without headers** for maximum downstream compatibility.

---

## License & attribution

Feel free to copy/modify the script within your project.  
If you publish results produced with this tool, please acknowledge that feature matrices were generated by *xrd_to_dataset.py (local XRD → common grid + L2 normalization)*.
