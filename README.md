# zero_offset.py 使用说明（XRD 零点补正 + 晶格常数拟合）

`zero_offset.py` 用于把 **用户给定晶系（数字 1–8）** 与 **若干已指数化的峰位 (h k l, 2θ)** 作为输入，完成以下工作：

- 读取 **2θ 原点**（`2THETA ORIGIN = δ`，单位：**度**）并按以下约定显示/计算：  
  - 显示用的 `2THETA_OBS = 原始 OBS − δ`  
  - 计算用的 `2THETA_CAL` **仅由晶格常数**计算（不再减 δ）  
  - 表中 **“CAL-OBS”** 一列打印 **`OBS_shown − CAL`**
- 按所选晶系在 **Q 空间**做线性回归，拟合晶格常数（立方/四方/六方/正交/单斜/三斜均支持）
- 提供两种工作模式：  
  - `fixed`：固定 δ，拟合晶格参数  
  - `scan`：在区间内**网格扫描** δ（如 −0.5…0.5，步长 0.1），按目标函数选出最佳 δ
- 输出 **urlap 风格报表**（标题、晶系、波长与 δ、峰表、Direct Cell 常数）

---

## 支持的晶系（在输入文件第 2 行首个数字指定）

| 代码 | 晶系 | 备注 |
|---:|---|---|
| 1 | triclinic（斜方晶） | a≠b≠c, α≠β≠γ |
| 2 | monoclinic(b) | 唯一直轴 b，β≠90° |
| 3 | monoclinic(c) | 唯一直轴 c，γ≠90° |
| 4 | orthorhombic | a≠b≠c，α=β=γ=90° |
| 5 | tetragonal | a=b≠c |
| 6 | cubic | a=b=c |
| 7 | trigonal | 以六方轴表示 |
| 8 | hexagonal | 六方 |

> 注：7（三方）与 8（六方）在 d 间距公式上同用六方轴设定：
> \( Q = \frac{4}{3}\frac{h^2+hk+k^2}{a^2} + \frac{l^2}{c^2} \)。

---

## 输入文件格式（RAW / urlap 风格）

```text
TITLE 你的标题
<晶系代码> 0 0 0 0 0 0 0 0
<波长(Å)>   <2THETA_ORIGIN(度)>
h k l  2theta1
h k l  2theta2
...
1000 0 0 0.0
```

- **第 1 行**：任意标题  
- **第 2 行**：首个数字为**晶系代码（1–8）**，其余给 0 即可  
- **第 3 行**：波长 λ（Å）与 `2THETA ORIGIN`（δ，单位**度**）  
- **第 4 行起**：每行一个峰：`h k l  2θ`（2θ 单位**度**）  
- **终止行**：`1000 0 0 0.0`

> ⚠️ **强烈注意**：δ 的单位是**度（deg）**，不是弧度（rad）。

---

## 安装与环境

```bash
python -V      # 建议 Python 3.8+
pip install numpy
```

将 `zero_offset.py` 放到你的项目或工作目录（也可以加入 CI 任务）。

---

## 命令行用法

### 1) 固定 δ，拟合并输出报表

```bash
python zero_offset.py \
  -i INPUT.TXT \
  --mode fixed \
  --out OUTPUT.txt
```

- 默认从输入第 3 行读取 λ 与 δ（可用 `--wavelength` / `--delta` 覆盖）。
- 输出为 urlap 风格 TXT：  
  - `2THETA_ORIGIN = δ`  
  - `2THETA_OBS = OBS_raw − δ`（显示值）  
  - `2THETA_CAL` 由拟合的晶格常数计算  
  - `CAL-OBS = OBS_shown − CAL`  
  - 末尾包含 `DIRECT CELL CONSTANT`（a,b,c 与不确定度）

### 2) 扫描 δ（网格点，例如 −0.5…0.5，步长 0.1）

```bash
python zero_offset.py \
  -i INPUT.TXT \
  --mode scan \
  --delta-range -0.5 0.5 \
  --step 0.1 \
  --select-by uncprod \
  --out OUTPUT_scan.txt
```

- `--select-by`（选择最优 δ 的目标函数）：  
  - `rms`：最小化角度域 RMS（`RMS(2THETA_CAL − 2THETA_OBS_shown)`）  
  - `uncprod`：**不确定度乘积最小**：  
    - 立方：`da³`  
    - 四方 / 六方 / 三方：`da² · dc`  
    - 其他（直方/单斜/三斜）：`da · db · dc`

该模式会在每个 δ 上拟合参数并计算目标函数，然后选**最优 δ**写入报表。

---

## 示例

### 例 1：四方晶（代码 5），固定 δ = −0.1°

**输入 `TETRA_INPUT.TXT`**：
```text
TITLE Fe(Se0.25Te0.75)0.92 10K
    5   0    0    0    0    0    0    0    0
   1.5402       -0.1
    2    0    0   30.535
    1    0    1   31.922
    2    2    0   43.788
    2    1    1   44.815
 1000    0    0    0.000000
```

**运行**：
```bash
python urlap_report_v3.py -i TETRA_INPUT.TXT --mode fixed --out TETRA_out.txt
```

### 例 2：同一文件，扫描 δ（−0.5…0.5，步长 0.1），以不确定度乘积最小为准

```bash
python urlap_report_v3.py \
  -i TETRA_INPUT.TXT \
  --mode scan --delta-range -0.5 0.5 --step 0.1 \
  --select-by uncprod \
  --out TETRA_scan_uncprod0p1.txt
```

---

## 输出字段说明

- **TITLE**：原样打印  
- **LATTICE SYSTEM**：按代码识别输出（`TRICLINIC / MONOCLINIC(b) / MONOCLINIC(c) / ORTHORHOMBIC / TETRAGONAL / CUBIC / TRIGONAL / HEXAGONAL`）  
- **WAVE LENGTH**：波长 λ（小数 5 位）  
- **2THETA ORIGIN**：选定的 δ（**度**）  
- **峰表（每峰一行）**：  
  - `H K L`：米勒指数  
  - `2THETA_OBS`：`OBS_raw − δ`  
  - `2THETA_CAL`：由晶格常数计算的布拉格角  
  - `CAL-OBS`：`OBS_shown − CAL`（注意与原 urlap 表头一致，但符号按此约定输出）  
- **DIRECT CELL CONSTANT**：  
  - `A, B, C` 与其不确定度 `DA, DB, DC`（由 Q 空间回归协方差通过**数值雅可比**传播得到）  
- **ALPHA/BETA/GAMMA**：按晶系输出（正交/四方/立方均为 90°；六方 γ=120°；单斜/三斜按拟合值）

> 如需同时输出 **Reciprocal 常数块**（A\*, DA\* 等），可以在后续版本追加（当前 v3 为简洁起见只输出 Direct 块）。

---

## 计算方法（简述）

- 使用 \(Q = (2\sin\theta/\lambda)^2 = 1/d^2\) 在 **Q 空间**做**过原点线性回归**。  
- 不同晶系对 \(Q\) 与 \(h,k,l\) 的关系均为**线性**：  
  - 立方：\(Q=(h^2+k^2+l^2)/a^2\)  
  - 四方：\(Q=(h^2+k^2)/a^2 + l^2/c^2\)  
  - 六方/三方：\(Q=\frac{4}{3}\frac{h^2+hk+k^2}{a^2} + \frac{l^2}{c^2}\)  
  - 正交：\(Q=h^2/a^2 + k^2/b^2 + l^2/c^2\)  
  - 单斜(b/c) 与 三斜：先线性回归参数，再还原 \(a,b,c,(\alpha,\beta,\gamma)\)  
- 协方差：\(\mathrm{cov}(m) = \sigma^2 (X^\top X)^{-1}\)，\(\sigma^2=\mathrm{RSS}/(n-p)\)  
- 误差传播：通过**数值雅可比** \(\mathbf{J}=\partial(a,b,c)/\partial m\)，近似  
  \(\mathrm{cov}(a,b,c) \approx \mathbf{J}\,\mathrm{cov}(m)\,\mathbf{J}^\top\)。

---

## 常见问题（FAQ / 排错）

**Q1. δ 的单位是什么？**  
A. 始终是 **度（deg）**。不要用弧度。若 UI 有其它单位，请先换算成度填入第 3 行。

**Q2. 为什么 “CAL-OBS” 列有时为负？**  
A. 我们按你的习惯输出 **`OBS_shown − CAL`**。正负号仅表示 `CAL` 与 `OBS_shown` 的相对偏差方向。

**Q3. 峰数不够会怎样？**  
A. 需要的最少峰数 ≥ 模型参数个数（例如：立方≥1，四方/六方≥2，正交≥3，单斜≥4，三斜≥6）；不足会导致病态或高不确定度。

**Q4. 扫描 δ 的推荐步长？**  
A. UI 里常用 **0.1°** 粗扫；若要更准，可在最佳点附近做细扫（v3 以简单为主，如需可加“细化”选项）。

**Q5. 不确定度为什么与其它软件略有差异？**  
A. 加权方式、残差模型、是否强制过原点等细节不同都会造成差异。v3 采用**过原点线性回归**与**数值雅可比传播**。

---

## 选项速查

| 选项 | 说明 |
|---|---|
| `-i, --input` | 输入文件路径（必填） |
| `--mode {fixed,scan}` | 运行模式：固定 δ 或 扫描 δ |
| `--wavelength <Å>` | 覆盖输入文件中的波长 λ |
| `--delta <deg>` | 覆盖输入文件中的 δ（度；仅 `fixed` 模式用） |
| `--delta-range <min> <max>` | 扫描 δ 的范围（度；默认 `-0.5 0.5`） |
| `--step <deg>` | 扫描步长（度；默认 0.1） |
| `--select-by {rms,uncprod}` | 选择最优 δ 的目标函数（默认 `uncprod`） |
| `--out <path>` | 输出报表路径（默认 `<input>.out.txt`） |



---

## 许可与贡献


# urlap_autopipe.py

一键把 **urlap 晶格拟合（自动寻找零点 δ）** 和 **XRD 光谱“零点补正 → 出图 → 导出 CSV”** 串起来。  
无需手动输入 δ：脚本先调用 `urlap_report_v3.py` 生成报表并解析 `2THETA ORIGIN = ...`，再把该 δ 传给绘图脚本 `plot_xrd_before_after_split.py`。

## 功能
- 两种模式：`scan`（网格扫描 δ，自动选最佳）/ `fixed`（使用输入或命令行 δ）
- 零点补正约定：`2θ_after = 2θ_before − δ`（°）
- 一次性输出：叠加图、分开图（Before/After）、差分（After−Before）、比值（After/Before）+ 对应 CSV

## 依赖与放置
```bash
pip install numpy pandas matplotlib
```
将三个脚本放在同一目录：`urlap_autopipe.py`, `urlap_report_v3.py`, `plot_xrd_before_after_split.py`

## 快速开始
```bash
python urlap_autopipe.py   --peaks LATTICE9.TXT   --spectrum before.csv   --mode scan --delta-range -0.5 0.5 --step 0.1 --select-by uncprod   --xlim 5 90 --compare both --separate --save-csv   --out-prefix sample1
```
固定 δ（若 RAW 第 3 行已有 δ 可以省略 `--delta`）：
```bash
python urlap_autopipe.py   --peaks TETRA_INPUT.TXT   --spectrum before.csv   --mode fixed   --xlim 5 90 --separate --compare both --save-csv
```

## 输出
以 `--out-prefix sample1` 为例：生成 `sample1_urlap.out.txt`，以及 `sample1_before_after.png / _before.png / _after.png / _compare_diff.png / _compare_ratio.png`，和 `sample1_before.csv / _after.csv / _compare.csv`。

## δ 的来源逻辑
1) `--delta` 2) `--delta-report <REPORT.txt>` 3) CSV 列 `delta_deg_used/delta/origin/twotheta_origin` 4) 同名报告 5) `--delta-default`（0.0）
（在管线中自动传入 `--delta-report`，无需手工输入 δ）

## 许可证
MIT（或按项目实际许可证）

欢迎提交 Issue / PR 以完善更多晶系细节、加权策略、细化扫描、**Reciprocal 常数块**输出、补正谱/UMAP 数据导出等功能。

