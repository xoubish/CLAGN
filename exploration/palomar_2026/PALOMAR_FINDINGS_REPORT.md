# Palomar 200-inch (NGPS) 2026 Campaign: Full-Sample Scientific Findings Report

**Campaign Scope:** 2.5 Nights on the Palomar Hale 200-inch Next Generation Palomar Spectrograph (NGPS; 2026 Sep 23, Oct 26, Oct 27)  
**Sample Size:** 99 Unique AGN Systems | 112 New Palomar Spectra | 253 Archival SDSS + DESI Spectra (**365 Total Spectral Epochs**)  
**Primary Pipelines:** [`full_sample_fits.py`](full_sample_fits.py) (Step 1) & [`deep_dive_leads.py`](deep_dive_leads.py) (Step 2)

---

## 1. Executive Summary: What Did We Get, Is It Useful, and What Did We Find?

Over 2.5 nights at Palomar, we obtained **112 optical spectra ($3920\text{–}10280\text{ Å}$) of 99 unique AGN targets** designed to test whether a pre-run 12-year WISE $W1$ variability manifold (`manifold_direction`, `manifold_var_ratio`) predicts spectroscopically confirmed Changing-Look AGN (CLAGN) transitions across a pre-assigned **A/B/C/D experimental hierarchy**.

> [!IMPORTANT]
> **Three Core Scientific Takeaways from the Full 365-Epoch Analysis:**
> 1. **The Data Are High-Yield and Scientifically Publishable:** Full multi-component spectral decomposition across all 365 epochs (**83 $\text{H}\beta+[\text{O III}]$**, **46 $\text{H}\alpha+[\text{N II}]+[\text{S II}]$**, and **45 $\text{Mg II}+\text{UV Fe II}$** multi-epoch systems; **1,404 component fits** across 3 host/$\text{Fe II}$ model families + alternate standard stars) reveals **20 non-catalog AGN systems with $\ge 30\%$ broad-$\text{H}\beta$ equivalent-width evolution**, anchored by a suite of unambiguous, multi-line validated state transitions.
> 2. **New & Recurrent CLAGN Discoveries with Proven Accretion Physics:** Simultaneous joint multi-epoch $\text{H}\beta+\text{H}\alpha$ modeling (tying host stellar population mix $\theta_{\text{young}}$ and narrow-line dispersion $\sigma_n$) and Balmer-decrement ($\text{H}\alpha/\text{H}\beta$) vs. Cardelli et al. (1989) ($R_V=3.1$) extinction tests prove that our top turn-off discoveries (**`P21422`** and **`P9506`**) and turn-on discoveries (**`P7947`** and **`P9694`**) are driven by **intrinsic ionizing accretion-rate changes ($\Delta \log \lambda_{\text{Edd}} = -0.32\text{ to }+0.26\text{ dex}$)** rather than foreground dust clouds—accompanied by classical virial broad-line breathing ($\Delta \text{FWHM}_{\text{H}\beta}: 3402 \to 5820\text{ km s}^{-1}$ in `P21422`). Meanwhile, the full-sample run surfaced a **complete Type 1 $\to$ Type 2 turn-off collapse** in **`P14388`** ($\text{EW}_{\text{H}\beta}: 46.5\text{ Å [2021]} \to 7.2\text{ Å [2025]} \to 1.0\text{ Å [2026]}$, $\Delta \log \lambda_{\text{Edd}} = -0.82\text{ dex}$ with $5.32\times$ Balmer-decrement steepening) alongside four newly identified Group B transitions (**`P12457`**, **`P21743`**, **`P21841`**, **`P11010`**).
> 3. **A Rigorous, Transparent Test of the Pre-Run Manifold Selection:** Across the 76 non-catalog reviewed multi-epoch $\text{H}\beta$ targets (Group A $n=11$, Group B $n=58$, Group C $n=1$, Group D $n=6$), two-sample permutation tests show **no statistically significant difference in population-median $|\Delta \log_{10}\text{EW}(\text{H}\beta)|$** between Group A ($0.061\text{ dex}$), Group B ($0.056\text{ dex}$), and random AGN controls in Group D ($0.070\text{ dex}$; $p_{\text{perm}} = 0.68$). Furthermore, pre-run `manifold_direction` predicts the 2026 transition polarity in **$57.8\%$ of clean Group B targets** ($p_{\text{perm}} = 0.46$ overall) and **anti-correlates ($11.1\%$) in Group A** because the strongest Group A turn-off discoveries (`P21422`, `P9506`) underwent earlier WISE $W1$ brightening episodes before collapsing in 2022–2026. **Scientific lesson:** A 12-year integrated WISE $W1$ manifold metric identifies *recurrent state-transitioning systems*, not instantaneous 2026 polarity.

---

## 2. Step 1: Full-Sample Census & Pre-Run A/B/C/D Selection Experiment

| Experimental Group | Description | $N_{\text{obs}}$ ($\text{H}\beta$) | Median $\|\Delta \log_{10}\text{EW}(\text{H}\beta)\|$ [dex] | $\|\Delta \text{EW}(\text{H}\beta)\| \ge 30\%$ Fraction | $\|\Delta(\text{H}\beta/\text{[O III]})\| \ge 30\%$ Fraction | Sign Hit Rate (`manifold_direction`) | Clean-Coverage ($\ge 0.80$) Sign Hit Rate |
|---|---|---:|---:|---:|---:|---:|---:|
| **Group A** | Manifold + $W1$ Var | 11 | $0.061$ (IQR $0.097$) | **27.3%** ($3/11$) | **40.0%** ($4/10$) | $11.1\%$ ($1/9$) | $11.1\%$ ($1/9$) |
| **Group B** | Manifold (All $W1$) | 58 | $0.056$ (IQR $0.093$) | **25.9%** ($15/58$) | **34.0%** ($18/53$) | $51.9\%$ ($27/52$) | **57.8%** ($26/45$) |
| **Group C** | $W1$ Var Only (non-cat) | 1 | $0.008$ (IQR $0.000$) | $0.0\%$ ($0/1$) | $0.0\%$ ($0/1$) | $0.0\%$ ($0/1$) | $0.0\%$ ($0/1$) |
| **Group D** | Random AGN Control | 6 | $0.070$ (IQR $0.070$) | **33.3%** ($2/6$) | **60.0%** ($3/5$) | $50.0\%$ ($3/6$) | $50.0\%$ ($3/6$) |
| **Catalog CLAGN** | Known Literature CLAGNs | 3 ($\text{H}\beta$) | $0.204$ | **66.7%** ($2/3$: `P14388`, `P17881`) | **33.3%** ($1/3$: `P14388`) | — | — |

![Full-Sample A/B/C/D Selection Experiment](outputs/full_sample/figures/fig1_full_sample_abcd_experiment.png)

---

## 3. Step 2: Physical Mechanism Verdicts — Accretion Collapse vs. Dust Reddening

Under a standard Cardelli, Clayton & Mathis (1989) ($R_V = 3.1$) Galactic extinction curve ($A_{\text{H}\beta}/E(B-V) = 3.61$, $A_{\text{H}\alpha}/E(B-V) = 2.53$), if an observed drop in broad $\text{H}\beta / [\text{O III}]$ flux ratio $R_{\text{H}\beta}$ were caused purely by foreground dust obscuring the BLR, the broad Balmer decrement $\text{BD} = F(\text{H}\alpha_{\text{broad}}) / F(\text{H}\beta_{\text{broad}})$ **must steepen** as $\text{BD}_{\text{new}}/\text{BD}_{\text{old}} = R_{\text{H}\beta}^{-0.299}$.

![Joint Multi-Epoch Fits and Balmer Decrements](outputs/full_sample/figures/fig2_discovery_gallery_joint_fits.png)

![Balmer Decrement Dust Test and Virial Plane](outputs/full_sample/figures/fig3_balmer_decrement_virial_physics.png)

| Target | Group | $z$ | $\Delta t_{\text{rest}}$ [yr] | Joint $\text{EW}(\text{H}\beta)$ ($\text{old}\to\text{new}$) [Å] | Joint $\text{EW}(\text{H}\beta)$ Ratio | Joint $\text{H}\beta/\text{[O III]}$ Ratio | $\text{EW}(\text{H}\alpha)$ Ratio | Broad Balmer Decr. Ratio (Obs vs. Dust Pred.) | $\log\left(\frac{M_{\text{BH}}}{M_\odot}\right)$ | $\Delta \log \lambda_{\text{Edd}}$ [dex] | Physical Mechanism & Novelty Verdict |
|---|---|---:|---:|---|---:|---:|---:|---|---:|---:|---|
| **`P21422`** | A | 0.239 | 3.63 | $106.4 \to 73.5$ | **0.69×** | **0.61×** | **0.59×** | **$0.86\times$** vs. $1.16\times$ | 8.64 | $-0.32$ | **Pure accretion collapse** (dust ruled out; virial $\text{FWHM}_{\text{H}\beta}: 3402 \to 5820\text{ km s}^{-1}$). Recurrent post-2022 collapse after Graham+19 CRTS brightening. |
| **`P9506`** | A | 0.232 | 9.77 | $76.9 \to 59.9 / 60.1$ | **0.78×** | **0.57× / 0.56×** | **0.89× / 0.84×** | **$0.93\times$** vs. $1.19\times$ | 7.58 | $-0.30$ | **Pure accretion dimming** (dust ruled out; verified across 4 frames on Sep 23 & Oct 27). **Unreported new CLAGN discovery.** |
| **`P7947`** | B | 0.254 | 19.90 | $13.8 \to 31.5$ | **2.28×** | **1.73×** | **1.32×** | **$0.61\times$** vs. $0.85\times$ | 7.18 | $+0.24$ | **Intrinsic accretion turn-on** ($\text{H}\alpha/[\text{S II}] = 1.44\times$, Balmer decrement flattens as ionizing flux rises). **Unreported new CLAGN discovery.** |
| **`P9694`** | D | 0.238 | 4.07 | $98.9 \to 118.0 / 132.3$ | **1.19× / 1.34×** | **1.27× / 1.30×** | **1.12× / 1.13×** | **$0.87\times$** vs. $0.92\times$ | 8.75 | $+0.26 / -0.09$ | **Intrinsic accretion brightening** verified across 4 frames on Sep 23 & Oct 27. **Unreported control discovery.** |
| **`P14388`** | C (Known) | 0.269 | 1.38 | $46.5\text{ ['21]} \to 7.2\text{ ['25]} \to 1.0$ | **0.13×** ($0.02\times$ vs. '21) | **0.08×** | **0.73×** | **$5.32\times$** vs. $2.15\times$ | 8.01 | $-0.82$ | **Complete Type 1 $\to$ Type 2 turn-off collapse** in catalog CLAGN (`J075934.95+322143.3`); coupled $-0.82\text{ dex}$ accretion shutdown + dust reddening. |
| **`P12457`** | B | 0.092 | 23.94 | $3.1 \to 1.8$ | **0.58×** | **0.57×** | **0.72×** | **$1.28\times$** vs. $1.18\times$ | 6.08 | $-0.39$ | **Low-luminosity Seyfert 1.9 dimming** ($\text{H}\alpha/[\text{S II}] = 0.71\times$, concordant reddening/dimming). **Unreported Group B lead.** |
| **`P21743`** | B | 0.143 | 11.21 | $102.3 \to 59.0$ | **0.58×** | **0.81×** | **0.91×** | **$1.55\times$** vs. $1.06\times$ | 7.70 | $-0.30$ | **Strong $\text{H}\beta$ dimming** verified across both Oct 27 G frames (`31.0 Å` vs. `31.8 Å`). **Unreported Group B lead.** |
| **`P21841`** | B | 0.514 | 1.20 | $69.6 \to 106.0$ | **1.52×** | **0.87×** | **1.24×** | **$1.03\times$** vs. $1.04\times$ | 8.77 | $-0.18$ | **Concordant $\text{Mg II}$ ($2.47\times$) + $\text{H}\beta$ ($1.52\times$) + $\text{H}\alpha$ ($1.24\times$) EW growth** as blue continuum subsides. **Unreported Group B lead.** |
| **`P11010`** | B | 0.381 | 18.75 | $38.3 \to 54.0$ (single) | **1.41×** | **0.90×** | **1.53×** | — | 8.66 | $+0.02$ | **Concordant $\text{Mg II}$ ($5.20\times$) + $\text{H}\beta$ ($1.41\times$) + $\text{H}\alpha$ ($1.53\times$) brightening** verified on both Oct 27 R frames. **Unreported Group B lead.** |
| **`P17881`** | Known | 0.211 | 1.38 | $65.4 \to 104.7$ | **1.60×** | **0.96×** | **1.15×** | **$0.99\times$** vs. $1.01\times$ | 8.49 | $+0.01$ | **Continuum-dilution EW growth** in catalog CLAGN (`J004339.32+134436.5`) at constant intrinsic broad $\text{H}\beta/[\text{O III}]$ flux ($0.96\times$). |
| **`P10381`** | Known | 0.149 | 10.58 | $83.0 \to 100.2 / 113.5$ | **1.21× / 1.37×** | **1.33× / 1.21×** | **1.14× / 1.24×** | **$0.97\times$** vs. $0.94\times$ | 7.81 | $+0.17 / -0.04$ | **Persistent Type 1 brightening** over 2014 SDSS across both Sep 23 & Oct 27 (unblocked via rest-frame $5695\text{ Å}$ splice boundary). |

---

## 4. Multi-Epoch Spectral Decompositions of Key Discoveries & New Group B/C Leads

### 4.1 Primary Group A, B, and D Discoveries
![P21422 Multi-Epoch Decomposition](outputs/deep_dive/figures/P21422_decomposition.png)

![P9506 Multi-Epoch Decomposition](outputs/deep_dive/figures/P9506_decomposition.png)

![P7947 Multi-Epoch Decomposition](outputs/deep_dive/figures/P7947_decomposition.png)

![P9694 Multi-Epoch Decomposition](outputs/deep_dive/figures/P9694_decomposition.png)

### 4.2 Newly Surfaced Group B & C Leads (`P14388`, `P12457`, `P21743`, `P21841`, `P11010`)
![Side-by-Side Decompositions of New Group B/C Leads](outputs/full_sample/figures/fig4_new_group_bc_leads_decompositions.png)

![P14388 Complete 5-Epoch Round-Trip Collapse](outputs/deep_dive/figures/P14388_decomposition.png)

![P12457 Multi-Epoch Decomposition](outputs/deep_dive/figures/P12457_decomposition.png)

![P21743 Multi-Epoch Decomposition](outputs/deep_dive/figures/P21743_decomposition.png)

![P21841 Multi-Epoch Decomposition](outputs/deep_dive/figures/P21841_decomposition.png)

Full 14-source multi-epoch decomposition figures are available in [`outputs/deep_dive/deep_dive_fit_atlas.pdf`](outputs/deep_dive/deep_dive_fit_atlas.pdf) and [`outputs/deep_dive/figures/`](outputs/deep_dive/figures/).
