# Thermocline Reconstruction

MATLAB code for the figures and tables in:

> Lakhani, K. Q., Lynch-Stieglitz, J., & Findley, B. (2024). Reconstructing the
> tropical thermocline from oxygen-isotopes in planktonic and benthic
> foraminifera. *Paleoceanography and Paleoclimatology*, 39, e2023PA004794.
> https://doi.org/10.1029/2023PA004794

The paper describes a regression method that reconstructs the upper 600 m of
the tropical water column δ¹⁸O profile from the δ¹⁸O of foraminifera that
calcify at different depths:

- *G. ruber albus* sets the mixed layer value.
- *G. tumida*, *N. dutertrei*, and *P. obliquiloculata* constrain the
  subsurface.
- Benthic *U. peregrina* from core KNR195-5 GGC43 (617 m) anchors the
  bottom of the profile.

Thermocline depth is read from each reconstructed profile as the depth where
δ¹⁸O = −0.66‰, the value that matches the 20 °C isotherm (Z20) in today's
tropical Pacific. The method is tested at 20 Holocene core sites across the
tropical Pacific and then applied to the Last Glacial Maximum (LGM) at ODP 849.

## Requirements

- MATLAB. The scripts use only core MATLAB functions plus `prctile` and
  `nanmean`/`nanstd` from the Statistics and Machine Learning Toolbox.
- The data files listed under [Data](#data).

## Repository layout

```
Thermocline_Reconstruction/
├── Figure_*.m, Supp_Figure_*.m, Table_2.m   One script per figure or table in the paper
├── Functions/                               The model, helper functions, and small data files
└── Data/                                    Input data (not included, see below)
```

### Figures and tables

| Script | Paper item |
| --- | --- |
| `Figure_1_review.m` | **Figure 1.** (a–d) Holocene δ¹⁸O for each planktonic species at the core sites, over the climatological δ¹⁸O at that species' apparent calcification depth (ACD). (e) Difference in predicted δ¹⁸O at 600 m relative to the benthic core site, with the benthic core marked by a star. |
| `Figure_1_old.m` | Earlier version of Figure 1 from before peer review, without climatology behind panels a–d. Kept for reference. |
| `Figure_2.m` | **Figure 2.** Model profiles at VM28-234 (western), ML1208-19GC (central), and ODP 849 (eastern Pacific), without (a–c) and with (d–f) benthic data, including Monte Carlo spread from ACD uncertainty. (g) Equatorial transect of climatological δ¹⁸O with regional mean regression profiles. |
| `Figure_3.m` | **Figure 3.** (a) Holocene Z20 from the regression at each core, over the climatological Z20-equivalent depth. (b) Equatorial (5°S–5°N) east–west transect of Z20 with Monte Carlo error bars. |
| `Figure_4.m` | **Figure 4.** Reconstructed versus climatological Z20 at every core, with 95% Monte Carlo intervals, linear fit, and correlation (R = 0.97). |
| `Figure_5.m` | **Figure 5.** Holocene (a) and LGM (b) profiles at ODP 849, and both together (c) after removing the glacial ice volume (1.0‰) and cooling (0.525‰) signals. |
| `Figure_6.m` | **Figure 6.** *G. ruber* minus each subsurface species, and minus the subsurface mean, against climatological Z20. |
| `Figure_7.m` | **Figure 7.** Each subsurface species, and the subsurface mean, minus benthic δ¹⁸O, against climatological Z20. |
| `Table_2.m` | **Table 2.** Climatological and reconstructed Z20 for each core, with Monte Carlo standard deviation and 2.5/97.5 percentiles. **Run `Figure_4.m` first**: this script reads variables from that script's workspace. |
| `Supp_Figure_1.m` | **Figure S1.** Sensitivity of the profile at ML1208-19GC to each model parameter (β₁, β₂, MLD) and to ACD uncertainty. |
| `Supp_Figure_2_3_4_5.m` | **Figures S2–S5.** Model profile at every core, grouped by region. Opens one figure per core and pauses after each region; press any key in the Command Window to continue. |
| `Supp_Figure_6.m` | **Figure S6.** Depth where climatological δ¹⁸O = −0.66‰ against the depth of the 20 °C isotherm, showing that this δ¹⁸O value tracks Z20. |

### Functions

| File | Purpose |
| --- | --- |
| `run_thermocline_model3.m` | The regression model (paper Section 2.2, Equation 1). The profile is constant at the *G. ruber* value through the mixed layer and follows `-β₁^-(z+β₂) + β₃` below it. β₃ is fixed by continuity at the base of the mixed layer. β₁, β₂, and the mixed layer depth are found by minimising RMS misfit to the data, first on a coarse grid and then on a finer grid around the best coarse value. Returns `[beta, MLD]`. |
| `depth_realization.m` | Draws a random calcification depth from a species' ACD distribution, within set bounds. Used for the Monte Carlo uncertainty estimates. |
| `Depth_Distributions.mat` | ACD distributions for each species from Lakhani et al. (2022), stored as `xi_tot`, `ksdens`, and `minmax`. Used by `depth_realization`. The data behind these distributions is in [lakhani1118/Foraminifera-Calcification](https://github.com/lakhani1118/Foraminifera-Calcification). |
| `cycle_data_3D.m` | Shifts a gridded field from −180–180° longitude to 30–390° so the Pacific is continuous across the date line. |
| `world_coastp.m`, `rgrd_WorldCstLinePacific.mat` | Draws Pacific-centred coastlines on map figures. |
| `colorblind_jet.m` | Colourblind-friendly diverging colormap ([ColorBrewer](https://colorbrewer2.org/#type=diverging&scheme=RdYlBu&n=9) RdYlBu). |
| `mindata.m`, `maxdata.m` | Clamp an array index to the valid range. |

## Data

The scripts load their inputs from a `Data/` folder in the repository root.
That folder is not in this repository.

| File | Contents | Source |
| --- | --- | --- |
| `Z20.mat` | `good_data`: one row per core, with columns latitude, longitude, and Holocene (0–6 ka) δ¹⁸O for *G. ruber*, *T. sacculifer*, *G. tumida*, *N. dutertrei*, and *P. obliquiloculata*. `Z20`: climatological 20 °C isotherm depth at each core. | Paper Table 1. The individual measurements are archived at NOAA NCEI: https://doi.org/10.25921/eywq-n941 |
| `d18O_paper.mat` | Gridded predicted δ¹⁸O of calcite (`d18O1` on grids `X3`, `Y3`, `Z3`), as in Lakhani et al. (2022). | World Ocean Atlas 2013 temperature (Locarnini et al., 2013), LeGrande & Schmidt (2006) seawater δ¹⁸O, and the Kim & O'Neil (1997) paleotemperature equation |
| `WOA_04_temp.mat` | World Ocean Atlas 0.25° annual mean temperature. Used by `Figure_4.m`. | World Ocean Atlas |
| `TP_80.mat` | Depth where temperature has fallen 80% of the way from the surface to the deep value (`TP_data`). Used by `Figure_4.m`. | Calculated from World Ocean Atlas temperature |
| `Legrande_Schmidt.mat` | Gridded seawater δ¹⁸O. Used by `Figure_4.m`. | LeGrande & Schmidt (2006) |

Constants written into the scripts:

| Value | Meaning |
| --- | --- |
| ACDs of 17–18 m (*G. ruber*), 210 m (*G. tumida*), 114 m (*N. dutertrei*), and 94–96 m (*P. obliquiloculata*) | Average apparent calcification depths from Lakhani et al. (2022) |
| Benthic δ¹⁸O of 1.75‰ | Holocene *U. peregrina* at KNR195-5 GGC43 after the 0.47‰ species offset |
| Benthic δ¹⁸O of 2.797‰ at about 500 m | LGM benthic value, with the depth moved up 121 m for lower glacial sea level |
| δ¹⁸O of −0.66‰ | Z20-equivalent δ¹⁸O value |

## Running the scripts

1. Put the data files in `Data/`.
2. In MATLAB, set the current folder to the repository root. Each script adds
   `Functions/` and `Data/` to the path itself.
3. Run a script, for example `Figure_4`.

Notes:

- Cores VM28-227 and VM28-229 are left out because their δ¹⁸O stratigraphy is
  uncertain. That leaves the 20 cores in the paper.
- *T. sacculifer* δ¹⁸O is in `good_data` but is not used by the model.
- The Monte Carlo steps draw random ACDs, so error bars and profile spreads
  vary a little from run to run and from the published numbers. To get the
  same output every time, call `rng(seed)` before running a script. The
  scripts that run 100–500 realizations per core (`Figure_3.m`, `Figure_4.m`,
  and `Figure_5.m`) can take several minutes.
- The version of the code used for the paper is release
  [`tag1`](https://github.com/lakhani1118/Thermocline_Reconstruction/releases/tag/tag1).

## Citation

If you use this code, please cite the paper above and the archived code
release:

> Lakhani, K. Q. (2024). lakhani1118/Thermocline_Reconstruction: Thermocline
> reconstruction v1.0 (version tag1) [Software]. Zenodo.
> https://doi.org/10.5281/zenodo.10806892

The foraminiferal δ¹⁸O data:

> Lakhani, K. Q., Lynch-Stieglitz, J., & Findley, B. (2024). NOAA/WDS
> Paleoclimatology – Tropical Pacific planktonic foraminiferal d18O records for
> 0–6 ka [Dataset]. NOAA National Centers for Environmental Information.
> https://doi.org/10.25921/eywq-n941

The ACD distributions:

> Lakhani, K. Q., Lynch-Stieglitz, J., & Monteagudo, M. M. (2022). Constraining
> calcification habitat using oxygen isotope measurements in tropical planktonic
> foraminiferal tests from surface sediments. *Marine Micropaleontology*, 170,
> 102074. https://doi.org/10.1016/j.marmicro.2021.102074

The data from that paper is at
[lakhani1118/Foraminifera-Calcification](https://github.com/lakhani1118/Foraminifera-Calcification).

## Contact

For questions, open an issue on this repository.
