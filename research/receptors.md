# receptors.json: sources and assumptions

## Olfactory

### `glomerulus_map`
- Built from DoOR.data `door_mappings.csv` (https://raw.githubusercontent.com/ropensci/DoOR.data/master/data/door_mappings.csv) and checked against Task et al. 2022 *eLife* 11:e72599.
- `connectome_ORN_type` = `ORN_<glomerulus>`. This is the hemibrain/FlyWire convention.
- Receptors that DoOR marks as larval-only are flagged. Ignore them for the adult: Or85c, Or45a/b, Or30a, Or1a, Or22c, Or24a, Or59a, Or74a, Or94a/b.

**Corrections to the brief:**
- **Ir75a → DP1l**, not DL2. DL2d/DL2v belong to Ir75b/Ir75c.
- **Or83c → DC3**, not VA6. VA6 is Or82a.

**Renaming trap:** Schlegel 2021 and Task 2022 renamed three glomeruli:
- VC3l → VC3 (Or35a)
- VC3m → VC5 (Ir41a)
- old VC5 → VM6 (Rh50, ammonia)

Hemibrain v1.2 type strings may still use the old names; FlyWire uses the new ones.

### `door_consensus`
- Source: DoOR 2.0 (Münch & Galizia 2016 *Sci Rep* 6:21841), matrix `door_response_matrix.csv` from the same GitHub data folder.
- For each of 41 odorants I give two sets of values: `consensus_raw` (0–1, spontaneous rate included) and `consensus_sfr_reset` (value minus the SFR row, where 0 = spontaneous and negative = inhibition). The second follows DoOR's own `resetSFR` convention.
- Only receptors with data are listed. A missing receptor means untested, not zero.
- Ir columns come from calcium imaging of a handful of odorants (Silbering 2011). "1.0" there means best of the few tested.

**Top receptors (SFR-reset):**

| Odorant | Receptors (glomerulus, value) |
|---|---|
| Ethyl acetate | Or42a (VM7d) 0.72, Or42b (DM1) 0.61, Or59b (DM4) 0.58, Or85e (VC1), Or43b (VM2) |
| Isoamyl acetate | Or98a (VM5v) 0.64, Or22a (DM2) 0.48, Or10a (DL1), Or47a (DM3), Or43b (VM2) |
| Acetic acid | Ir64a-DC4 1.0, Ir75a (DP1l) 0.56, Ir64a-DP1m 0.48 |
| 2-Phenylethanol | Or67b (VA3) 0.83, Or67a (DM6), Or69a (D) |
| Isoamyl alcohol | Or35a (VC3) 0.67, Or67b (VA3) |
| Hexanal | Or85b (VM5d), Or35a (VC3), Or7a (DL5) |
| Acetoin | Or9a (VM3) 0.59, Or19a (DC1) |
| Eugenol | Or71a (VC2) 0.73, Or69a (D) |
| Diacetyl | Or92a (VA2) 0.69, Ir64a, Or42a |
| CO2 | Gr21a/Gr63a (V) 0.86 |
| Ethanol | weak everywhere (≤0.11) |

### `hallem2006_spikes`
- Source: Hallem & Carlson 2006 *Cell* 125:143. 24 ORs tested in the empty-neuron system at 10⁻² dilution.
- **Important:** DoOR's `Hallem.2006.EN` column stores **absolute** firing rates (net response plus spontaneous rate), even though its metadata says spontaneous is subtracted.
  - I checked this against the Hallem columns reproduced in Knaden et al. 2012 Table S1. The two match exactly after subtracting spontaneous.
  - The file therefore gives **net** spikes/s plus the spontaneous rates.
- Or42b and Or92a are not in this dataset.
- Responses below about 30–50 spikes/s are weak.

### `dose_response`
**From the literature:**
- **Population rule:** one Hill function with n ≈ 1.42 for calcium signals; EC50s follow a power-law distribution (Si et al. 2019 *Neuron*).
- **Or22a EC50s:** Pelz 2006, taken from DoOR columns `Pelz.2006.AntEC50/ALEC50`. I read the values as −log10 of the EC50 dilution (e.g. isoamyl acetate about 10⁻⁴, ethyl hexanoate about 10⁻⁶·⁷). That reading of the column is my inference.
- **ORN → PN transform:** Olsen, Bhandawat & Wilson 2010. PN = Rmax·ORN^1.5/(ORN^1.5 + σ^1.5 + (m·LFP)^1.5), with Rmax 165, σ 12, m 10.63; DM1 has σ 44.8.
- **ORN dynamics:** filter peak 56 ms, Weber-Fechner adaptation (Martelli 2013; Nagel & Wilson 2011; Gorur-Shandilya 2017).
- **CO2:** about half-maximal at 5% in the empty neuron (Kwon 2007). Behavioral avoidance starts at +0.1% (Suh 2004; Faucher 2006). Above about 5%, avoidance runs through Ir64a/DC4 (Ai 2010).

**Estimates (no fitted value was accessible):**
- EC50 for Or42b to ethyl acetate: about 10⁻⁵·⁵ dilution.
- EC50 for Or92a to diacetyl: about 10⁻⁵.
- EC50 for Ir64a/Ir75a to acetic acid: about 0.3–1%.
- EC50 for native ab1C to CO2: about 0.5–2%.

**Converting dilutions to ppm:** use mole fraction × Psat × carrier dilution. Treat the result as an **estimate** good to about ±1 log unit.

## Gustatory (`gustatory`)

**Sugar**
- Anchor points: 100 mM sucrose gives 58.9 spikes/s (Cameron 2010); 50 mM gives 43 spikes/s (Dweck 2022).
- I fitted a Hill curve through those two points: EC50 ≈ 59 mM, Rmax ≈ 94 spikes/s, n = 1. The two points come from different labs, so this fit is **derived**.
- Glucose and labellar fructose EC50s are **estimates**.
- Gr43a (internal fructose sensor) responds from about 5 mM, strongly at 25–100 mM. Hemolymph fructose is 1.9 mM and rises to 6–19 mM after a meal (Miyamoto 2012).

**PER and bitter**
- The sucrose concentration giving about 50–60% PER is 800 mM fed, 300 mM after 1 day starved and 200 mM after 2 days (Inagaki 2014).
- Bitter: standard test concentrations only. Caffeine EC50 of about 5 mM is an **estimate**.
- There are **no Drosophila taste data for tannins or polyphenols**. I model them as a weak Gr66a agonist with a threshold of 0.1–1 mM (**estimate**). Ripe pulp (about 0.03 mM) sits below that threshold; unripe pulp (about 1 mM) sits at it.

**Water (ppk28)**
- Water alone gives 12 spikes/s.
- Heterologous (HEK-cell) data show activation below about 216 mOsm.
- Half-inhibition at about 150 mOsm is an **estimate**. Banana juice (about 1000+ mOsm) should silence these neurons.

**Salt**
- Attraction at 1–100 mM, peaking at 50 mM; rejection at ≥200 mM (Zhang 2013).
- Cell classes from Jaeger 2018 and Dweck 2022 (Ir56b in sweet GRNs); Ir7c from McDowell 2022.
- Banana is essentially Na-free (0.6 mM).

**Acids**
- Bitter GRNs start firing at pH ≤5 (Charlu 2013).
- Two-choice acetic acid taste curve: ≤0.5% indifferent, 1% mildly preferred, ≥5% avoided via Ir7a (Rimal 2019).
- OtopLA carries HCl rejection (Ganguly 2021).
- Tarsal Ir25a/Ir76b neurons sense acid (Chen & Amrein 2017).
- Lactic acid acts on sweet GRNs via IR25a (Stanley 2021).

**Amino acids, fats, glycerol, carbonation, ethanol**
- Amino acids and yeast: Ganguly 2017, Aryal 2022, Steck 2018 (taste pegs: yeast response > 500 mM sucrose).
- Fatty acids: Ahn 2017.
- Glycerol: Gr64e (Wisotsky 2011). The EC50 is an **estimate**.
- Carbonation: taste pegs (Fischler 2007). No numbers found.
- Ethanol taste: aversive (no PER at 0.1–40%), yet flies drink 15% ethanol food (Devineni & Heberlein 2009).

## Internal state (`internal_state_modulation`)

**Fold changes as reported**

| Effect | Magnitude | Source |
|---|---|---|
| Sugar sensitivity, 2 days starved | ×4.6 | Inagaki 2014 |
| Bitter sensitivity, 2 days starved | ÷4.8 | Inagaki 2014 |
| Sugar + bitter mixture acceptance | ×10.2 | Inagaki 2014 |
| Dopamine on sugar GRN calcium signal | ×1.3–1.4 | Inagaki 2012 |
| TH-VUM dopamine neuron firing | 1 → 25 Hz | Marella 2012 |
| DM5 (Or85a) activation threshold, starved | ×4 | Ko 2015 |
| Acetic acid calcium signal, sugar GRNs (fed → starved) | ×1.42 | Devineni 2019 |
| Acetic acid calcium signal, bitter GRNs (fed → starved) | ×1.89 | Devineni 2019 |
| Yeast feeding, amino-acid deprived | ×1.6, first long yeast visit ~3× sooner | Corrales-Carvajal 2016 |

Also from Devineni 2019: PER to 10% acetic acid is at water level in fed flies and 86% in starved flies.

**Given only in figures:**
- DM1 facilitation under starvation (Root 2011). I used about ×1.75 as an **estimate**.
- Yeast-GRN gain after protein deprivation (Steck 2018). I used about ×1.5 as an **estimate**.
- Effect sizes in Ribeiro & Dickson 2010 and Walker 2015.

**Conflict:** Devineni 2019 saw no hunger effect on bitter GRN output, unlike Inagaki 2014 and LeDue 2016 (which found octopamine-mediated depotentiation).

`suggested_state_gains_for_model_ESTIMATE` turns these numbers into multiplicative gains for the simulation.

## GRN type → modality (`grn_type_modality`)

**Sources**
- Primary: Tastekin et al. 2025 bioRxiv 10.1101/2025.08.25.671814, the taste-feeding connectome. It defines the LB/tpGRN/PhG/LgAG/LgLG names by matching GAL4 morphology.
- Cross-checked against the FlyWire v783 `cell_sub_class` annotations (Schlegel 2024; github flyconnectome/flywire_annotations).
- Shiu et al. 2024 root IDs for sugar, water, bitter and Ir94e neurons were mapped onto FlyWire types (github philshiu/Drosophila_brain_model).
- Also Engert 2022.

**Confident assignments**

| Type | Modality |
|---|---|
| LB1a–d | bitter |
| LB1e | Ir94e (mildly aversive) |
| LB3 (FlyWire, unsplit) | sugar + water |
| LB3a | water |
| LB3b | sugar + low salt (Ir56b) |
| LB3c | sugar |
| LB3d | high salt / aversive |
| claw_tpGRN | carbonation (Ir56d) |
| dorsal_tpGRN | amino acids / yeast |
| PhG1a–c | glycerol (Gr64e) |
| LgAG1 | bitter |
| LgAG2 | sugar |
| LgLG types | ppk23/ppk25/Ir52 pheromone neurons; LgLG3/4 are sugar |

**Uncertain**
- **LB2a–c:** FlyWire's "low-salt" label rests on morphology alone. Tastekin places them in an unknown, aversive cluster, so do not treat them as attractive.
- **LB2d:** the sources conflict.
- **LB4a/b:** probably appetitive, modality unknown.
- **Most PhG and LgAG3–9 types:** assigned by receptor match or clustering only.

**Second-order neurons** (FlyWire names):
- Sugar: G2N-1 = CB0616, FMIn = CB0366, Usnea = CB0008, Phantom = CB0062.
- Water: Fudog = DNg67.
- Bitter: the main shared partner is the GABAergic CB0159. It is probably GNG016 (my inference).
- Aversive chain from Tastekin: GNG016 → GNG510 → Scapula (GNG087).
