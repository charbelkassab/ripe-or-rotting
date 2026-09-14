# chemistry.json: sources and assumptions

The parameter file is `chemistry.json`. Every number in it has a source and a confidence tag next to it:

- **high**: I checked the number in a fetched primary paper or database.
- **medium**: from a checked abstract or secondary summary, or a single study.
- **low**: from a search snippet or a different cultivar/fruit.
- **ESTIMATE** / **derived**: my own estimate or calculation. Always replace these when better data turn up.

## 1. Banana pulp time course (`banana_pulp_time_course`)

**Sources**

| Topic | Source | Notes |
|---|---|---|
| Base composition | USDA FoodData Central "Bananas, raw" (FDC 173944 / SR Legacy 09040) | Water 74.9 g, sucrose 2.39, glucose 4.98, fructose 4.85, starch 5.38, protein 1.09 g/100 g; K 358 mg, Na 1 mg. It is a composite with no ripeness stage, so treat it as generic "ripe". |
| Ripening trajectory | Phillips et al. 2021 *PLoS ONE* 16:e0253366 | Retail Cavendish, slightly ripe → ripe → overripe. Starch 4.5 → 2.5 → 0.5 g/100 g. Sucrose peaks at 4.6 and falls to 1.9. Hexoses reach about 13–14 g/100 g. |
| Tannins / phenolics | Bashmil et al. 2021 *Antioxidants* 10:1521 | Ripe pulp tannin 0.02 mg CE/g against 0.66 unripe, a roughly 30× drop. Peel keeps 0.19–0.50 mg/g. |
| Respiration CO2 | Faucher et al. 2006 *J Exp Biol* 209:2739 | Yellow banana about 200 µL CO2/min; black banana about 30 µL/min. |

**Conversions and estimates**

- mM values assume 0.749 L of pulp water per kg of fresh weight. For ripe pulp that gives glucose about 369 mM, fructose about 359 mM, sucrose about 93 mM, K⁺ about 122 mM and Na⁺ about 0.6 mM.
- The first-order rate constants for starch and sucrose were fitted by eye to the Phillips trajectory, so they are **estimates**.

**Weak spots**

- **pH and acids:** pH 4.5–5.2 in ripe pulp and malic acid as the main acid both come from secondary sources (Wyman & Palmer 1964). The g/100 g values for malic and citric acid are estimates.
- **Free amino acids:** there is no accessible total. Alsmairat, Engelgau & Beaudry 2018 (*J Am Soc Hortic Sci*) measured them, but only in figures. I entered about 50–150 mg/100 g as an **estimate**.
- **Volatiles (biggest gap):** I could not get absolute µg/kg values. Jordán et al. 2001 and Pino & Febles 2013 are paywalled, and their abstracts give no numbers.
  - The only anchor is total esters of about 20–24 mg/kg in Dwarf Cavendish.
  - The ranking of compounds does come from the literature: isoamyl acetate, 2-pentyl acetate, isoamyl butanoate and isobutyl acetate dominate. C6 aldehydes fall with ripening, while esters and ethanol rise.
  - The per-compound mg/kg defaults are therefore **order-of-magnitude estimates**, scaled so that ripe pulp adds up to about 10–25 mg/kg of esters.
  - The ethanol series (250–375 → >600 ppm) comes from a search summary whose primary paper I could not identify.

## 2. Fermentation (`fermentation`)

**Stoichiometry.** Glucose → 2 ethanol + 2 CO2 (theoretical yields 0.511 and 0.489 g/g). Ethanol → acetic acid (theoretical 1.304 g/g).

**Kinetic parameters**

- *S. cerevisiae* (anaerobic, 30 °C): μmax 0.31 h⁻¹, Ks 0.099 g/L, biomass yield 0.10 g/g (Verduyn et al. 1990).
- Aerobic, respiratory growth: biomass yield 0.50 g/g. Above the Crabtree threshold: 0.16 g/g (Postma et al. 1989).
- Values at 25 °C are scaled from 30 °C assuming Q10 ≈ 2 (**estimate**).
- There is no numeric μmax for *Hanseniaspora uvarum* or Acetobacter; I used about 0.4 h⁻¹ and about 0.1 h⁻¹ as **estimates**.

**Ethanol measured in natural fruit**

| Fruit / system | Ethanol | Source |
|---|---|---|
| Astrocaryum palm | 0.9% ripe, 4.5% overripe | Dudley 2004 |
| Spondias fruit dropped by spider monkeys | 1.6–1.9%, with a regression against Brix | Campbell et al. 2022 |
| Chimpanzee diet fruits | 0.31% | Maro et al. 2025 |
| Bertam palm nectar | 0.5–3.8% | Wiens et al. 2008 |
| Fruits that produced *D. melanogaster* | about 60% had little or none; the rest mostly 1–4% | McKechnie & Morgan 1982, via summary |

**Banana-specific ferments**

- Spontaneous banana juice ("Tonto"): 0 → 4.7% v/v ethanol and pH 4.9 → 3.6 within 60 h (Atwine et al. 2026).
- Banana wine made into vinegar: about 49 g/L acetic acid (Tanaka et al. 2016 preprint).

**Estimates**

- There is no direct measurement of ethanol or acetic acid in naturally rotting banana. The "rotting banana" day-by-day defaults for ethanol, acetic acid and pH are **estimates** bounded by the data above.
- Ethyl acetate, acetoin and 2-phenylethanol ranges come from wine and yeast studies (Hanseniaspora ethyl acetate 30–118 mg/L). They are liquid-culture values, and solid fruit probably reaches lower concentrations.

**Model sketch.** `suggested_kinetic_model_ESTIMATE` is a standard Monod-type model: starch → sucrose → hexose, yeast growth with ethanol inhibition, surface acetic acid bacteria, and volatile loss through an air-side mass-transfer term. Its defaults come from the rows above.

## 3. Physical constants (`volatile_physical_constants`)

**Henry's law.** From the Sander 2023 compilation (henrys-law.org), using review or measured values rather than QSPR estimates.
- Kaw = 1/(Hscp·R·T) at 298.15 K.
- Acetic acid: the Kaw applies to the undissociated acid. Multiply by 1/(1+10^(pH−4.76)).
- Acetoin is the least certain constant (about ×3).
- I replaced DoOR's HLC column with Sander values. For acetic acid the DoOR figure is about 2.5× lower.

**Diffusion coefficients.**
- Measured values come from Tang et al. 2015 (*ACP*), converted as D[cm²/s] = D[Torr cm²/s]/760. This covers ethanol, acetic acid, ethyl acetate, butyl and isobutyl acetate, and ethyl butyrate.
- Other compounds use the Fuller (1966) method, corrected with measured analogs. Isoamyl acetate uses n-pentyl acetate; isoamyl alcohol uses 1-pentanol. Hexanal, acetoin, 2-heptanone and 2-phenylethanol get Fuller −8 to −10%. These are **estimates** good to about ±10%.

**Vapor pressures.** From DoOR.data `odor.csv` (VP.25). They match standard values for ethanol (59.3), ethyl acetate (93.2), acetic acid (15.7) and isoamyl acetate (5.6 mmHg).

**Headspace examples.** The equilibrium headspace values are **derived upper bounds** from C_air = Kaw·C_aq. For example, 2% ethanol in pulp gives about 3000 ppmv, and ripe-level isoamyl acetate gives about 14 ppmv.

## 4. Behavioral anchors

This file keeps only a short list. The quantitative validation targets are in `behavior.json`, and the receptor and taste physiology is in `receptors.json`.

Two corrections to citations in the original brief:
- Joseph et al. 2009 is in *PNAS*, not *Science*.
- Zhu, Park & Baker 2003 studied overripe **mango**.

## Using the file in a simulation

- Treat pulp as a well-mixed aqueous phase at pH 4.9 that falls to about 3.6–3.8 under fermentation. Volatiles partition by Kaw at the surface, and only the undissociated fraction of acetic acid does.
- Tissue diffusion of volatiles is not given. It is roughly water-like, D ≈ 1×10⁻⁵ cm²/s (**estimate**), and this slows emission from intact pulp a great deal.
- Surface area and fruit mass are rough **estimates** (120–150 g fruit, about 200–250 cm² peel).
