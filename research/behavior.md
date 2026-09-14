# behavior.json: sources and assumptions

This file lists behavioral targets for checking the smell and taste model against real flies.

Confidence tags used below and in the JSON:
- **VERIFIED**: number read in the paper's full text, abstract or supplement.
- **SUMMARIZER**: number taken from a page summary.
- **FIGURE-ONLY**: the paper only shows the value in a figure. The value is not given in the file.
- **ESTIMATE**: my own number.

## 1. Knaden et al. 2012, *Cell Reports* 1:392 (`knaden2012_valence`)

**Where the data came from.** Tables S1 and S2 were downloaded from MPG PuRe, which is open access:
- Table S1: https://pure.mpg.de/rest/items/item_1449055_9/component/file_1525217/content
- Table S2: https://pure.mpg.de/rest/items/item_1449055_9/component/file_1525216/content

`all_110_odorants` contains the median attraction index (AI), the Wilcoxon p value and the valence class for every odorant.

**Assay (it is not a T-maze).**
- Two traps, one odor and one solvent, left for 24 h.
- 50 flies per run, starved for 24 h, 10 replicates.
- AI = (O − C)/50, where O and C are the flies caught in the odor and control traps.
- Dilution was about 1:100 in water. The µ symbol is garbled in the PDF, so this ratio is my interpretation.

**Results for the banana and fermentation compounds.**

| Result | Odorants (AI) |
|---|---|
| Attractive | acetic acid +0.24, 2-phenylethanol +0.44, phenylethyl acetate +0.54, isoamyl alcohol +0.31, hexanal +0.26, diacetyl +0.62, ethyl propionate +0.45, methyl/propyl acetate +0.24, eugenol +0.33 |
| Neutral | ethyl acetate +0.18 (p = 0.06), isoamyl acetate −0.11, butyl acetate −0.10, ethanol +0.05, 2-heptanone +0.07, ethyl butyrate +0.08, ethyl hexanoate +0.02 |
| Aversive | benzaldehyde −0.53, linalool −0.42, 1-octen-3-ol −0.28 |
| Not tested | acetoin, geosmin, CO2 |

**Glomerular valence (projection neurons).**
- DM4, DM5 and DM2 respond more to attractive odors.
- D, DA4, DC3, DL1, DL4 and DL5 respond almost only to aversive odors.
- DA4 and DC3 are aversive-specific at both the sensory-neuron and projection-neuron levels.

**Caveats.**
- Every odorant was tested at one concentration, and the trap assay integrates over 24 h.
- DM5 comes out "attractive" here but "aversive" in Semmelhack & Wang 2009, which used an immediate olfactometer assay. Weight the Knaden AIs as trap-entry tendencies, not as fixed valences.

## 2. Concentration dependence (`concentration_dependence`)

**Apple cider vinegar**
- Semmelhack & Wang 2009 (VERIFIED):
  - Preference index 75% at 3 ppm, slightly higher at 12 ppm, 9% at 32 ppm.
  - At 32 ppm with DM5 silenced the index is 87%.
  - Silencing DM1 drops it to −4%; silencing VA2 drops it to 50%.
- Ko 2015: starved flies' appetitive index rises from 23 to 60 over 0.5–25% vinegar, and starved flies score above fed flies at every concentration.
- Álvarez-Salvado 2018 (walking flies): upwind velocity follows a Hill curve with Kd 0.072% vinegar and n 1.03, plus a fitted navigation model.

**Acetic acid**
- Joseph 2009 (a *PNAS* paper, not *Science*):
  - At 5%: oviposition index +0.82 while position index is −0.33.
  - At 0.25%: oviposition +0.34, position −0.03.
- Ai 2010: flies avoid acid through DC4. Vinegar at pH 2.5 is avoided, but the same vinegar neutralized is attractive.
- Devineni 2019 (hunger switch): fed flies show water-level PER to acetic acid; starved flies show 86% PER to 10%.
- Rimal 2019 (taste two-choice): up to 0.5% indifferent, 1% mildly preferred, 5% and above avoided.

**Ethanol**
- Ogueta 2010: flies prefer food with up to 5% ethanol and avoid 23%.
- Azanchi 2013: egg-laying preference is significant from 3% and peaks at 5%.
- Schneider 2012: in a trap assay, preference index 0.87 for food plus ethanol against ethanol alone.
- van Breugel 2018: ethanol is attractive in every activity state.
- The relative-attraction curve in `model_curve_ESTIMATE` combines these sources and is my ESTIMATE.

**CO2**
- Avoidance starts around +0.1% (Suh 2004; Faucher 2006). Bräcker 2013 reports avoidance above +0.02%.
- Hunger combined with food odor suppresses avoidance (Bräcker 2013; Lewis 2015).
- Active or flying flies are attracted instead. The attraction depends on Ir25a, and 5% works best (van Breugel 2018; Wasserman 2013 in tethered flight).
- Most avoidance and attraction indices are FIGURE-ONLY.

## 3. Contact feeding (`contact_feeding_preferences`)

**Ethanol food**
- Devineni & Heberlein 2009: in the CAFE assay, preference for 15% ethanol food builds over 4–5 days. Preference indices are FIGURE-ONLY.
- Pohl 2012: no preference at 1%, preference at 5–15%, and none against a calorie-matched alternative.

**Yeast versus sucrose**
- Corrales-Carvajal 2016 (VERIFIED): amino-acid deprivation raises yeast feeding 1.6-fold. The first long yeast visit comes at 4.4 min instead of 12.4 min. Mated females shift feeding toward yeast.
- Lee 2008: flies self-select a protein:carbohydrate ratio of about 1:4.
- Vargas 2010, Ribeiro & Dickson 2010 and Walker 2015 report the direction of the effect in text; the size is FIGURE-ONLY.

**Sugar**
- Dus 2011: sugar-blind mutants still come to prefer nutritive sugar after 22 h of starvation.
- PER dose-response data are in `receptors.json`.

## 4. Navigation anchors

| Anchor | Value | Source |
|---|---|---|
| Upwind surge after plume contact | 190 ms | van Breugel & Dickinson 2014 |
| Start of casting after plume loss | 450 ms | van Breugel & Dickinson 2014 |
| Walking in intermittent plumes | encounter-rate-dependent upwind bias of −21.6°/Hz | Demir 2020 |
| Dispersal | up to about 12 km in one flight, ~1 m/s groundspeed | Leitch 2021 |

No study measures directly how far away a fly can detect a fruit. Derive detection distance from plume physics, the source strengths in `chemistry.json`, and the receptor thresholds in `receptors.json`.

## 5. Suggested validation tests

`suggested_validation_tests_ESTIMATE` lists 8 qualitative predictions the simulation should reproduce, each with its source:
- Ripe versus rotting banana.
- The vinegar dose sweep with DM5 recruitment.
- Loss of attraction when DM1 is silenced.
- Geosmin overriding attraction.
- Ethanol preference peaking near 5%.
- CO2 avoidance in resting flies.
- The shift to yeast in mated, protein-deprived females.
- The hunger switch for acetic acid taste.
