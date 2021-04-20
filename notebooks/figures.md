---
title: "CESM2-MARBL"
tags: presentation
slideOptions:
  theme: white
  transition: fade #none/fade/slide/convex/concave/zoom
  progress: true
  center: true  
  minScale: 0.2
  maxScale: 1
---

<style type="text/css">
  .reveal h1 {
    text-align: left;    
    font-size: 150%;
    text-shadow: none;
    font-weight: 500;
  }
  .reveal h2 {
    text-align: left;
    font-size: 120%;
    font-weight: 500;
  }
  .reveal h3 {
    text-align: left;
    font-size: 100%;
    font-weight: 500;
  }  
  .reveal ul {
    display: block;
    text-align: left;
  }
  .reveal ol {
    display: block;
  }
  td, th {
     border: none
  }  
</style>

# The Marine Biogeochemistry Library (MARBL)

Matthew C. Long, Keith Lindsay, J. Keith Moore, 
Michael Levy, Kristen Krumhardt, 
Jessica Luo, Scott C. Doney, [Your name here?]

---

## Paper outline
- Review history of BEC-->MARBL development
- Description of model formulations
    - Phytoplankton growth, zooplankton grazing
    - Detrital organic pools
    - Nitrogen cycle
    - Iron cycle
    - Riverine forcing and sediments
    - Dissolved oxygen (the kludge)
    - MARBL features not enabled in CMIP6
- **Results: Analyze results from CESM2-MARBL CMIP6 coupled integrations**

---

## Results
- Nutrient simulation
- Surface chlorophyll
- Limiting nutrient
- Nitrogen cycle balance
- Ocean circulation biases (<sup>14</sup>C, CFCs, C<sub>ant</sub>)
- Surface CO<sub>2</sub> flux 
- Transient response: NPP, export production, air-sea CO<sub>2</sub> flux

---

[//]: # (AUTO-GENERATED CONTENT BELOW)

### Ventilation biases: C<sub>ant</sub>

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/za-obs-comparison-Cant.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### Ventilation biases: <sup>14</sup>C

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/za-obs-comparison-Del14C.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### Ventilation biases: pCFC-12

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/za-obs-comparison-pCFC12.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### Ventilation biases: pCFC-11

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/za-obs-comparison-pCFC11.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### MOC

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/moc.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### Mixed layer depth

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/mld-maps.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### Surface chlorophyll

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/chl-surface-maps.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### Nutrient limitation

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/nutrient-limitation-maps.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### Biological pump: NPP & Export

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/bio-pump-NPP-Export-e_ratio-N2P.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### Biological pump: Mineral fluxes

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/bio-pump-CaCO2-SiO2.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### DOM

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/DOM-ocean-obs-PDF.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### DOM maps

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/DOM-concentration-maps.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### Surface nutrients

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/nutrients-surface-maps.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### Nutrient sections: NO3

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/nutrients-sections-NO3.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### Nutrient sections: SiO3

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/nutrients-sections-SiO3.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### Nutrient sections: PO4

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/nutrients-sections-PO4.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### Nutrient profiles

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/nutrients-global-profiles.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### N Cycle

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/Ncycle.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### Fe Cycle

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/iron-budget-maps.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### Fe obs comparison: PDF

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/iron-global-ocean-obs-PDF.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### Fe obs comparison: maps

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/iron-concentration-maps.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### Transient response of the biological pump

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/global-timeseries-bio-pump.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---

### Transient CO<sub>2</sub> uptake

<img src="https://raw.githubusercontent.com/marbl-ecosys/cesm2-marbl/master/notebooks/figures/global-timeseries-FG_CO2.png" 
     style="background:none; border:none; box-shadow:none;" 
     width=1000px>

---
