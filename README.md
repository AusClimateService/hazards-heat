# ACS Hazard Team on Heatwaves and Extreme Heat

## Description
GitHub repository for the ACS Heat Team - heatwaves and extreme heat. 

## Indices considered by the hazard team:
- Excess Heat Factor (EHF)
- Heatwave frequency (HWF)
- Heatwave duration (HWD)
- Heatwave number (HWN)
- Heatwave peak temperature (HWAtx)
- Annual mean daily maximum temperature (TXm)
- Annual maximum daily maximum temperature (TXx)
- Days above 35C (TXge35)
- Days above 40C (TXge40)
- Days above 45C (TXge45)
- Days below 2C (TNle02)

## Products:
Status of the NCRA deliverables. 

The two dots (in order from first/top/left to last/bottom/right) represent the datasets used to compute indices:
- Dot 1: Pre-processed BARPA/CCAM – downscaled but NOT bias-corrected, 5 km (deliverable for 30 June)
- Dot 2: Bias-corrected BARPA/CCAM – downscaled AND bias-corrected, 5 km (deliverable for 31 July)
Where only one dot is in the cell the format type does not apply to the metric.
 
In terms of the colors:
- :green_circle: The data is available in its final official form
- :yellow_circle: The data creation is currently in progress and available soon
- :red_circle: The data processing has not yet started
- :white_circle: Not intended for delivery/not applicable

| Index/metric | variable_id | temporal resolution | reference | time series (ts) | GWLs ts | GWLs 2D | MME 2D | MME 2D change | Notes | Data<br>location | Last update
| -----        | -----       | -----               | -----     | :-:              |:-:      |:-:      |:-:     |:-:            |-----  |-----             |-----
|Heatwave frequency|HWF|annual|[Nairn and Fawcett 2015](https://www.mdpi.com/1660-4601/12/1/227)|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:||`/g/data/ia39/ncra/heat/<variable_id>`|30/06/2024
|Heatwave duration|HWD|annual|[Nairn and Fawcett 2015](https://www.mdpi.com/1660-4601/12/1/227)|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:||`/g/data/ia39/ncra/heat/<variable_id>`|30/06/2024
|Heatwave number|HWN|annual|[Nairn and Fawcett 2015](https://www.mdpi.com/1660-4601/12/1/227)|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:||`/g/data/ia39/ncra/heat/<variable_id>`|30/06/2024
|Heatwave peak temperature|HWAtx|annual|[Nairn and Fawcett 2015](https://www.mdpi.com/1660-4601/12/1/227)|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:||`/g/data/ia39/ncra/heat/<variable_id>`|30/06/2024
|Annual mean daily maximum temperature|TXm|annual|[ETCCDI Climate Change Indices](http://etccdi.pacificclimate.org/list_27_indices.shtml)|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:||`/g/data/ia39/ncra/heat/<variable_id>`|28/06/2024
|Annual maximum daily maximum temperature|TXx|annual|[ETCCDI Climate Change Indices](http://etccdi.pacificclimate.org/list_27_indices.shtml)|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:||`/g/data/ia39/ncra/heat/<variable_id>`|28/06/2024|
|Days above 35C|TXge35|annual|[ETCCDI Climate Change Indices](http://etccdi.pacificclimate.org/list_27_indices.shtml)|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:||`/g/data/ia39/ncra/heat/<variable_id>`|28/06/2024
|Days above 40C|TXge40|annual|[ETCCDI Climate Change Indices](http://etccdi.pacificclimate.org/list_27_indices.shtml)|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:||`/g/data/ia39/ncra/heat/<variable_id>`|28/06/2024
|Days above 45C|TXge45|annual|[ETCCDI Climate Change Indices](http://etccdi.pacificclimate.org/list_27_indices.shtml)|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:||`/g/data/ia39/ncra/heat/<variable_id>`|28/06/2024
|Days below 2C|TNle02|annual|[ETCCDI Climate Change Indices](http://etccdi.pacificclimate.org/list_27_indices.shtml)|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:|:green_circle:<br>:green_circle:||`/g/data/ia39/ncra/heat/<variable_id>`|28/06/2024

## Sample Results

Hazard metrics have been computed and results are available on NCI (see directory paths in above table).

In addition to the data files (in netcdf format), results have also been presented in the form of maps (also available on NCI) and, for some hazard metrics, tables summarising regional changes in the aggregated hazard indice(s).

![Encemble central estimate for the change in TXx between GWL 3.0 and GWL1.2](figures/TXx_AGCD-05i_MME50_ssp370_v1-r1-ACS-QME-AGCD-1960-2022_GWL30-GWL12-change.png)
*Figure 1. Sample map showing the change in **TXx** (hottest day of the year) between GWL 3.0 and GWL 1.2. Here, results are presented as the ensemble median (i.e., the ensemble central estimate) for the 13 ACS regional climate model simulations used in the analysis.*

![Regional changes in TXx between GWL 3.0 and GWL1.2](figures/TXx_regional_summary_GWL30-GWL12.png)
*Figure 2. Sample heatmap showing the change in **TXx** (hottest day of the year) between GWL 3.0 and GWL 1.2, averaged accross the different NCRA study regions. Here, results are presented for each of the 13 ACS regional climate model simulations, as well as the ensemble median.*

![Encemble central estimate for the change in HWF between GWL 3.0 and GWL1.2](figures/HWF_AGCD-05i_MME50_ssp370_v1-r1-ACS-QME-AGCD-1960-2022_GWL30-GWL12-change.png)
*Figure 3. Sample map showing the change in **HWF** (days per year experiencing heatwave conditions) between GWL 3.0 and GWL 1.2. Here, results are presented as the ensemble median (i.e., the ensemble central estimate) for the 13 ACS regional climate model simulations used in the analysis.*

![Regional changes in HWF between GWL 3.0 and GWL1.2](figures/HWF_regional_summary_GWL30-GWL12.png)
*Figure 4. Sample heatmap showing the change in **HWF** (days per year experiencing heatwave conditions) between GWL 3.0 and GWL 1.2, averaged accross the different NCRA study regions. Here, results are presented for each of the 13 ACS regional climate model simulations, as well as the ensemble median.*

## Authors and acknowledgment
Hazard team:
- [ ] Mitchell Black (BOM, lead)
- [ ] Matt Beaty (ABS, alternate lead)
- [ ] Cass Rogers (BOM, contributor)
- [ ] David Hoffmann (BOM, contributor)
- [ ] Jessica Bhardwaj (BOM, contributor)
- [ ] Steph Jacobs (BOM, contributor)
- [ ] Aurel Griesser (BOM, contributor)
- [ ] Naomi Benger (BOM, contributor)

