# InSpectro-Gadget

InSpectro-Gadget is a tool to obtain reference gene expression values from brain regions of interest for a set of neurotransmitter receptors and receptor subunits. Values are estimates for each voxel within the specified region and then presented graphically in a PDF report. Values are also returned as a Python dictionary for additional analysis.

The receptors and receptor subunits included are:

| Receptor           | Genes                                  |
|--------------------|----------------------------------------|
| GABA<sub>A</sub>   | GABRA1, GABRA2, GABRA3, GABRA4, GABRA5 |
| GABA<sub>B</sub>   | GABRB1, GABRB2, GABRB3                 |
|                    | GABRG1, GABRG2, GABRG3                 |
|                    | GABBR1, GABBR2                         |
| NMDA               | GRIN1                                  |
|                    | GRIN2A, GRIN2B, GRIN2C                 |
| AMPA               | GRIA1, GRIA2, GRIA3, GRIA4             |
| Kainate            | GRIK1, GRIK2, GRIK3, GRIK4, GRIK5      |
| mGlu               | GRM1, GRM5                             |
|                    | GRM2, GRM3, GRM4                       |
|                    | GRM5, GRM6, GRM7, GRM8                 |
| mACh               | CHRM1, CHRM2, CHRM3, CHRM4, CHRM5      |
| nACh               | CHRNA2, CHRNA3, CHRNA4, CHRNA5, CHRNA6, CHRNA7, CHRNA9, CHRNA10 |
|                    | CNRNB2, CHRNB3 |
| &alpha;-adrenergic | ADRA1D, ADRA1B, ADRA1A |
|                    | ADRA2A, ADRA2B, ADRA2C |
| &beta;-adrenergic  | ADRB1, ADRB2, ADRB3 |
| Dopamine | D1, D5 |
| | D2, D3, D4 |
| 5-HT | 5HT1A, 5HT1B, 5HT1D, 5HT1E, 5HT1F |
| | 5HT2A, 5HT2B, 5HT2C |
| | 5HT3A, 5HT4, 5HT5A, 5HT6, 5HT7 |

Gene expression data is obtained from the Allen Human Brain Atlas (https://human.brain-map.org/).

## Installation

The package can be installed via pip as follows:

```
pip install https://github.com/lizmcmanus/Inspectro-Gadget/archive/refs/heads/main.zip
```

## Usage
For each analysis, enter the path to your region or regions of interest. Label(s) for the regions and a specific output directory can also be specified but are not required.

Regions of interest should be binary images in 2 mm isotropic MNI152NLIN6Asym space. 

### Single region
To obtain receptor estimates for a single region enter that region's mask. A PDF report will be produced in the current working directory or in the output directory specified.

```
from inspectro_gadget.gadget import gadget

gadget_out =  gadget(['/home/data/region1_mask.nii.gz'],
              mask_labels=['Region 1'],
              out_root='/home/data/results')
```


![alt text](images/image.jpg?raw=true "Title")

### Two region comparison
To compare receptor estimates between two regions, enter masks for each. A PDF report will be produced in the current working directory or in the output directory specified.

```
gadget_out = gadget([['/home/data/region1_mask.nii.gz'],
                    ['/home/data/region2_mask.nii.gz']],
                    mask_labels=['Region 1', 'Region 2'],
                    out_root='/home/data/results')
```

### Multiple participant comparison
To analyse the similarity of receptor estimates across a group of participants, enter masks for each one. A participant overlap nifti image will be created in the current working directory or in the output directory specified, along with a PDF report. 

```
gadget_out = gadget(['/home/data/sub-01_region1_mask.nii.gz',
                     '/home/data/sub-02_region1_mask.nii.gz',
                     '/home/data/sub-03_region1_mask.nii.gz',
                     '/home/data/sub-04_region1_mask.nii.gz'],
                    mask_labels=['sub-01', 'sub-02', 'sub-03', 'sub-04'],
                    out_root='/home/data/results')
```

The report will include a variety of information:

- An image showing the overlap between participant masks

![Brain slices showing overlap between the participant masks.](images/multi-overlap.png?raw=true "Multiple participant mask overlap image")

- Estimates excitation/inhibition ratios for each participant's mask (and a group average value)

- ![Excitation/inhibition ratio estimates for each participant.](images/multi-eiratio.png?raw=true "Multiple participant excitation/inhibition ratios")

- Radarplots showing the receptor fingerprints for each participant's mask. Separate ones are plotted for GABA+glutanate genes and for neuromodulators

![Radarplots showing receptor fingerprints for participants.](images/multi-radars.png?raw=true "Multiple participant receptor fingerprints")

- Violinplots showing the range of gene expression within each participant's mask for each of the included genes. How far each participant mask's distribution median is from the group average is also provided

![Example violinplots showing expression levels for each of the included genes. A set of AMPA receptor genes are shown for illustration.](images/multi-violins.png?raw=true "Individual expression estimates for a set of AMPA receptor genes")

### Important limitations
The Allen Human Brain Atlas data comes from six people, five of whom were male. Their ages ranged from 24 to 57 years. As brain morphology and gene expression can differ between the sexes and across the lifespan, the dataset may not be fully representative of the human population. As such, the results produced by this tool should be interpreted with these limitations in mind. 

## Citing
When reporting results obtained with the tool please cite the following work:

_Inspectro-Gadget_

McManus, E., Muhlert, N., Duncan, N.W. InSpectro-Gadget: A tool for estimating neurotransmitter and neuromodulator receptor distributions for MRS voxels. _bioRxiv_ 2023. doi: 10.1101/2023.11.02.565296 

_Allen Human Brain Atlas_

Hawrylycz, M.J., Lein, E.S., Guillozet-Bongaarts, A.L., Shen, E.H., Ng, L., Miller, J.A., van de Lagemaat, L.N., Smith, K.A., Ebbert, A., Riley, Z.L., Abajian, C., Beckmann, C.F., Bernard, A., Bertagnolli, D., Boe, A.F., Cartagena, P.M., Chakravarty, M.M., Chapin, M., Chong, J., Dalley, R.A. et al. An anatomically comprehensive atlas of the adult human brain transcriptome. _Nature_, 489(7416):391–399, 2012. doi: 10.1038/nature11405

_Whole-brain gene expression images_

Gryglewski, G., Seiger, R., James, G.M., Godbersen, G.M., Komorowski, A., Unterholzner, J., Michenthaler, P., Hahn, A., Wadsak, W., Mitterhauser, M., Kasper, S. and Lanzenberger, R. Spatial analysis and high resolution mapping of the human whole-brain transcriptome for integrative analysis in neuroimaging. _NeuroImage_, 176:259–267, 2018. doi: 10.1016/j.neuroimage.2018.04.068.