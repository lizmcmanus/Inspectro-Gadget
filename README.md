# Inspectro-Gadget

A tool to obtain gene expression values from a defined brain region for a set of neurotransmitter receptor subunits.

Gene expression data is obtained from the Allen Human Brain Atlas (https://human.brain-map.org/).

## Installation

The package can be installed via pip as follows:

```
pip install https://github.com/lizmcmanus/Inspectro-Gadget/archive/refs/heads/main.zip
```


### Important limitations


## Usage
### Single region

```
from inspectro_gadget.gadget import gadget

gadget_out =  gadget(['/home/data/region1_mask.nii.gz'],
              mask_labels=['Region 1'],
              out_root='/home/data/results')
```

### Two region comparison

```
gadget_out = gadget([['/home/data/region1_mask.nii.gz'],
                    ['/home/data/region2_mask.nii.gz']],
                    mask_labels=['Region 1', 'Region 2'],
                    out_root='/home/data/results')
```

### Multiple subject comparison

```
gadget_out = gadget(['/home/data/sub-01_region1_mask.nii.gz',
                     '/home/data/sub-02_region1_mask.nii.gz',
                     '/home/data/sub-03_region1_mask.nii.gz',
                     '/home/data/sub-04_region1_mask.nii.gz'],
                    mask_labels=['sub-01', 'sub-02', 'sub-03', 'sub-04'],
                    out_root='/home/data/results')
```

## Citing
When reporting results obtained with the tool please cite the following work:

_Inspectro-Gadget_

_Allen Human Brain Atlas_

Hawrylycz, M.J., Lein, E.S., Guillozet-Bongaarts, A.L., Shen, E.H., Ng, L., Miller, J.A., van de Lagemaat, L.N., Smith, K.A., Ebbert, A., Riley, Z.L., Abajian, C., Beckmann, C.F., Bernard, A., Bertagnolli, D., Boe, A.F., Cartagena, P.M., Chakravarty, M.M., Chapin, M., Chong, J., Dalley, R.A. et al. An anatomically comprehensive atlas of the adult human brain transcriptome. Nature, 489(7416):391–399, 2012. doi: 10.1038/nature11405

_Whole-brain gene expression images_

Gryglewski, G., Seiger, R., James, G.M., Godbersen, G.M., Komorowski, A., Unterholzner, J., Michenthaler, P., Hahn, A., Wadsak, W., Mitterhauser, M., Kasper, S. and Lanzenberger, R. Spatial analysis and high resolution mapping of the human whole-brain transcriptome for integrative analysis in neuroimaging. NeuroImage, 176:259–267, 2018. doi: 10.1016/j.neuroimage.2018.04.068.