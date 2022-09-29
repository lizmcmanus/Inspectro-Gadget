# Inspectro-Gadget

A tool to obtain gene expression values from a defined brain region for a set of neurotransmitter receptor subunits.

Gene expression data is obtained from the Allen Human Brain Atlas (https://human.brain-map.org/).

## Installation

The package can be installed via pip as follows:

```
python -m pip install https://github.com/lizmcmanus/Inspectro-Gadget/archive/main.tar.gz
```



## Usage
### Single region

```
gadget_out =  gadget(['/home/data/region1_mask.nii.gz'],
              mask_labels=['Region 1'],
              out_root='/home/data/results')
```

### Two region comparison

```
gadget_out = gadget([['/home/data/region1_mask.nii.gz'],
                    ['/home/data/region2_mask.nii.gz']],
                    mask_lables=['Region 1', 'Region 2'],
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

_Whole-brain gene expression images_
