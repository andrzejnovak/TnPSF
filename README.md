# W TagAndProbe rhalphalib implementation

### Setup combine and rhalphalib
Follow official combine instruction or setup with `conda` using [this branch](https://github.com/andrzejnovak/HiggsAnalysis-CombinedLimit/tree/root6.22-compat).

Setup [rhalphalib](https://github.com/nsmith-/rhalphalib) and then

```bash
git clone https://github.com/andrzejnovak/TnPSF.git
cd TnPSF
```

### Generate variations 
For each root file generate variations (only matched - catp2)
```bash
python scalesmear.py -i templates/wfit_nskim17_n2/wtag_pass.root  --plot
python scalesmear.py -i templates/wfit_nskim17_n2/wtag_fail.root  --plot
```

### Generate combine/rhalphalib workspace and fit

```bash
python sf.py --fit single -t templates/ref18/wtemplates_n2cvb_var.root -o FitSingle
cd FitSingle
```

or for two-cut setup:
```
python sf.py --fit double -t templates/ref18/wtemplates_n2cvb_var.root --t2 templates/ref18/wtemplates_cvl_var.root -o FitDouble
cd FitDouble
```
and run the fit
```
combine -M FitDiagnostics --expectSignal 1 -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 --saveShapes --saveWithUncertainties --rMin 0.5 --rMax 1.5
```

To make plots from `FitDiagnostics` output, run within the fit folder:
```
python ../results.py --year 2018
```

