# W TagAndProbe rhalphalib implementation

### Setup combine and rhalphalib
Follow official combine instruction or setup with `conda` using [this branch](https://github.com/andrzejnovak/HiggsAnalysis-CombinedLimit/tree/root6.22-compat).

Setup [rhalphalib](https://github.com/nsmith-/rhalphalib)


### Generate variations 
Assumes 1 root file per category containing matched, unmatched and data templates.
For each root file generate variations (only matched - catp2)
```
python scalesmear.py -i templates/wfit_nskim17_n2/wtag_pass.root  --plot
python scalesmear.py -i templates/wfit_nskim17_n2/wtag_fail.root  --plot
```

### Generate combine/rhalphalib workspace and fit
```
# python sf.py
python sf.py --fit single --tp templates/wfit_nskim17_n2/wtag_var_pass.root --tf templates/wfit_nskim17_n2/wtag_var_fail.root
cd fitdir
bash build.sh
bash t2w.sh # Set SF as POI
```

and run the fit
```
combine -M FitDiagnostics --expectSignal 1 -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 --justFit
```

