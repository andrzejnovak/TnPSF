import ROOT as r
import uproot
import json
import os
import pprint
import numpy as np
import matplotlib
import warnings
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mplhep as hep

from matplotlib.offsetbox import AnchoredText
import argparse

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def convTH1(TH1):
    vals = TH1.values()
    # edges = TH1.edges
    edges = TH1.axes[0].edges()
    variances = TH1.variances()
    vals = vals * np.diff(edges)
    variances = variances * np.diff(edges)
    return vals, edges, variances


def covnTGA(tgasym):
        # https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/wiki/nonstandard
        # Rescale density by binwidth for actual value
        # _binwidth = tgasym._fEXlow + tgasym._fEXhigh
        # _x = tgasym._fX
        # _y = tgasym._fY * _binwidth
        # _xerrlo, _xerrhi = tgasym._fEXlow, tgasym._fEXhigh
        # _yerrlo, _yerrhi = tgasym._fEYlow * _binwidth, tgasym._fEYhigh * _binwidth
        # return _x, _y, [_yerrlo, _yerrhi], [_xerrlo, _xerrhi]
        _x, _y = tgasym.values(axis='both')
        _xerrlo, _xerrhi = tgasym.errors("low", axis='x'), tgasym.errors("high", axis='x')
        _yerrlo, _yerrhi = tgasym.errors("low", axis='y'), tgasym.errors("high", axis='y')
        _binwidth = _xerrlo + _xerrhi
        _y = _y * _binwidth
        _yerrlo = _yerrlo * _binwidth
        _yerrhi = _yerrhi * _binwidth
        return _x, _y, [_yerrlo, _yerrhi], [_xerrlo, _xerrhi]
        

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", default='', help="Model/Fit dir")
parser.add_argument("-i",
                    "--input",
                    default='fitDiagnosticsTest.root',
                    help="Input shapes file")
parser.add_argument("--fit",
                    default='fit_s',
                    choices={"prefit", "fit_s"},
                    dest='fit',
                    help="Shapes to plot")
parser.add_argument("-o",
                    "--output-folder",
                    default='plots',
                    dest='output_folder',
                    help="Folder to store plots - will be created if it doesn't exist.")
parser.add_argument("--year",
                    default=None,
                    choices={"2016", "2017", "2018"},
                    type=str,
                    help="year label")
parser.add_argument("--sflabels",
                    default="N2+CvB,CvL",
                    type=str,
                    help="labels")
parser.add_argument('-f',
                    "--format",
                    type=str,
                    default='png',
                    choices={'png', 'pdf'},
                    help="Plot format")
parser.add_argument("--scale",
                    default=1,
                    type=float,
                    help="Scale value, as used in template generation. (datacard scaling is parsed automatically)")
parser.add_argument("--smear",
                    default=0.5,
                    type=float,
                    help="Smear value, as used in template generation. (datacard scaling is parsed automatically)")

args = parser.parse_args()
if args.output_folder.split("/")[0] != args.dir:
    args.output_folder = os.path.join(args.dir, args.output_folder)

rd = r.TFile.Open(os.path.join(args.dir, args.input))
fd = uproot.open(os.path.join(args.dir, args.input))
with open(os.path.join(args.dir,'config.json')) as cfg_file:
    cfg = json.load(cfg_file)

if not os.path.exists(os.path.join(args.dir, 'plots')):
    os.mkdir(os.path.join(args.dir, 'plots'))

par_names = rd.Get('fit_s').floatParsFinal().contentsString().split(',')
par_names = [p for p in par_names if 'smear' in p or 'scale' in p  or 'eff' in p]

out = {}
for pn in par_names:
    out[pn] = {}
    out[pn]['val'] = round(rd.Get('fit_s').floatParsFinal().find(pn).getVal(), 3)
    out[pn]['unc'] = round(rd.Get('fit_s').floatParsFinal().find(pn).getError(), 3)

if np.isclose(abs(out['CMS_scale']['val']), 1):
    warnings.warn("CMS_scale is 1 (combine hist interpolation boundary)."
                  "Try inject additional scaling into datacards with --scale `>1`"
                  "when running sf.py")
if np.isclose(abs(out['CMS_smear']['val']), 1):
    warnings.warn("CMS_smear is 1 (combine hist interpolation boundary)."
                  "Try inject additional scaling into datacards with --smear `>1`"
                  "when running sf.py")


out['shift_SF'] = cfg['scale'] * out['CMS_scale']['val'] * args.scale  # (template shape)
out['shift_SF_ERR'] = cfg['scale'] * out['CMS_scale']['unc'] * args.scale  # (template shape)
out['smear_SF'] = 1 + cfg['smear'] * out['CMS_smear']['val'] * args.smear  # (template shape)
out['smear_SF_ERR'] = cfg['smear'] * out['CMS_smear']['unc'] * args.smear # (template shape)

if 'effSF' in out.keys():
    out['V_SF'] = out['effSF']['val']
    out['V_SF_ERR'] = out['effSF']['unc']
if 'effwSF' in out.keys():
    out['W_SF'] = out['effwSF']['val']
    out['W_SF_ERR'] = out['effwSF']['unc']

pprint.pprint(out)

plt.style.use([hep.style.ROOT])


shapetype = 'shapes_{}'.format(args.fit)
regions = [r.decode('utf').replace(";1", '') for r in fd[shapetype].keys() if "/" not in r]

lumi = {
    2016 : 36.33,
    2017 : 41.53,
    2018 : 59.74,
}

for i, reg in enumerate(regions):
    fig, (ax, rax) = plt.subplots(2,1, figsize=(10, 10), gridspec_kw = {'height_ratios':[3, 1]}, sharex=True)
    fig.subplots_adjust(hspace=0)
    ho1 = convTH1(fd[shapetype+'/'+reg+'/qcd'])
    ho2 = convTH1(fd[shapetype+'/'+reg+'/wqq'])
    tgo = covnTGA(fd[shapetype+'/'+reg+'/data'])
    
    ax.errorbar(tgo[0], tgo[1], yerr=tgo[2], xerr=tgo[3], fmt='o', color='black', label='Data')
    hep.histplot([ho1[0], ho2[0]], ho1[1], stack=True, ax=ax, label=['Unmatched', "Matched"], histtype='fill')

    ax.set_xlim(40, 136.5)
    ax.set_ylabel('Events', y=1, ha='right')
    ax.legend()

    from uncertainties import unumpy
    data = unumpy.uarray(tgo[1], tgo[2])
    mc1 =  unumpy.uarray(ho1[0], ho1[2])
    mc2 =  unumpy.uarray(ho2[0], ho2[2])
    ratio = data/(mc1 + mc2)

    rax.errorbar(tgo[0], unumpy.nominal_values(ratio)[0], unumpy.std_devs(ratio), xerr=tgo[3], fmt='o', color='black', label='Data')
    rax.hlines(1, 40, 150, linestyle='--', color='k', alpha=0.7)
    rax.set_ylim(0.4, 1.6)
    rax.set_ylabel("Data/MC")
    rax.set_xlabel("Jet $m_{SD}$", x=1, ha='right');
    
    hep.cms.label(ax=ax, data=True, year=args.year, lumi=lumi[int(args.year)])
    
    if i == len(regions) - 1:
        ax.set_ylim(None, ax.get_ylim()[-1]*1.5)    
        if "Secondary" in reg:
            sfl1, sfl2 = args.sflabels.split(",")
            sfstr = ("SF ({}) = {:.3f} $\pm$ {:.3f}".format(sfl1, out['effwSF']['val'], out['effwSF']['unc'])
                + "\nSF ({}) = {:.3f} $\pm$ {:.3f}".format(sfl2, out['effSF']['val'], out['effSF']['unc'])
                + "\nScale = {:.3f} $\pm$ {:.3f}".format(out['shift_SF'], out['shift_SF_ERR'])
                + "\nSmear = {:.3f} $\pm$ {:.3f}".format(out['smear_SF'], out['smear_SF_ERR'])
            )
        else:
            sfl1, sfl2 = args.sflabels.split(",")
            sfstr = (
                "SF ({}) = {:.3f} $\pm$ {:.3f}".format(sfl1, out['effSF']['val'], out['effSF']['unc'])
                + "\nScale = {:.3f} $\pm$ {:.3f}".format(out['shift_SF'], out['shift_SF_ERR'])
                + "\nSmear = {:.3f} $\pm$ {:.3f}".format(out['smear_SF'], out['smear_SF_ERR'])
            )
            
        at = AnchoredText(sfstr,
                    loc='upper left', frameon=False, prop=dict(size=20)
                    )
        ax.add_artist(at)

    fig.savefig('{}/{}.{}'.format(args.output_folder, shapetype+'_'+reg, 'pdf'),
                bbox_inches="tight")
    fig.savefig('{}/{}.{}'.format(args.output_folder, shapetype+'_'+reg, 'png'),
                bbox_inches="tight")
