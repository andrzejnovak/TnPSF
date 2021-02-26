from __future__ import print_function, division
import sys
import os
import rhalphalib as rl
import numpy as np
import scipy.stats
import pickle
import uproot


def get_templ(f, sample, syst=None, sumw2=True):
    hist_name = sample
    if syst is not None:
        hist_name += "_" + syst
    h_vals = f[hist_name].values
    h_edges = f[hist_name].edges
    h_variances = f[hist_name].variances
    h_key = 'msd'
    if not sumw2:
        return (h_vals, h_edges, h_key)
    else:
        return (h_vals, h_edges, h_key, h_variances)


def test_sfmodel(tmpdir, fittype='single', scale=1, smear=0.1, template_dir={}):
    sys_scale = rl.NuisanceParameter('CMS_scale', 'shape')
    sys_smear = rl.NuisanceParameter('CMS_smear', 'shape')
    lumi = rl.NuisanceParameter('CMS_lumi', 'lnN')
    jecs = rl.NuisanceParameter('CMS_jecs', 'lnN')
    pu = rl.NuisanceParameter('CMS_pu', 'lnN')
    effSF = rl.IndependentParameter('effSF', 1., -20, 100)
    effSF2 = rl.IndependentParameter('effSF_un', 1., -20, 100)
    effwSF = rl.IndependentParameter('effwSF', 1., -20, 100)
    effwSF2 = rl.IndependentParameter('effwSF_un', 1., -20, 100)

    msdbins = np.linspace(40, 201, 24)
    msd = rl.Observable('msd', msdbins)
    model = rl.Model("sfModel")

    if fittype == 'single':
        regions = ['pass', 'fail']
    elif fittype == 'double':
        regions = ['pass', 'passfail', 'fail']
    else:
        raise NotImplementedError

    for region in regions:
        ch = rl.Channel("wsf{}".format(region))
        wqq_templ = get_templ(template_dir[region], 'catp2')
        wqq_sample = rl.TemplateSample("{}_wqq".format(ch.name), rl.Sample.SIGNAL, wqq_templ)
        wqq_sample.setParamEffect(sys_scale, 
                                  get_templ(template_dir[region], 'catp2', 'scaleUp', False),
                                  get_templ(template_dir[region], 'catp2', 'scaleDown', False),
                                  scale=scale,
                                  )
        wqq_sample.setParamEffect(sys_smear, 
                                  get_templ(template_dir[region], 'catp2', 'smearUp', False),
                                  get_templ(template_dir[region], 'catp2', 'smearDown', False),
                                  scale=smear,
                                  )
        wqq_sample.setParamEffect(lumi, 1.023)
        wqq_sample.setParamEffect(jecs, 1.02)
        wqq_sample.setParamEffect(pu, 1.05)
        wqq_sample.autoMCStats(lnN=True)
        ch.addSample(wqq_sample)

        qcd_templ = get_templ(template_dir[region], 'catp1')
        qcd_sample = rl.TemplateSample("{}_qcd".format(ch.name), rl.Sample.BACKGROUND, qcd_templ)
        qcd_sample.setParamEffect(lumi, 1.023)
        qcd_sample.setParamEffect(jecs, 1.02)
        qcd_sample.setParamEffect(pu, 1.05)
        qcd_sample.autoMCStats(lnN=True)
        ch.addSample(qcd_sample)

        data_obs = get_templ(template_dir[region], 'data_obs')[:-1]
        #mask = data_obs[0] == 0
        #print(mask)
        #ch.mask = mask
        ch.setObservation(data_obs)

        model.addChannel(ch)

    if fittype == 'single':
        for sample, SF in zip(['wqq', 'qcd'], [effSF, effSF2]):
            pass_sample1 = model['wsfpass'][sample]
            fail_sample = model['wsffail'][sample]
            pass_fail = pass_sample1.getExpectation(nominal=True).sum() / fail_sample.getExpectation(nominal=True).sum()
            pass_sample1.setParamEffect(SF, 1.0 * SF)
            fail_sample.setParamEffect(SF, (1 - SF) * pass_fail + 1)

    elif fittype == 'double':
        for sample, SF in zip(['wqq', 'qcd'], [effSF, effSF2]):
            pass_sample1 = model['wsfpass'][sample]
            pass_sample2 = model['wsfpassfail'][sample]
            fail_sample = model['wsffail'][sample]
            pass_fail = (pass_sample1.getExpectation(nominal=True).sum() + 
                         pass_sample2.getExpectation(nominal=True).sum()) / fail_sample.getExpectation(nominal=True).sum()
            pass_sample1.setParamEffect(SF, 1.0 * SF)
            pass_sample2.setParamEffect(SF, 1.0 * SF)
            fail_sample.setParamEffect(SF, (1 - SF) * pass_fail + 1)

        for sample, SF in zip(['wqq', 'qcd'], [effwSF, effwSF2]):
            pass_sample = model['wsfpass'][sample]
            fail_sample = model['wsfpassfail'][sample]
            pass_fail = pass_sample.getExpectation(nominal=True).sum() / fail_sample.getExpectation(nominal=True).sum()
            pass_sample.setParamEffect(SF, 1.0 * SF)
            fail_sample.setParamEffect(SF, (1 - SF) * pass_fail + 1)

    model.renderCombine('fitdir')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    def str2bool(v):
        if isinstance(v, bool):
            return v
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    parser.add_argument("--fit", type=str, choices={'single', 'double'}, default='double',
                        help="Fit type")
    parser.add_argument("--scale", type=float, default='1',
                        help="Datacard magnitude for scale shift.")
    parser.add_argument("--smear", type=float, default='0.1',
                        help="Datacard magnitude for smear shift.")

    parser.add_argument("--tp", "--template-pass", dest='tp', type=str, 
                        #default='templates/wfit_nskim17_n2cvb/wtag_var_pass.root',
                        default='templates/wfit_nskim17_cvl/wtag_var_pass.root',
                        help="Pass(Pass/Pass) templates")

    parser.add_argument("--tpf", "--template-passfail", dest='tpf', type=str, 
                        default='templates/wfit_nskim17_cvl/wtag_var_fail.root',
                        help="Pass/Fail templates, only for `fit=double`")  

    parser.add_argument("--tf", "--template-fail", dest='tf', type=str, 
                        default='templates/wfit_nskim17_n2cvb/wtag_var_fail.root',
                        help="Fail templates")

    args = parser.parse_args()
    print("Running with options:")
    print("    ", args)

    if not os.path.exists('fitdir'):
        os.mkdir('fitdir')

    if args.fit == 'single':
        regions = ['pass', 'fail']
    elif args.fit == 'double':
        regions = ['pass', 'passfail', 'fail']
    else:
        raise NotImplementedError

    templates = {}    
    for region, path in zip(regions, [args.tp, args.tpf, args.tf,]):
        try:
            templates[region] = uproot.open(path)
        except:
            raise ValueError("Expected a root file at", path )

    test_sfmodel('fitdir', 
        fittype=args.fit,
        scale=args.scale,
        smear=args.smear,
        template_dir=templates,
        )
