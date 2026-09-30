# Toxicity model interpretation

This is an exact decision-tree-path decomposition of the two archived random forests.
Each prediction equals the average root probability plus feature contributions along
the selected paths. The residual is checked below 1e-12. This is path-dependent,
not SHAP, a counterfactual intervention, or evidence of a biological mechanism.
Correlated physicochemical/composition/motif features can trade attribution.
The forests differ in training data and tree count (700 versus 900), so the
comparison describes the fitted models, not an isolated controlled causal effect.
JNP anchors were used for adaptation; this analysis is in-sample interpretation.

## REI-26

Probability 0.296644 -> 0.105449; baseline change -0.010314.

| Feature | Public path contribution | Adapted path contribution | Change |
|---|---:|---:|---:|
| dipep_IS | +0.034961 | +0.004701 | -0.030260 |
| length | +0.033781 | +0.009377 | -0.024404 |
| dipep_GR | +0.013354 | +0.002724 | -0.010630 |
| basic_frac | -0.015407 | -0.025897 | -0.010490 |
| dipep_LP | +0.019771 | +0.009602 | -0.010170 |
| aa_E | -0.015626 | -0.005777 | +0.009849 |
| dipep_KR | -0.003184 | -0.012107 | -0.008923 |
| boman | -0.043728 | -0.035328 | +0.008400 |

Negative changes lower the fitted toxicity probability; positive changes raise it.
Only the largest changes are shown; the complete CSV includes all 432 features.

## TRR-26

Probability 0.266620 -> 0.104980; baseline change -0.010314.

| Feature | Public path contribution | Adapted path contribution | Change |
|---|---:|---:|---:|
| length | +0.036013 | +0.012754 | -0.023258 |
| boman | -0.040491 | -0.022869 | +0.017622 |
| hmoment | +0.017936 | +0.002602 | -0.015334 |
| dipep_KR | -0.001042 | -0.014374 | -0.013332 |
| dipep_IL | +0.011624 | +0.000832 | -0.010791 |
| dipep_AL | +0.019129 | +0.008910 | -0.010219 |
| basic_frac | -0.014007 | -0.024040 | -0.010033 |
| dipep_RT | +0.001451 | -0.008006 | -0.009458 |

Negative changes lower the fitted toxicity probability; positive changes raise it.
Only the largest changes are shown; the complete CSV includes all 432 features.

## LRS-21

Probability 0.226964 -> 0.081752; baseline change -0.010314.

| Feature | Public path contribution | Adapted path contribution | Change |
|---|---:|---:|---:|
| dipep_IS | +0.031297 | +0.005971 | -0.025326 |
| dipep_SP | +0.004480 | -0.005554 | -0.010033 |
| dipep_SI | +0.015102 | +0.005430 | -0.009672 |
| alpha | -0.000165 | -0.007440 | -0.007276 |
| charge | +0.004646 | -0.002030 | -0.006676 |
| aa_S | +0.004353 | -0.002220 | -0.006573 |
| dipep_SR | +0.001462 | -0.005011 | -0.006474 |
| apv | -0.006552 | -0.012966 | -0.006414 |

Negative changes lower the fitted toxicity probability; positive changes raise it.
Only the largest changes are shown; the complete CSV includes all 432 features.

## GFD-30

Probability 0.261923 -> 0.105537; baseline change -0.010314.

| Feature | Public path contribution | Adapted path contribution | Change |
|---|---:|---:|---:|
| length | +0.035143 | +0.018551 | -0.016592 |
| dipep_GY | +0.018951 | +0.007565 | -0.011386 |
| dipep_SP | +0.002826 | -0.006625 | -0.009450 |
| dipep_KR | +0.000577 | -0.008852 | -0.009429 |
| dipep_IL | +0.008317 | -0.000983 | -0.009299 |
| aa_R | +0.017511 | +0.009490 | -0.008021 |
| basic_frac | -0.005969 | -0.013888 | -0.007919 |
| aa_F | -0.001009 | -0.008514 | -0.007505 |

Negative changes lower the fitted toxicity probability; positive changes raise it.
Only the largest changes are shown; the complete CSV includes all 432 features.

## TFK-18

Probability 0.276969 -> 0.726842; baseline change -0.010314.

| Feature | Public path contribution | Adapted path contribution | Change |
|---|---:|---:|---:|
| aa_R | -0.028645 | -0.003121 | +0.025524 |
| boman | -0.027637 | -0.005117 | +0.022520 |
| hydrophobic_frac | -0.013785 | +0.006572 | +0.020357 |
| dipep_NG | +0.009267 | +0.029338 | +0.020071 |
| dipep_RN | -0.003655 | +0.014705 | +0.018360 |
| dipep_FV | -0.013574 | +0.001570 | +0.015144 |
| basic_frac | +0.000004 | +0.014101 | +0.014097 |
| dipep_FL | +0.019441 | +0.032909 | +0.013468 |

Negative changes lower the fitted toxicity probability; positive changes raise it.
Only the largest changes are shown; the complete CSV includes all 432 features.

## LL-37

Probability 0.137395 -> 0.678480; baseline change -0.010314.

| Feature | Public path contribution | Adapted path contribution | Change |
|---|---:|---:|---:|
| aa_R | -0.023880 | +0.002019 | +0.025899 |
| boman | -0.035601 | -0.010520 | +0.025081 |
| hydrophobic_frac | -0.024271 | -0.002187 | +0.022085 |
| aa_E | -0.019135 | +0.002434 | +0.021568 |
| dipep_FL | +0.014519 | +0.031979 | +0.017459 |
| dipep_DF | +0.003447 | +0.020565 | +0.017118 |
| dipep_KE | -0.005484 | +0.011266 | +0.016751 |
| dipep_EF | +0.000913 | +0.014810 | +0.013897 |

Negative changes lower the fitted toxicity probability; positive changes raise it.
Only the largest changes are shown; the complete CSV includes all 432 features.

## melittin

Probability 0.920571 -> 0.962324; baseline change -0.010314.

| Feature | Public path contribution | Adapted path contribution | Change |
|---|---:|---:|---:|
| dipep_IS | +0.041020 | +0.018816 | -0.022204 |
| dipep_AV | +0.008470 | +0.022835 | +0.014365 |
| dipep_LP | +0.029255 | +0.018841 | -0.010414 |
| aa_F | -0.017110 | -0.007478 | +0.009632 |
| dipep_QQ | +0.007085 | +0.016668 | +0.009583 |
| dipep_TT | +0.013745 | +0.022659 | +0.008914 |
| dipep_KV | +0.007559 | +0.016318 | +0.008759 |
| length | +0.030227 | +0.022088 | -0.008139 |

Negative changes lower the fitted toxicity probability; positive changes raise it.
Only the largest changes are shown; the complete CSV includes all 432 features.

TFK-18 refers to the preserved 33-residue archive record; see the sequence caveat in JNP_DOMAIN_ADAPTATION.md.
