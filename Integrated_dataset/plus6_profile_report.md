# plus6 integrated milestone report

## Overview

- Final plus6 cells: `5,128,904`
- Final plus6 GSEs: `31`
- Final plus6 Leiden clusters: `25`

## Input compatibility

| alias                   | source_gse_id          | n_cells | n_genes_before_align | n_samples | n_paired_ab | n_paired_gd | sorted_gdt_flag | count_matrix_source                    | n_genes_after_align | n_missing_vs_base | output_h5ad                                                                                                                                                   |
| ----------------------- | ---------------------- | ------- | -------------------- | --------- | ----------- | ----------- | --------------- | -------------------------------------- | ------------------- | ----------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| GDT_2020AUG_woCOV       | GDT_2020AUG_woCOV      | 25904   | 14122                | 4         | 0           | 15066       | True            | X_integer_like                         | 27413               | 13478             | /home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/high_speed_temp/Integrated_dataset/plus6_inputs/GDT_2020AUG_woCOV_plus6_prepped.h5ad       |
| GDTlung2023july_7p      | GDTlung2023july_7p     | 15175   | 33538                | 14        | 0           | 2927        | True            | X_integer_like                         | 27413               | 7123              | /home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/high_speed_temp/Integrated_dataset/plus6_inputs/GDTlung2023july_7p_plus6_prepped.h5ad      |
| GSE144469               | GSE144469              | 107068  | 19109                | 39        | 65616       | 10556       | False           | layers[counts]                         | 27413               | 9169              | /home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/high_speed_temp/Integrated_dataset/plus6_inputs/GSE144469_plus6_prepped.h5ad               |
| GSE206325               | GSE206325              | 501012  | 13303                | 288       | 0           | 0           | False           | layers[counts]                         | 27413               | 14464             | /home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/high_speed_temp/Integrated_dataset/plus6_inputs/GSE206325_plus6_prepped.h5ad               |
| HRA005041               | HRA005041              | 766639  | 32602                | 321       | 547970      | 5471        | False           | pseudo_counts_from_seurat_lognormalize | 27413               | 7215              | /home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/high_speed_temp/Integrated_dataset/plus6_inputs/HRA005041_plus6_prepped.h5ad               |
| MalteGDT                | MalteGDT               | 7800    | 12655                | 2         | 0           | 5650        | True            | X_integer_like                         | 27413               | 14838             | /home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/high_speed_temp/Integrated_dataset/plus6_inputs/MalteGDT_plus6_prepped.h5ad                |
| current_integrated_base | multiple_existing_gses | 3705306 | 27413                | 767       | 651958      | 0           | False           | X                                      | 27413               | 0                 | /home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/high_speed_temp/Integrated_dataset/plus6_inputs/current_integrated_base_plus6_prepped.h5ad |

## Phase 3

| n_cells | n_genes | n_batches | n_gses | n_leiden | hvg_method | runtime_seconds   |
| ------- | ------- | --------- | ------ | -------- | ---------- | ----------------- |
| 5128904 | 27413   | 1422      | 31     | 25       | seurat_v3  | 20087.35547542572 |

## Phase 4 score summary

| score_name                   | n_cells | mean                | std                | min                 | median              | p95                 | p99                | max                |
| ---------------------------- | ------- | ------------------- | ------------------ | ------------------- | ------------------- | ------------------- | ------------------ | ------------------ |
| phase4_tra_score             | 5128904 | -0.0515778698027133 | 0.0270633324980735 | -0.1815230697393417 | -0.0526680946350097 | -0.005500815808773  | 0.0164763629436492 | 0.4353005886077881 |
| phase4_trb_score             | 5128904 | -0.0551135949790477 | 0.0369601733982563 | -0.2288609445095062 | -0.0552418753504753 | 0.0052126571536064  | 0.0342838913202285 | 0.3020846247673034 |
| phase4_trab_score            | 5128904 | -0.0464685671031475 | 0.0265924483537673 | -0.1913400739431381 | -0.0462781563401222 | -0.0032932832837104 | 0.0160468555986881 | 0.2511661648750305 |
| phase4_trd_score             | 5128904 | 0.002657011616975   | 0.1302253752946853 | -0.1869708895683288 | -0.0340102314949035 | 0.2840672433376312  | 0.6334869265556335 | 1.5339179039001465 |
| phase4_trd_minus_trab        | 5128904 | 0.0491255708038806  | 0.1377832442522049 | -0.3822924494743347 | 0.0106742661446332  | 0.3480985462665558  | 0.7021002173423767 | 1.592238426208496  |
| phase4_trd_score_scaled      | 5128904 | 0.1101917326450347  | 0.075673297047615  | 0.0                 | 0.0888846814632415  | 0.2737179398536682  | 0.4767639636993408 | 1.0                |
| phase4_trab_score_scaled     | 5128904 | 0.3273885250091553  | 0.0600950941443443 | 0.0                 | 0.3278189301490783  | 0.4249584674835205  | 0.4686644077301025 | 1.0                |
| phase4_trd_minus_trab_scaled | 5128904 | -0.2171967625617981 | 0.1052103713154792 | -0.967548966407776  | -0.2375347167253494 | -0.0060873925685882 | 0.2032447755336761 | 0.7029247879981995 |

## Top GSE composition

| source_gse_id     | n_cells | phase4_trab_score_mean | phase4_trd_score_mean | phase4_trd_minus_trab_mean | phase4_trd_minus_trab_median |
| ----------------- | ------- | ---------------------- | --------------------- | -------------------------- | ---------------------------- |
| HRA005041         | 766639  | -0.051022034           | 0.012618273           | 0.063640304                | 0.009900527                  |
| GSE243013         | 759436  | -0.050323747           | -0.009930654          | 0.04039309                 | 0.015378488                  |
| GSE206325         | 501012  | -0.039145168           | 0.010443118           | 0.049588285                | 0.0054617003                 |
| GSE254249         | 406494  | -0.052943505           | 0.029713355           | 0.08265685                 | 0.017107306                  |
| GSE243905         | 350438  | -0.039652947           | -0.012097839          | 0.02755511                 | 0.0029981947                 |
| GSE287301         | 340227  | -0.057241328           | -0.040121302          | 0.017120024                | 0.005696323                  |
| GSE235863         | 241766  | -0.054777894           | -0.018749679          | 0.036028218                | 0.016731333                  |
| GSE212217         | 213706  | -0.03125158            | 0.012095077           | 0.043346662                | 0.00070753973                |
| GSE243572         | 188217  | -0.05227368            | -0.000679727          | 0.051593956                | 0.030236356                  |
| GSE228597         | 149812  | -0.030700447           | 0.038330663           | 0.06903111                 | 0.009437133                  |
| GSE188620         | 123218  | -0.05278486            | -0.011571252          | 0.04121361                 | 0.020554872                  |
| GSE308075         | 112337  | -0.03196189            | -0.033403225          | -0.0014413347              | -0.008720733                 |
| GSE287541         | 110874  | -0.04169737            | -0.016335618          | 0.025361756                | 0.006719536                  |
| GSE144469         | 107068  | -0.030282063           | 0.0004595502          | 0.030741613                | -0.026478382                 |
| GSE254176         | 96953   | -0.041343518           | 0.0018213745          | 0.043164894                | 0.0047306344                 |
| GSE267645         | 92501   | -0.057526663           | -0.0053662057         | 0.052160457                | 0.02440473                   |
| GSE311112         | 90599   | -0.057418093           | -0.009060316          | 0.04835778                 | 0.027409064                  |
| GSE162498         | 75880   | -0.03311588            | -0.014638874          | 0.018477006                | 0.0020527001                 |
| GSE178882         | 67755   | -0.04918713            | -0.0029068857         | 0.046280246                | 0.026218463                  |
| GSE125527         | 46396   | -0.030399946           | 0.013178277           | 0.043578226                | 0.007789281                  |
| GSE252762         | 45132   | -0.043470643           | 0.042292915           | 0.08576356                 | 0.009451585                  |
| GSE301528         | 44105   | -0.042616334           | 0.011575231           | 0.054191563                | 0.011502199                  |
| GSE211504         | 37244   | -0.02661628            | -0.026478376          | 0.00013790399              | -0.003032012                 |
| GSE171037         | 31937   | -0.026983105           | -0.028015286          | -0.0010321827              | -0.0072620567                |
| GDT_2020AUG_woCOV | 25904   | -0.055878177           | 0.52281415            | 0.5786924                  | 0.59703964                   |

## Tissue composition

| tissue                              | n_cells | phase4_trab_score_mean | phase4_trd_score_mean | phase4_trd_minus_trab_mean | phase4_trd_minus_trab_median |
| ----------------------------------- | ------- | ---------------------- | --------------------- | -------------------------- | ---------------------------- |
| blood                               | 1924475 | -0.045511767           | 0.0008124143          | 0.04632418                 | 0.012410499                  |
| tumor                               | 1145163 | -0.05161832            | -0.019048631          | 0.032569688                | 0.011855906                  |
| Tumor                               | 236892  | -0.033132803           | -0.0067534503         | 0.02637935                 | -0.0021884423                |
| Normal                              | 227320  | -0.044730064           | 0.027651185           | 0.07238125                 | 0.013801876                  |
| heart                               | 169851  | -0.03312401            | 0.032803122           | 0.06592713                 | 0.00936424                   |
| csf                                 | 128340  | -0.030479452           | -0.017576583          | 0.012902871                | -0.004758222                 |
| rectum                              | 124404  | -0.047390576           | -0.011475451          | 0.035915125                | 0.013297651                  |
| duodenum                            | 102207  | -0.049317084           | 0.0511396             | 0.10045668                 | 0.0111657195                 |
| liver_tumor                         | 95702   | -0.053346194           | -0.025974624          | 0.027371572                | 0.017392913                  |
| stomach                             | 93885   | -0.048773434           | -0.0026983062         | 0.046075124                | 0.0065023005                 |
| peripheral blood                    | 83669   | -0.050580118           | 0.113503374           | 0.1640835                  | 0.013465103                  |
| Descending, sigmoidal colon; rectum | 80500   | -0.03095097            | 0.0018986122          | 0.03284958                 | -0.025604423                 |
| jejunum                             | 72733   | -0.05381615            | 0.0717087             | 0.12552485                 | 0.013412576                  |
| pancreas                            | 67484   | -0.04652217            | -0.020285033          | 0.026237138                | 0.002673043                  |
| esophagus                           | 53012   | -0.054410458           | -0.004119788          | 0.05029067                 | 0.013102457                  |
| sigmoid colon                       | 48557   | -0.048652984           | 0.01347121            | 0.06212419                 | 0.007950645                  |
| bone_marrow                         | 44105   | -0.042616334           | 0.011575231           | 0.054191563                | 0.011502199                  |
| spleen                              | 37330   | -0.052676987           | 0.014472248           | 0.06714924                 | 0.009717016                  |
|                                     | 36800   | -0.043349613           | 0.01484507            | 0.058194686                | 0.006016005                  |
| lung                                | 32691   | -0.057779483           | 0.06628624            | 0.12406573                 | 0.024782773                  |
| liver                               | 23706   | -0.05862912            | 0.022624282           | 0.08125341                 | 0.015726887                  |
| ileum                               | 22355   | -0.04532739            | 0.015181027           | 0.060508415                | 0.009998187                  |
| balf                                | 21039   | -0.06362223            | -0.018964438          | 0.044657793                | 0.030014269                  |
| trachea                             | 20936   | -0.058760513           | -0.018248508          | 0.040512003                | 0.022668388                  |
| mesenteric lymph node               | 20733   | -0.045513444           | -0.011316866          | 0.034196578                | 0.007119205                  |

## Annotation cluster summary

| leiden | n_cells | n_gses | paired_ab_fraction | paired_gd_fraction | phase4_trd_score_median | phase4_trab_score_median | phase4_trd_minus_trab_median | gdt_marker_mean    | treg_marker_mean   | nk_marker_mean     | tcell_marker_mean  | cd4_marker_mean    | cd8_marker_mean    | gdt_rule | nk_rule | treg_rule | cd8_rule | cd4_rule | treg_cd4_ge_cd8 | treg_nk_low | simple_annotation_plus6 |
| ------ | ------- | ------ | ------------------ | ------------------ | ----------------------- | ------------------------ | ---------------------------- | ------------------ | ------------------ | ------------------ | ------------------ | ------------------ | ------------------ | -------- | ------- | --------- | -------- | -------- | --------------- | ----------- | ----------------------- |
| 15     | 107880  | 29     | 0.0415461624026696 | 0.0008713385242862 | -0.0376454330980777     | -0.0774349197745323      | 0.0397311374545097           | 0.0                | 0.0                | 0.0                | 0.2811846733093261 | 0.6077302098274231 | 0.0                | False    | False   | False     | False    | True     | True            | True        | CD4_T                   |
| 6      | 191877  | 30     | 0.2532820504802556 | 0.0347149475966374 | -0.028009183704853      | -0.0509499497711658      | 0.0201619807630777           | 0.0                | 0.0                | 0.9150388240814208 | 1.128160834312439  | 1.3305789232254028 | 1.2520802021026611 | False    | False   | False     | False    | True     | True            | False       | CD4_T                   |
| 8      | 327887  | 30     | 0.0682796207229929 | 0.000914949357553  | -0.039940271526575      | -0.0505489744246006      | 0.0085182189941406           | 0.0                | 0.0                | 0.0                | 0.2646889984607696 | 0.5293631553649902 | 0.2349841445684433 | False    | False   | False     | False    | True     | True            | True        | CD4_T                   |
| 4      | 354509  | 30     | 0.3735194311004798 | 0.0025020521340784 | -0.0360088795423507     | -0.0440567210316658      | 0.007143847644329            | 0.0                | 0.0                | 0.0                | 1.0958869457244873 | 1.1187341213226318 | 0.5604546070098877 | False    | False   | False     | False    | True     | True            | True        | CD4_T                   |
| 5      | 164505  | 29     | 0.148542597489438  | 0.0002127594905929 | -0.0364428050816059     | -0.0439124256372451      | 0.0065759271383285           | 0.0                | 0.0                | 0.0                | 1.1021606922149658 | 0.915513038635254  | 0.0                | False    | False   | False     | False    | True     | True            | True        | CD4_T                   |
| 7      | 223837  | 30     | 0.520901370193489  | 0.0003440003216626 | -0.0384599231183528     | -0.0450483039021492      | 0.0059772692620754           | 0.0                | 0.0                | 0.1510497480630874 | 1.1071398258209229 | 1.1998845338821411 | 0.9333791732788086 | False    | False   | False     | False    | True     | True            | True        | CD4_T                   |
| 10     | 465413  | 30     | 0.3481896723984933 | 0.00093680236693   | -0.0329542234539985     | -0.0376505926251411      | 0.0042670555412769           | 0.0                | 0.0                | 0.0                | 1.082517385482788  | 1.393742322921753  | 0.0                | False    | False   | False     | False    | True     | True            | True        | CD4_T                   |
| 11     | 615707  | 30     | 0.4926044368506448 | 0.0003004675925399 | -0.0345148183405399     | -0.037157941609621       | 0.0022608982399106           | 0.0                | 0.0                | 0.0                | 1.0888311862945557 | 1.5266317129135132 | 0.0                | False    | False   | False     | False    | True     | True            | True        | CD4_T                   |
| 13     | 223369  | 30     | 0.1051130640330574 | 0.0393608781881102 | 0.1233443990349769      | -0.0696019679307937      | 0.1970373839139938           | 0.0                | 0.0                | 1.5982145071029663 | 0.9029217958450316 | 0.0                | 1.2226529121398926 | False    | False   | False     | True     | False    | False           | False       | CD8_T                   |
| 9      | 453368  | 30     | 0.033696687900337  | 0.0040408674630763 | -0.0218251664191484     | -0.0669590830802917      | 0.0473387762904167           | 0.0                | 0.0                | 3.059014320373535  | 0.8614311218261719 | 0.0                | 1.669079303741455  | False    | False   | False     | True     | False    | False           | False       | CD8_T                   |
| 18     | 72152   | 30     | 0.1052223084599179 | 0.0971144251025612 | -0.0270117111504077     | -0.059513472020626       | 0.0310385972261428           | 0.0                | 0.0                | 1.0673092603683472 | 1.0068585872650146 | 0.0                | 1.1955344676971436 | False    | False   | False     | True     | False    | False           | False       | CD8_T                   |
| 19     | 34596   | 28     | 0.0238177824025898 | 0.0026303618915481 | -0.0301446430385112     | -0.0490427315235137      | 0.0164943244308233           | 0.0                | 0.3038396835327148 | 0.433466762304306  | 1.427048921585083  | 0.6630271673202515 | 1.5095372200012207 | False    | False   | False     | True     | False    | False           | False       | CD8_T                   |
| 2      | 90576   | 30     | 0.2881116410528175 | 0.0069113230877936 | -0.0512070097029209     | -0.0698814839124679      | 0.0162204056978225           | 0.0                | 0.0                | 0.5917186737060547 | 1.3488811254501345 | 0.0                | 0.7720804810523987 | False    | False   | False     | True     | False    | False           | False       | CD8_T                   |
| 17     | 21283   | 30     | 0.0617864022929098 | 0.0005638302870835 | -0.0387365408241748     | -0.0511512830853462      | 0.0103308409452438           | 0.0                | 0.0                | 0.4271779954433441 | 1.0454158782958984 | 0.5175121426582336 | 0.9557891488075256 | False    | False   | False     | True     | False    | False           | False       | CD8_T                   |
| 3      | 530147  | 30     | 0.2044885663787591 | 0.0051872405200821 | -0.0335335806012153     | -0.0420755892992019      | 0.0070427693426609           | 0.0                | 0.0                | 0.6590979695320129 | 1.162418007850647  | 0.3763015568256378 | 1.851546049118042  | False    | False   | False     | True     | False    | False           | False       | CD8_T                   |
| 22     | 68137   | 30     | 0.3056635895328529 | 0.000425613103013  | -0.0337272323668003     | -0.0381954908370971      | 0.0034085325896739           | 0.0                | 0.0                | 1.7217072248458862 | 1.2292749881744385 | 0.485409140586853  | 1.2959725856781006 | False    | False   | False     | True     | False    | False           | False       | CD8_T                   |
| 12     | 304966  | 30     | 0.3013909747316094 | 0.0015772250021313 | -0.0329069197177886     | -0.0340136624872684      | -0.0003819051198661          | 0.0                | 0.0                | 2.129648208618164  | 1.297743558883667  | 0.0                | 2.506847858428955  | False    | False   | False     | True     | False    | False           | False       | CD8_T                   |
| 16     | 262628  | 31     | 0.0614290936229191 | 0.003107056368704  | -0.0400705225765705     | -0.0411427207291126      | -0.0007164422422647          | 0.0                | 0.3627763390541076 | 1.0201274156570437 | 1.2737575769424438 | 0.0                | 2.44394326210022   | False    | False   | False     | True     | False    | False           | False       | CD8_T                   |
| 21     | 25443   | 29     | 0.0723578194395315 | 0.0016114451912117 | -0.0493887923657894     | -0.1010156720876693      | 0.0496662929654121           | 0.0                | 0.0                | 0.4400752484798431 | 0.0                | 0.1318529695272445 | 0.0                | False    | True    | False     | False    | False    | True            | False       | NK_cell                 |
| 26     | 21      | 1      | 0.0                | 0.0                | -0.037295363843441      | -0.0810733810067176      | 0.0442795604467392           | 0.0                | 0.0                | 0.613910973072052  | 0.0                | 0.758529543876648  | 0.3969292044639587 | False    | True    | False     | False    | False    | True            | False       | NK_cell                 |
| 0      | 231266  | 31     | 0.1044079112364117 | 0.0001816090562382 | -0.025305863469839      | -0.0534407347440719      | 0.0288165584206581           | 0.0                | 0.0                | 0.6111088395118713 | 0.0                | 0.0                | 0.2717048823833465 | False    | True    | False     | False    | False    | False           | False       | NK_cell                 |
| 1      | 319540  | 29     | 0.3384020779871065 | 0.0003035613694686 | -0.0397040285170078     | -0.0324412174522876      | -0.0080878566950559          | 0.0                | 0.7245356440544128 | 0.0                | 1.173225998878479  | 1.0262006521224976 | 0.0                | False    | False   | True      | False    | True     | True            | True        | Treg                    |
| 14     | 19705   | 29     | 0.1257041360060898 | 0.4093884800811976 | 0.4813214540481567      | -0.0535392239689827      | 0.5390191674232483           | 0.3345482349395752 | 0.0                | 0.0                | 1.9430385828018188 | 0.8090500831604004 | 0.0                | True     | False   | False     | False    | True     | True            | True        | gdT_cell                |
| 23     | 6136    | 18     | 0.0247718383311603 | 0.0032594524119947 | -0.0585503950715065     | -0.096307024359703       | 0.0333710461854934           | 0.0                | 0.0                | 0.0                | 0.0                | 0.0                | 0.0                | False    | False   | False     | False    | False    | True            | True        | other                   |
| 20     | 13956   | 27     | 0.1905273717397535 | 0.0209229005445686 | -0.0256383065134286     | -0.052821010351181       | 0.0295810122042894           | 0.0                | 0.0                | 0.0                | 0.0                | 0.0                | 0.0                | False    | False   | False     | False    | False    | True            | True        | other                   |

## Annotation by GSE

| source_gse_id      | simple_annotation_plus6 | n_cells |
| ------------------ | ----------------------- | ------- |
| GDT_2020AUG_woCOV  | gdT_cell                | 14877   |
| GDT_2020AUG_woCOV  | CD4_T                   | 9843    |
| GDT_2020AUG_woCOV  | CD8_T                   | 1137    |
| GDT_2020AUG_woCOV  | Treg                    | 28      |
| GDT_2020AUG_woCOV  | NK_cell                 | 17      |
| GDT_2020AUG_woCOV  | other                   | 2       |
| GDTlung2023july_7p | CD8_T                   | 12219   |
| GDTlung2023july_7p | CD4_T                   | 2854    |
| GDTlung2023july_7p | NK_cell                 | 76      |
| GDTlung2023july_7p | Treg                    | 22      |
| GDTlung2023july_7p | gdT_cell                | 3       |
| GDTlung2023july_7p | other                   | 1       |
| GSE125527          | CD4_T                   | 29344   |
| GSE125527          | CD8_T                   | 12966   |
| GSE125527          | Treg                    | 2691    |
| GSE125527          | NK_cell                 | 762     |
| GSE125527          | other                   | 607     |
| GSE125527          | gdT_cell                | 26      |
| GSE144469          | CD8_T                   | 49393   |
| GSE144469          | CD4_T                   | 47775   |
| GSE144469          | Treg                    | 5720    |
| GSE144469          | other                   | 3355    |
| GSE144469          | NK_cell                 | 799     |
| GSE144469          | gdT_cell                | 26      |
| GSE145926          | NK_cell                 | 11052   |
| GSE145926          | CD8_T                   | 8098    |
| GSE145926          | CD4_T                   | 1010    |
| GSE145926          | other                   | 651     |
| GSE145926          | Treg                    | 224     |
| GSE145926          | gdT_cell                | 4       |
| GSE162498          | CD4_T                   | 40197   |
| GSE162498          | CD8_T                   | 18104   |
| GSE162498          | NK_cell                 | 8858    |
| GSE162498          | Treg                    | 6779    |
| GSE162498          | other                   | 1902    |
| GSE162498          | gdT_cell                | 40      |
| GSE171037          | CD4_T                   | 20043   |
| GSE171037          | CD8_T                   | 9196    |
| GSE171037          | NK_cell                 | 1517    |
| GSE171037          | Treg                    | 1055    |
| GSE171037          | gdT_cell                | 95      |
| GSE171037          | other                   | 31      |
| GSE178882          | CD4_T                   | 36230   |
| GSE178882          | CD8_T                   | 24168   |
| GSE178882          | NK_cell                 | 6090    |
| GSE178882          | Treg                    | 1163    |
| GSE178882          | other                   | 79      |
| GSE178882          | gdT_cell                | 25      |
| GSE188620          | CD4_T                   | 51707   |
| GSE188620          | CD8_T                   | 41691   |
| GSE188620          | NK_cell                 | 27863   |
| GSE188620          | Treg                    | 1670    |
| GSE188620          | other                   | 233     |
| GSE188620          | gdT_cell                | 54      |
| GSE190870          | NK_cell                 | 16735   |
| GSE190870          | CD4_T                   | 10      |
| GSE190870          | other                   | 5       |
| GSE190870          | CD8_T                   | 3       |
| GSE206325          | CD8_T                   | 326160  |
| GSE206325          | CD4_T                   | 143715  |

## Annotation by tissue

| tissue_corrected                                           | simple_annotation_plus6 | n_cells |
| ---------------------------------------------------------- | ----------------------- | ------- |
|                                                            | CD8_T                   | 20313   |
|                                                            | CD4_T                   | 14908   |
|                                                            | Treg                    | 1514    |
|                                                            | NK_cell                 | 46      |
|                                                            | gdT_cell                | 11      |
|                                                            | other                   | 8       |
| Ascending, Descending, sigmoidal colon; rectum             | CD4_T                   | 2506    |
| Ascending, Descending, sigmoidal colon; rectum             | CD8_T                   | 2006    |
| Ascending, Descending, sigmoidal colon; rectum             | other                   | 269     |
| Ascending, Descending, sigmoidal colon; rectum             | Treg                    | 169     |
| Ascending, Descending, sigmoidal colon; rectum             | NK_cell                 | 26      |
| Ascending, Descending, sigmoidal colon; rectum             | gdT_cell                | 1       |
| Ascending, transverse, descending, sigmoidal colon; rectum | CD8_T                   | 5287    |
| Ascending, transverse, descending, sigmoidal colon; rectum | CD4_T                   | 2546    |
| Ascending, transverse, descending, sigmoidal colon; rectum | Treg                    | 672     |
| Ascending, transverse, descending, sigmoidal colon; rectum | other                   | 271     |
| Ascending, transverse, descending, sigmoidal colon; rectum | NK_cell                 | 36      |
| Ascending, transverse, descending, sigmoidal colon; rectum | gdT_cell                | 4       |
| Descending, sigmoidal colon; rectum                        | CD4_T                   | 38021   |
| Descending, sigmoidal colon; rectum                        | CD8_T                   | 34975   |
| Descending, sigmoidal colon; rectum                        | Treg                    | 4199    |
| Descending, sigmoidal colon; rectum                        | other                   | 2648    |
| Descending, sigmoidal colon; rectum                        | NK_cell                 | 640     |
| Descending, sigmoidal colon; rectum                        | gdT_cell                | 17      |
| Normal                                                     | CD8_T                   | 165255  |
| Normal                                                     | CD4_T                   | 57668   |
| Normal                                                     | Treg                    | 3471    |
| Normal                                                     | NK_cell                 | 858     |
| Normal                                                     | other                   | 38      |
| Normal                                                     | gdT_cell                | 30      |
| Sigmoidal colon                                            | CD4_T                   | 1822    |
| Sigmoidal colon                                            | CD8_T                   | 946     |
| Sigmoidal colon                                            | Treg                    | 301     |
| Sigmoidal colon                                            | other                   | 6       |
| Sigmoidal colon                                            | NK_cell                 | 2       |
| Sigmoidal colon                                            | gdT_cell                | 1       |
| Sigmoidal colon, rectum                                    | CD8_T                   | 6179    |
| Sigmoidal colon, rectum                                    | CD4_T                   | 2880    |
| Sigmoidal colon, rectum                                    | Treg                    | 379     |
| Sigmoidal colon, rectum                                    | other                   | 161     |
| Sigmoidal colon, rectum                                    | NK_cell                 | 95      |
| Sigmoidal colon, rectum                                    | gdT_cell                | 3       |
| Tumor                                                      | CD8_T                   | 140592  |
| Tumor                                                      | CD4_T                   | 71139   |
| Tumor                                                      | Treg                    | 23825   |
| Tumor                                                      | NK_cell                 | 1097    |
| Tumor                                                      | other                   | 174     |
| Tumor                                                      | gdT_cell                | 65      |
| adjacent_normal_liver                                      | CD8_T                   | 5509    |
| adjacent_normal_liver                                      | CD4_T                   | 763     |
| adjacent_normal_liver                                      | Treg                    | 9       |
| adjacent_normal_liver                                      | NK_cell                 | 6       |
| adjacent_normal_liver                                      | other                   | 1       |
| adjacent_normal_tissue                                     | CD4_T                   | 1684    |
| adjacent_normal_tissue                                     | CD8_T                   | 779     |
| adjacent_normal_tissue                                     | NK_cell                 | 332     |
| adjacent_normal_tissue                                     | Treg                    | 184     |
| adjacent_normal_tissue                                     | other                   | 27      |
| adjacent_normal_tissue                                     | gdT_cell                | 2       |
| ascending colon                                            | CD4_T                   | 11172   |

## γδ-focused candidate statistics

| metric                                             | value   |
| -------------------------------------------------- | ------- |
| total_cells                                        | 5128904 |
| sorted_gdT_true                                    | 48879   |
| paired_TRG_TRD_not_doublet                         | 33987   |
| TRD_minus_TRAB_gt_0p4_no_productive_TRA_TRB_not_NK | 182796  |
| doublet_proxy_paired_ab_and_paired_gd              | 5683    |
| union_any_of_three_criteria                        | 206524  |
| gdT_cell_annotation_total                          | 19705   |
| gdT_cell_annotation_with_paired_TRG_TRD            | 8067    |
| gdT_cell_annotation_meeting_any_three              | 15448   |

## γδ-focused overlap breakdown

| category                                      | n_cells |
| --------------------------------------------- | ------- |
| sorted_gdT_only                               | 15840   |
| paired_TRG_TRD_not_doublet_only               | 2970    |
| TRD_minus_TRAB_gt_0p4_only                    | 147301  |
| sorted_gdT_and_paired_TRG_TRD_only            | 4918    |
| sorted_gdT_and_TRD_minus_TRAB_gt_0p4_only     | 9396    |
| paired_TRG_TRD_and_TRD_minus_TRAB_gt_0p4_only | 7374    |
| all_three_criteria                            | 18725   |

## gdT cells with paired gdTCR by tissue

| tissue                                                     | gdt_cells | gdt_with_paired_TRG_TRD | fraction_paired_TRG_TRD_within_gdt |
| ---------------------------------------------------------- | --------- | ----------------------- | ---------------------------------- |
| cord blood                                                 | 10397     | 5129                    | 0.4933153794363759                 |
| peripheral blood                                           | 4772      | 2929                    | 0.613788767812238                  |
| Descending, sigmoidal colon; rectum                        | 17        | 4                       | 0.2352941176470588                 |
| pancreas                                                   | 58        | 3                       | 0.0517241379310344                 |
| duodenum                                                   | 95        | 1                       | 0.0105263157894736                 |
| lung                                                       | 22        | 1                       | 0.0454545454545454                 |
| thymus                                                     | 2127      | 0                       | 0.0                                |
| blood                                                      | 978       | 0                       | 0.0                                |
| tumor                                                      | 351       | 0                       | 0.0                                |
| csf                                                        | 307       | 0                       | 0.0                                |
| heart                                                      | 115       | 0                       | 0.0                                |
| spleen                                                     | 89        | 0                       | 0.0                                |
| Tumor                                                      | 65        | 0                       | 0.0                                |
| jejunum                                                    | 40        | 0                       | 0.0                                |
| bone marrow                                                | 36        | 0                       | 0.0                                |
| Normal                                                     | 30        | 0                       | 0.0                                |
| stomach                                                    | 29        | 0                       | 0.0                                |
| esophagus                                                  | 28        | 0                       | 0.0                                |
| bone_marrow                                                | 20        | 0                       | 0.0                                |
| sigmoid colon                                              | 20        | 0                       | 0.0                                |
| kidney                                                     | 15        | 0                       | 0.0                                |
| rectum                                                     | 15        | 0                       | 0.0                                |
| liver                                                      | 14        | 0                       | 0.0                                |
| skin                                                       | 11        | 0                       | 0.0                                |
| unknown                                                    | 11        | 0                       | 0.0                                |
| liver_tumor                                                | 6         | 0                       | 0.0                                |
| muscle                                                     | 5         | 0                       | 0.0                                |
| Ascending, transverse, descending, sigmoidal colon; rectum | 4         | 0                       | 0.0                                |
| balf                                                       | 4         | 0                       | 0.0                                |
| mixed_duodenum_pbmc                                        | 4         | 0                       | 0.0                                |
| Sigmoidal colon, rectum                                    | 3         | 0                       | 0.0                                |
| bladder                                                    | 3         | 0                       | 0.0                                |
| mesenteric lymph node                                      | 3         | 0                       | 0.0                                |
| adjacent_normal_tissue                                     | 2         | 0                       | 0.0                                |
| lymph node                                                 | 2         | 0                       | 0.0                                |
| tracheal_aspirate                                          | 2         | 0                       | 0.0                                |
| Ascending, Descending, sigmoidal colon; rectum             | 1         | 0                       | 0.0                                |
| Sigmoidal colon                                            | 1         | 0                       | 0.0                                |
| ascending colon                                            | 1         | 0                       | 0.0                                |
| ileum                                                      | 1         | 0                       | 0.0                                |

## gdT cells meeting at least one of the three criteria by tissue

| tissue                                                     | gdt_cells | gdt_meeting_any_three | fraction_meeting_any_three_within_gdt |
| ---------------------------------------------------------- | --------- | --------------------- | ------------------------------------- |
| cord blood                                                 | 10397     | 10397                 | 1.0                                   |
| peripheral blood                                           | 4772      | 4582                  | 0.960184409052808                     |
| blood                                                      | 978       | 189                   | 0.1932515337423312                    |
| tumor                                                      | 351       | 36                    | 0.1025641025641025                    |
| thymus                                                     | 2127      | 31                    | 0.0145745181006111                    |
| duodenum                                                   | 95        | 29                    | 0.3052631578947368                    |
| csf                                                        | 307       | 28                    | 0.0912052117263843                    |
| heart                                                      | 115       | 23                    | 0.2                                   |
| bone marrow                                                | 36        | 16                    | 0.4444444444444444                    |
| spleen                                                     | 89        | 15                    | 0.1685393258426966                    |
| pancreas                                                   | 58        | 15                    | 0.2586206896551724                    |
| jejunum                                                    | 40        | 13                    | 0.325                                 |
| Tumor                                                      | 65        | 7                     | 0.1076923076923077                    |
| esophagus                                                  | 28        | 7                     | 0.25                                  |
| bone_marrow                                                | 20        | 7                     | 0.35                                  |
| Normal                                                     | 30        | 6                     | 0.2                                   |
| lung                                                       | 22        | 6                     | 0.2727272727272727                    |
| stomach                                                    | 29        | 5                     | 0.1724137931034483                    |
| skin                                                       | 11        | 5                     | 0.4545454545454545                    |
| sigmoid colon                                              | 20        | 3                     | 0.15                                  |
| kidney                                                     | 15        | 3                     | 0.2                                   |
| rectum                                                     | 15        | 3                     | 0.2                                   |
| liver                                                      | 14        | 3                     | 0.2142857142857142                    |
| muscle                                                     | 5         | 3                     | 0.6                                   |
| bladder                                                    | 3         | 3                     | 1.0                                   |
| Descending, sigmoidal colon; rectum                        | 17        | 2                     | 0.1176470588235294                    |
| unknown                                                    | 11        | 2                     | 0.1818181818181818                    |
| balf                                                       | 4         | 2                     | 0.5                                   |
| mixed_duodenum_pbmc                                        | 4         | 2                     | 0.5                                   |
| lymph node                                                 | 2         | 2                     | 1.0                                   |
| liver_tumor                                                | 6         | 1                     | 0.1666666666666666                    |
| adjacent_normal_tissue                                     | 2         | 1                     | 0.5                                   |
| ileum                                                      | 1         | 1                     | 1.0                                   |
| Ascending, transverse, descending, sigmoidal colon; rectum | 4         | 0                     | 0.0                                   |
| Sigmoidal colon, rectum                                    | 3         | 0                     | 0.0                                   |
| mesenteric lymph node                                      | 3         | 0                     | 0.0                                   |
| tracheal_aspirate                                          | 2         | 0                     | 0.0                                   |
| Ascending, Descending, sigmoidal colon; rectum             | 1         | 0                     | 0.0                                   |
| Sigmoidal colon                                            | 1         | 0                     | 0.0                                   |
| ascending colon                                            | 1         | 0                     | 0.0                                   |

