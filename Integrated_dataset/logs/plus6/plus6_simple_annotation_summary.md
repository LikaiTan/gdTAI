# plus6 simple annotation summary

- gdT markers used in decision logic: `TRDC, TRDV1, TRDV2, TRGV9`
- gdT markers excluded from decision logic: `TRGC1`, `TRGC2`
- clusters labelled: `25`
- changed clusters after Treg/CD8 correction: `2`
- Treg clusters violating `cd4_marker_mean >= cd8_marker_mean`: `0`
- Treg clusters violating `nk_marker_mean < 0.35`: `0`

## Changed clusters

- Leiden 16: `Treg` -> `CD8_T` | n=262,628 | treg=0.363 | cd4=0.000 | cd8=2.444 | nk=1.020
- Leiden 19: `Treg` -> `CD8_T` | n=34,596 | treg=0.304 | cd4=0.663 | cd8=1.510 | nk=0.433

## Cluster labels

- Leiden 15: `CD4_T` | n=107,880 | paired_gd=0.001 | paired_ab=0.042 | TRD-TRAB median=0.040 | treg=0.000 | cd4=0.608 | cd8=0.000 | nk=0.000
- Leiden 6: `CD4_T` | n=191,877 | paired_gd=0.035 | paired_ab=0.253 | TRD-TRAB median=0.020 | treg=0.000 | cd4=1.331 | cd8=1.252 | nk=0.915
- Leiden 8: `CD4_T` | n=327,887 | paired_gd=0.001 | paired_ab=0.068 | TRD-TRAB median=0.009 | treg=0.000 | cd4=0.529 | cd8=0.235 | nk=0.000
- Leiden 4: `CD4_T` | n=354,509 | paired_gd=0.003 | paired_ab=0.374 | TRD-TRAB median=0.007 | treg=0.000 | cd4=1.119 | cd8=0.560 | nk=0.000
- Leiden 5: `CD4_T` | n=164,505 | paired_gd=0.000 | paired_ab=0.149 | TRD-TRAB median=0.007 | treg=0.000 | cd4=0.916 | cd8=0.000 | nk=0.000
- Leiden 7: `CD4_T` | n=223,837 | paired_gd=0.000 | paired_ab=0.521 | TRD-TRAB median=0.006 | treg=0.000 | cd4=1.200 | cd8=0.933 | nk=0.151
- Leiden 10: `CD4_T` | n=465,413 | paired_gd=0.001 | paired_ab=0.348 | TRD-TRAB median=0.004 | treg=0.000 | cd4=1.394 | cd8=0.000 | nk=0.000
- Leiden 11: `CD4_T` | n=615,707 | paired_gd=0.000 | paired_ab=0.493 | TRD-TRAB median=0.002 | treg=0.000 | cd4=1.527 | cd8=0.000 | nk=0.000
- Leiden 13: `CD8_T` | n=223,369 | paired_gd=0.039 | paired_ab=0.105 | TRD-TRAB median=0.197 | treg=0.000 | cd4=0.000 | cd8=1.223 | nk=1.598
- Leiden 9: `CD8_T` | n=453,368 | paired_gd=0.004 | paired_ab=0.034 | TRD-TRAB median=0.047 | treg=0.000 | cd4=0.000 | cd8=1.669 | nk=3.059
- Leiden 18: `CD8_T` | n=72,152 | paired_gd=0.097 | paired_ab=0.105 | TRD-TRAB median=0.031 | treg=0.000 | cd4=0.000 | cd8=1.196 | nk=1.067
- Leiden 19: `CD8_T` | n=34,596 | paired_gd=0.003 | paired_ab=0.024 | TRD-TRAB median=0.016 | treg=0.304 | cd4=0.663 | cd8=1.510 | nk=0.433
- Leiden 2: `CD8_T` | n=90,576 | paired_gd=0.007 | paired_ab=0.288 | TRD-TRAB median=0.016 | treg=0.000 | cd4=0.000 | cd8=0.772 | nk=0.592
- Leiden 17: `CD8_T` | n=21,283 | paired_gd=0.001 | paired_ab=0.062 | TRD-TRAB median=0.010 | treg=0.000 | cd4=0.518 | cd8=0.956 | nk=0.427
- Leiden 3: `CD8_T` | n=530,147 | paired_gd=0.005 | paired_ab=0.204 | TRD-TRAB median=0.007 | treg=0.000 | cd4=0.376 | cd8=1.852 | nk=0.659
- Leiden 22: `CD8_T` | n=68,137 | paired_gd=0.000 | paired_ab=0.306 | TRD-TRAB median=0.003 | treg=0.000 | cd4=0.485 | cd8=1.296 | nk=1.722
- Leiden 12: `CD8_T` | n=304,966 | paired_gd=0.002 | paired_ab=0.301 | TRD-TRAB median=-0.000 | treg=0.000 | cd4=0.000 | cd8=2.507 | nk=2.130
- Leiden 16: `CD8_T` | n=262,628 | paired_gd=0.003 | paired_ab=0.061 | TRD-TRAB median=-0.001 | treg=0.363 | cd4=0.000 | cd8=2.444 | nk=1.020
- Leiden 21: `NK_cell` | n=25,443 | paired_gd=0.002 | paired_ab=0.072 | TRD-TRAB median=0.050 | treg=0.000 | cd4=0.132 | cd8=0.000 | nk=0.440
- Leiden 26: `NK_cell` | n=21 | paired_gd=0.000 | paired_ab=0.000 | TRD-TRAB median=0.044 | treg=0.000 | cd4=0.759 | cd8=0.397 | nk=0.614
- Leiden 0: `NK_cell` | n=231,266 | paired_gd=0.000 | paired_ab=0.104 | TRD-TRAB median=0.029 | treg=0.000 | cd4=0.000 | cd8=0.272 | nk=0.611
- Leiden 1: `Treg` | n=319,540 | paired_gd=0.000 | paired_ab=0.338 | TRD-TRAB median=-0.008 | treg=0.725 | cd4=1.026 | cd8=0.000 | nk=0.000
- Leiden 14: `gdT_cell` | n=19,705 | paired_gd=0.409 | paired_ab=0.126 | TRD-TRAB median=0.539 | treg=0.000 | cd4=0.809 | cd8=0.000 | nk=0.000
- Leiden 23: `other` | n=6,136 | paired_gd=0.003 | paired_ab=0.025 | TRD-TRAB median=0.033 | treg=0.000 | cd4=0.000 | cd8=0.000 | nk=0.000
- Leiden 20: `other` | n=13,956 | paired_gd=0.021 | paired_ab=0.191 | TRD-TRAB median=0.030 | treg=0.000 | cd4=0.000 | cd8=0.000 | nk=0.000
