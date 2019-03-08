#!/bin/bash
# Data files for pptc-pdx-oncoprints
# INPUT files
# Approved hugo symbols file = 2019-02-14-Hugo-Symbols-approved.txt
# clinical file = 2019-02-09-pdx-clinical-final-for-paper.txt
# MAF file = 2019-02-14-allpdx-clean-maf-240.rda
# Expression matrix = 2019-02-14-PPTC_FPKM_matrix_withModelID-244.rda
# focal CN matrix = 2019-02-27-short_cn_matrix_fpkm1.txt
# genelists
# neuroblastoma specific leisons = 2018-12-06-neuroblastoma-specific-lesions.txt
# classifier scores = classifier_scores_with_clinical_and_alterations.tsv
# arm lesions = 2018-10-09-arm-lesions.txt
# read collapsed fusion file = DriverFusions_Collapsed.txt
# renal specific leisons = renal-specific-lesions.txt

cd ~
mkdir -p pptc-pdx-oncoprints/data/
cd ~/pptc-pdx-oncoprints/data/

# 1. Approved Hugo Symbols file
wget --output-document='2019-02-14-Hugo-Symbols-approved.txt' https://ndownloader.figshare.com/files/14460317

# 2. Clinical file
wget --output-document='pptc-pdx-clinical-web.txt' https://ndownloader.figshare.com/files/14508536

# 3. Download DNA MAF file
wget --output-document='2019-02-14-allpdx-clean-maf-240.rda' https://ndownloader.figshare.com/files/14414198

# 4. Download Expression matrix
wget --output-document='2019-02-14-PPTC_FPKM_matrix_withModelID-244.rda' https://ndownloader.figshare.com/articles/7751825/versions/2

# 5. focal CN matrix
wget --output-document='short_cn_matrix_fpkm1.txt' https://ndownloader.figshare.com/files/14545979

# 6. download gene lists
wget --output-document='leukemia-goi-list.txt' https://ndownloader.figshare.com/files/14541185
wget --output-document='osteosarcoma-goi-list.txt' https://ndownloader.figshare.com/files/14541191
wget --output-document='renal-goi-list.txt' https://ndownloader.figshare.com/files/14541197
wget --output-document='brain-goi-list.txt' https://ndownloader.figshare.com/files/14541182
wget --output-document='neuroblastoma-goi-list.txt' https://ndownloader.figshare.com/files/14541188
wget --output-document='rare-goi-list.txt' https://ndownloader.figshare.com/files/14541194
wget --output-document='sarcoma-goi-list.txt' https://ndownloader.figshare.com/files/14541200

# 7. Neuroblastoma specific leisons
wget --output-document='neuroblastoma-specific-lesions.txt' https://ndownloader.figshare.com/files/14545973

# 8. Classifier scores
wget --output-document="classifier_scores_with_clinical_and_alterations.tsv" https://ndownloader.figshare.com/files/14545970

# 9. arm leisons
wget --output-document='arm-lesions.txt' https://ndownloader.figshare.com/files/14545967

# 10. Driverfusions collapsed
wget --output-document='DriverFusions_Collapsed.txt' https://ndownloader.figshare.com/files/14545976

#11. renal specific leisons
wget --output-document='renal-specific-lesions.txt' https://ndownloader.figshare.com/files/14550614
