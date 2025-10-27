mixcr buildLibrary --debug \
    --v-genes-from-fasta HomoSapiens_IGHV.fasta \
    --v-gene-feature VRegion \
    --j-genes-from-fasta HomoSapiens_IGHJ.fasta \
    --d-genes-from-fasta HomoSapiens_IGHD.fasta \
    --c-genes-from-fasta HomoSapiens_AB_C.fasta \
    --chain IGH \
    --taxon-id 9606 \
    --species hsapiens \
    IGH.json.gz

mixcr buildLibrary --debug \
    --v-genes-from-fasta HomoSapiens_IGKV.fasta \
    --v-gene-feature VRegion \
    --c-genes-from-fasta HomoSapiens_AB_C.fasta \
    --j-genes-from-fasta HomoSapiens_IGKJ.fasta \
    --chain IGK \
    --taxon-id 9606 \
    --species hsapiens \
    IGK.json.gz

mixcr buildLibrary --debug \
    --v-genes-from-fasta HomoSapiens_IGLV.fasta \
    --v-gene-feature VRegion \
    --j-genes-from-fasta HomoSapiens_IGLJ.fasta \
    --c-genes-from-fasta HomoSapiens_AB_C.fasta \
    --chain IGL \
    --taxon-id 9606 \
    --species hsapiens \
    IGL.json.gz

mixcr buildLibrary --debug \
    --v-genes-from-fasta HomoSapiens_TRAV.fasta \
    --v-gene-feature VRegion \
    --j-genes-from-fasta HomoSapiens_TRAJ.fasta \
    --c-genes-from-fasta HomoSapiens_AB_C.fasta \
    --chain TRA \
    --taxon-id 9606 \
    --species hsapiens \
    TRA.json.gz

mixcr buildLibrary --debug \
    --v-genes-from-fasta HomoSapiens_TRGV.fasta \
    --v-gene-feature VRegion \
    --j-genes-from-fasta HomoSapiens_TRGJ.fasta \
    --chain TRG \
    --taxon-id 9606 \
    --species hsapiens \
    TRG.json.gz

mixcr buildLibrary --debug \
    --v-genes-from-fasta HomoSapiens_TRBV.fasta \
    --v-gene-feature VRegion \
    --j-genes-from-fasta HomoSapiens_TRBJ.fasta \
    --c-genes-from-fasta HomoSapiens_AB_C.fasta \
    --d-genes-from-fasta HomoSapiens_TRBD.fasta \
    --chain TRB \
    --taxon-id 9606 \
    --species hsapiens \
    TRB.json.gz

mixcr buildLibrary --debug \
    --v-genes-from-fasta HomoSapiens_TRDV.fasta \
    --v-gene-feature VRegion \
    --c-genes-from-fasta HomoSapiens_AB_C.fasta \
    --j-genes-from-fasta HomoSapiens_TRDJ.fasta \
    --d-genes-from-fasta HomoSapiens_TRDD.fasta \
    --chain TRD \
    --taxon-id 9606 \
    --species hsapiens \
    TRD.json.gz

mixcr mergeLibrary \
    IGH.json.gz \
    IGK.json.gz \
    IGL.json.gz \
    TRA.json.gz \
    TRB.json.gz \
    TRG.json.gz \
    TRD.json.gz \
    bd_lib.json.gz
