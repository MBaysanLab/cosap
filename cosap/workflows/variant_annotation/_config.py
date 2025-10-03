class VariantAnnotationConfig:
    """Configuration class for variant annotation mappings and thresholds."""

    VEP_COLUMN_MAPPING = {
        "variant_id": "Uploaded_variation",
        "chr": "#CHROM",
        "pos": "POS",
        "ref": "REF",
        "alt": "ALT",
        "gene_id": "Gene",
        "gene_symbol": "SYMBOL",
        "feature": "Feature",
        "feature_type": "Feature_type",
        "consequence": "Consequence",
        "cdna_position": "cDNA_position",
        "cds_position": "CDS_position",
        "rs_id": "Existing_variation",
        "impact": "IMPACT",
        "hgvsg": "HGVSg",
        "hgvsc": "HGVSc",
        "hgvsp": "HGVSp",
        "vep_pick": "PICK",
        "mane": "MANE",
        "canonical": "CANONICAL",
        "gnomadg_af": "gnomADg_AF",
        "gnomadg_afr_af": "gnomADg_AFR_AF",
        "gnomadg_ami_af": "gnomADg_AMI_AF",
        "gnomadg_amr_af": "gnomADg_AMR_AF",
        "gnomadg_asj_af": "gnomADg_ASJ_AF",
        "gnomadg_eas_af": "gnomADg_EAS_AF",
        "gnomadg_fin_af": "gnomADg_FIN_AF",
        "gnomadg_mid_af": "gnomADg_MID_AF",
        "gnomadg_nfe_af": "gnomADg_NFE_AF",
        "gnomadg_oth_af": "gnomADg_OTH_AF",
        "gnomadg_sas_af": "gnomADg_SAS_AF",
        "aminoacid_change": "Amino_acids",
        "sift_score": "SIFT",
        "polyphen_score": "PolyPhen",
        "clinical_significance": "CLIN_SIG",
        "alphamissense_class": "am_class",
        "alphamissense_pathogenicity": "am_pathogenicity",
        "revel_score": "REVEL",
        "spliceai_pred_dp_ag": "SpliceAI_pred_DP_AG",
        "spliceai_pred_dp_al": "SpliceAI_pred_DP_AL",
        "spliceai_pred_dp_dg": "SpliceAI_pred_DP_DG",
        "spliceai_pred_dp_dl": "SpliceAI_pred_DP_DL",
        "spliceai_pred_ds_ag": "SpliceAI_pred_DS_AG",
        "spliceai_pred_ds_al": "SpliceAI_pred_DS_AL",
        "spliceai_pred_ds_dg": "SpliceAI_pred_DS_DG",
        "spliceai_pred_ds_dl": "SpliceAI_pred_DS_DL",
        "spliceai_pred_symbol": "SpliceAI_pred_SYMBOL",
        "pubmed": "PUBMED",
    }

    SPLICEAI_SCORE_COLUMNS = [
        "spliceai_pred_ds_ag",
        "spliceai_pred_ds_al",
        "spliceai_pred_ds_dg",
        "spliceai_pred_ds_dl",
    ]

    SPLICEAI_THRESHOLDS = {
        "high_impact": 0.8,
        "moderate_impact": 0.5,
        "low_impact": 0.2,
    }

    REVEL_THRESHOLDS = {"pathogenic": 0.75, "likely_pathogenic": 0.5, "uncertain": 0.25}
