
GOI = ["psbA2", "rnpB", "pilA1", "psaA", "cmpC", "pixG", "pilB", "pixG;  pisG;  taxP1;  rer1", "kpsM"]

# This maps file names to the corresponding DESeq2 design formula
DESIGNS = {
    "rf": ("~ replicate + fraction", "fraction"),
}

# Just a condition map for nicer Plots
CONDIMAP = {
    "M": "Membrane",
    "TC": "Total Cell",
    "F": "no Puromycin",
    "T": "Pyromycin",
    "C": "Cytosol",
    "L": "Light",
    "D": "Dark"
}

rule all:
    input:
        ld = expand(
            "PipelineData/Plots/Vulcano/set{mode}_design_{design}_cond{condition}_vs_base{baseline}.html", zip,
            condition=["M", "M", "M"], baseline=["TC", "TC", "TC"], design=["rf", "rf", "rf"], mode=["Light", "LightDark", "LightDarkLight"]
        ),
        ld3= expand(
            "PipelineData/Plots/Vulcano/set{mode}_design_{design}_cond{condition}_vs_base{baseline}.html",zip,
            condition=["C", "C", "C"],baseline=["TC", "TC", "TC"],design=["rf", "rf", "rf"],mode=["Light", "LightDark",
                                                                                              "LightDarkLight"]
        ),
        ld2 = expand(
            "PipelineData/Plots/Vulcano/set{mode}_design_{design}_cond{condition}_vs_base{baseline}.html",zip,
            condition=["M", "M", "M"],baseline=["C", "C", "C"],design=["rf", "rf", "rf"],mode=["Light", "LightDark",
                                                                                                  "LightDarkLight"]
            ),
        enrich1 = expand(
            "PipelineData/Plots/Enrichment/{enrichment}_set{mode}_design_{design}_cond{condition}_vs_base{baseline}_ud_{updown}.html" ,
            condition=["M"],baseline=["C", "TC"],design=["rf"], mode=["Light", "LightDark", "LightDarkLight"], enrichment=["GO", "KEGG"], updown=["upregulated", "downregulated"]
            ),
        enrich2= expand(
            "PipelineData/Plots/Enrichment/{enrichment}_set{mode}_design_{design}_cond{condition}_vs_base{baseline}_ud_{updown}.html",
                condition=["C"],baseline=["TC"],design=["rf"],mode=["Light", "LightDark", "LightDarkLight"],enrichment=["GO", "KEGG"],updown=["upregulated", "downregulated"]
        ),


rule setup:
    """Installs R dependencies"""
    output:
        install_file = "InstallFinished.txt"
    script: "setup.R"


rule GOSetup:
    input:
        locustag2GO = "20200113_locusTags_GOterms.tsv",
        locustag2Symbol = "locusTagtoGeneName.csv",
        setup = rules.setup.output.install_file
    output:
        finished_file = "PipelineData/IntermediateData/AnnotationBuild.txt",
        annodb = temporary(directory("PipelineData/AnnotationDB"))
    script: "createAnnotationFile.R"



rule joinDataFrames:
    input:
        counts = "Data/All_CDS.csv",
        spike_ins = "Data/All_ERCC.csv"
    output:
        all_counts = "PipelineData/Data/JoinedCounts.csv",
        annotation = "PipelineData/Data/annotation.csv"
    run:
        import pandas as pd
        df = pd.read_csv(input.counts, index_col=0)
        spike_ins = pd.read_csv(input.spike_ins, index_col=0)
        df = pd.concat((df, spike_ins), axis=0)
        cols = []
        data = {
            "name": [],
            "puromycin": [],
            "fraction": [],
            "replicate": [],
            "sample_group": [],
        }
        # Constructing the annotation File for DESeq and renaming columns for the Light experiment
        for col in df.columns:
            name = "_".join(col.replace("-", "minus_").split("_")[1:])
            data["name"].append(name)
            mp, p, frac, rep = name.split("_")
            p = True if mp == "plus" else False
            sgroup = f"{str(p)[0]}{frac}"
            data["puromycin"].append(p)
            data["fraction"].append(frac)
            data["replicate"].append(rep)
            data["sample_group"].append(sgroup)
        df.columns = data["name"]
        df.index.name = None
        df.index = df.index.str.replace("'", "")
        #df = df[~df.index.str.contains("UTR")]
        df.to_csv(output.all_counts, sep="\t")
        annotation = pd.DataFrame(data)
        annotation.to_csv(output.annotation, index=False, sep="\t")



rule clean_light_dark:
    input:
        table = "Data/CDS_LightDark.csv"
    output:
        table = "PipelineData/Data/cleaned_LightDark.csv",
        annotation =  "PipelineData/Data/annotation_LightDark.csv"
    run:
        import pandas as pd
        df = pd.read_csv(input.table, sep=",", index_col=0)
        data = {
            "name": [],
            "puromycin": [],
            "fraction": [],
            "replicate": [],
            "light_dark": [],
            "sample_group": [],
        }
        # Constructing the annotation File for DESeq and renaming columns for the Light Dark experiment

        for col in df.columns:
            name = "_".join(col.split("_")[1:5])
            name = name.replace("plusP", "plus_P")
            name = name.replace("minP", "minus_P")
            data["name"].append(name)
            mp, p, frac, ld, rep = name.split("_")
            rep = f"RD{rep[-1]}"
            p = True if mp == "plus" else False

            sgroup = f"{str(p)[0]}_{frac}_{ld}"
            data["puromycin"].append(p)
            data["fraction"].append(frac)
            data["replicate"].append(rep)
            data["light_dark"].append(ld)
            data["sample_group"].append(sgroup)
        df.columns = data["name"]
        df.index.name = None
        df.index = df.index.str.replace("'","")
        #df = df[~df.index.str.contains("UTR")]
        df.to_csv(output.table, sep="\t")
        annotation = pd.DataFrame(data)
        annotation.to_csv(output.annotation, index=False, sep="\t")




rule dropPuromycinandDark:
    input:
        counts = rules.joinDataFrames.output.all_counts,
        annotation = rules.joinDataFrames.output.annotation,
        counts_ld = rules.clean_light_dark.output.table,
        annotation_ld = rules.clean_light_dark.output.annotation
    output:
        counts = "PipelineData/Data/JoinedCountsNoPD.csv",
        annotation = "PipelineData/Data/annotationNoPD.csv",
        counts_ld = "PipelineData/Data/cleaned_LightDarkNoPD.csv",
        annotation_ld ="PipelineData/Data/annotation_LightDarkNoPD.csv"
    run:
        import pandas as pd
        df = pd.read_csv(input.counts,index_col=0, sep="\t")
        annotation = pd.read_csv(input.annotation,index_col=0, sep="\t")
        df = df.loc[:, annotation["puromycin"] != True]
        annotation = annotation[annotation["puromycin"] != True]
        df.to_csv(output.counts, sep="\t")
        annotation.to_csv(output.annotation, sep="\t")

        df = pd.read_csv(input.counts_ld,index_col=0,sep="\t")
        annotation = pd.read_csv(input.annotation_ld,index_col=0,sep="\t")
        df = df.loc[:, annotation["puromycin"] != True]
        annotation = annotation[annotation["puromycin"] != True]
        df = df.loc[:, annotation["light_dark"] != "D"]
        df = df.loc[:, ~((annotation["replicate"] == "RD2") & (annotation["fraction"] == "M"))]
        annotation = annotation.loc[annotation["light_dark"] != "D"]
        annotation = annotation.loc[~((annotation["replicate"] == "RD2") & (annotation["fraction"] == "M"))]
        df.to_csv(output.counts_ld,sep="\t")
        annotation.to_csv(output.annotation_ld,sep="\t")

rule joinLightDarkLight:
    input:
        ldcounts = rules.dropPuromycinandDark.output.counts_ld,
        ldanno = rules.dropPuromycinandDark.output.annotation_ld,
        lcounts = rules.dropPuromycinandDark.output.counts,
        lanno = rules.dropPuromycinandDark.output.annotation,
    output:
        counts = "PipelineData/Data/LightDarkLightJoinedCountsNoPD.csv",
        annotation = "PipelineData/Data/LightDarkLightannotationNoPD.csv",
    run:
        import pandas as pd
        lddf = pd.read_csv(input.ldcounts,index_col=0,sep="\t")
        ldannotation = pd.read_csv(input.ldanno,index_col=0,sep="\t")
        ldannotation["dataset"] = "LD"
        ldannotation = ldannotation.drop(("light_dark"), axis=1)
        ldf = pd.read_csv(input.lcounts,index_col=0,sep="\t")
        lannotation = pd.read_csv(input.lanno,index_col=0,sep="\t")
        lannotation["dataset"] = "L"
        df = pd.concat((lddf, ldf), axis=1)
        annotation = pd.concat((ldannotation, lannotation), axis=0)
        df = df[~df.index.str.contains("ERCC")]
        df.to_csv(output.counts, sep="\t")
        annotation.to_csv(output.annotation, sep="\t")


rule DESeq2:
    """
    This rule just runs DESeq for a given design formula for the Light experiment.
    It uses spike ins for normalization
    """

    input:
        counts = rules.dropPuromycinandDark.output.counts,
        annotation = rules.dropPuromycinandDark.output.annotation,
        setup = rules.setup.output.install_file
    output:
        pca_data = "PipelineData/IntermediateData/setLight_design_{design}_normed_pca.csv",
        heatmap = "PipelineData/Plots/Correlation/setLight_design_{design}_correltation.pdf",
        dispest = "PipelineData/Plots/Dispersion/setLight_design_{design}_dispersionestimates.png",
        deseq_result = "PipelineData/IntermediateData/setLight_design_{design}_deseq_res_obj.RData"
    params:
        design = lambda wildcards: DESIGNS[wildcards.design][0]
    script:
        "deseq.R"


rule extract_result_Light:
    """
    This basically calls results using a given comparison on the DESeq object
    """
    input:
        deseq_result = "PipelineData/IntermediateData/setLight_design_{design}_deseq_res_obj.RData"
    params:
        design = lambda wildcards: DESIGNS[wildcards.design][1]
    output:
        result_table = "PipelineData/Tables/DESeq/setLight_design_{design}_cond{condition}_vs_base{baseline}.tsv"
    script: "extractDESeqResult.R"



rule least_significant_Light:
    """
    Identifies least significantly different regulated genes of the Light experiment using
    only Cytosol vs Membrane. This set is used as control genes to normalize the LightDark set
    """
    input:
        tcm = "PipelineData/Tables/DESeq/setLight_design_rf_condM_vs_baseC.tsv"
    output:
        least_significant_list = "PipelineData/IntermediateData/least_sig_NormLightDark.tsv"
    run:
        import pandas as pd
        file = input.tcm
        df = pd.read_csv(file, sep="\t", index_col=0)
        df["absChange"] = df["log2FoldChange"].abs()
        subdf = df
        subdf = df[df["absChange"] < 0.1]
        subdf = subdf[subdf["baseMean"] > 100]
        subdf = subdf[subdf["padj"] >= 0.05]
        subdf = subdf.sort_values(by=["absChange", "padj"], ascending=[True, False])
        outdf = subdf[0:100]
        outdf.to_csv(output.least_significant_list, sep="\t")


rule LightDarkDESeq:
    input:
        counts = rules.dropPuromycinandDark.output.counts_ld,
        annotation = rules.dropPuromycinandDark.output.annotation_ld,
        norm_genes = rules.least_significant_Light.output.least_significant_list,
        setup = rules.setup.output.install_file
    output:
        pca_data = "PipelineData/IntermediateData/setLightDark_design_{design}_normed_pca.csv",
        heatmap = "PipelineData/Plots/Correlation/setLightDark_design_{design}_correltation.pdf",
        dispest = "PipelineData/Plots/Dispersion/setLightDark_design_{design}_dispersionestimates.png",
        deseq_result = "PipelineData/IntermediateData/setLightDark_design_{design}_deseq_res_obj.RData"
    params:
        design = lambda wildcards: DESIGNS[wildcards.design][0]
    script:
        "deseqLightDark.R"

rule LightDarkLightDESeq:
    input:
        counts = rules.joinLightDarkLight.output.counts,
        annotation = rules.joinLightDarkLight.output.annotation,
        norm_genes = rules.least_significant_Light.output.least_significant_list,
        setup = rules.setup.output.install_file
    output:
        pca_data = "PipelineData/IntermediateData/setLightDarkLight_design_{design}_normed_pca.csv",
        heatmap = "PipelineData/Plots/Correlation/setLightDarkLight_design_{design}_correltation.pdf",
        dispest = "PipelineData/Plots/Dispersion/setLightDarkLight_design_{design}_dispersionestimates.png",
        condresult = "PipelineData/IntermediateData/setLightDarkLight_design_{design}_deseq_res_obj.RData"
    params:
        design = lambda wildcards: DESIGNS[wildcards.design][0]
    script:
        "deseqLightDark.R"

rule extract_result_LightDark:
    input:
        deseq_result = "PipelineData/IntermediateData/setLightDark_design_{design}_deseq_res_obj.RData",
        setup = rules.setup.output.install_file
    params:
        design = lambda wildcards: DESIGNS[wildcards.design][1]
    output:
        result_table = "PipelineData/Tables/DESeq/setLightDark_design_{design}_cond{condition}_vs_base{baseline}.tsv"
    script: "extractDESeqResult.R"


rule extract_result_LightDarkLight:
    input:
        deseq_result = "PipelineData/IntermediateData/setLightDarkLight_design_{design}_deseq_res_obj.RData",
        setup= rules.setup.output.install_file
    params:
        design = lambda wildcards: DESIGNS[wildcards.design][1]
    output:
        result_table = "PipelineData/Tables/DESeq/setLightDarkLight_design_{design}_cond{condition}_vs{baseline}.tsv"
    script: "extractDESeqResult.R"


rule plotFoldPval:
    input:
        table = "PipelineData/Tables/DESeq/set{mode}_design_{design}_cond{condition}_vs_base{baseline}.tsv",
        tag2name = "locusTagtoGeneName.csv"
    output:
        plot = "PipelineData/Plots/Vulcano/set{mode}_design_{design}_cond{condition}_vs_base{baseline}.html"
    run:
        import pandas as pd
        import plotly.graph_objs as go
        import numpy as np
        df = pd.read_csv(input.table, sep="\t", index_col=0)
        tag2name = pd.read_csv(input.tag2name, index_col=0)
        tag2name.index = tag2name.index.str.replace("'","")
        df = pd.concat((df, tag2name), axis=1)
        df = df.dropna()

        df["-log10padj"] = -1 * np.log10(df["padj"])
        sig = df[df["gene_name"].isin(GOI)]
        nsig = df[~df["gene_name"].isin(GOI)]
        utrs = df[df.index.str.contains("UTR")]
        fig = go.Figure()
        sig["gene_name"] = sig["gene_name"].str.cat(sig.index, sep="-")
        nsig["gene_name"] = nsig["gene_name"].str.cat(nsig.index, sep="-")
        utrs["gene_name"] = utrs["gene_name"].str.cat(utrs.index, sep="-")
        gpval = nsig[(nsig["padj"]<= 0.05) & ((nsig["log2FoldChange"] >= 1)|(nsig["log2FoldChange"] <= -1))]
        gpval2 = sig[(sig["padj"]<= 0.05) & ((sig["log2FoldChange"] >= 1)| (sig["log2FoldChange"] <= -1))]

        nrsig = len(gpval.index) + len(gpval2.index)

        fig.add_trace(go.Scatter(
            x=nsig["log2FoldChange"],y=nsig["-log10padj"], mode="markers", marker=dict(color="blue"),
            hovertemplate=
            '<i>Y</i>: %{y:.2f}' +
            '<br><b>X</b>: %{x}<br>' +
            '<b>%{text}</b>',
            text=nsig["gene_name"],
            name="Not so Interesting Genes"

        ))
        fig.add_trace(go.Scatter(
            x=sig["log2FoldChange"], y=sig["-log10padj"], mode="markers", marker=dict(color="red", size=10),
            hovertemplate=
            '<i>Y</i>: %{y:.2f}' +
            '<br><b>X</b>: %{x}<br>' +
            '<b>%{text}</b>',
            text=sig["gene_name"],
            name="Interesting Genes"
        ))
        fig.add_trace(go.Scatter(
            x=utrs["log2FoldChange"],y=utrs["-log10padj"],mode="markers",marker=dict(color="green"),
            hovertemplate=
            '<i>Y</i>: %{y:.2f}' +
            '<br><b>X</b>: %{x}<br>' +
            '<b>%{text}</b>',
            text=utrs["gene_name"],
            name="UTRs",
            visible="legendonly"
        ))
        ptitle = f"Set:{wildcards.mode} Condition:{CONDIMAP[wildcards.condition]} vs Base:{CONDIMAP[wildcards.baseline]} significant:{nrsig}"
        fig.update_layout(
            title=ptitle,
            xaxis_title="Log2FoldChange",
            yaxis_title="-Log10(pval)",
            legend_title="Interesting or not",
            font=dict(
                size=18,
            )
        )
        fig.add_hline(-1* np.log10(0.05))
        fig.write_html(output.plot)

rule PlotPCA:
    input:
        pca_data = "PipelineData/IntermediateData/set{mode}_design_{design}_normed_pca.csv",
    output:
        pca_plot = "PipelineData/Plots/PCA/set{mode}_design_{design}_pca_plot.html"
    run:
        import pandas as pd
        import plotly.express as px
        df = pd.read_csv(input.pca_data, index_col=0)
        fig = px.scatter(df, x="PC1", y="PC2", color="sample_group")
        fig.write_html(output.pca_plot)



rule KEGGEnrichment:
    input:
        defile = "PipelineData/Tables/DESeq/set{mode}_design_{design}_cond{condition}_vs_{baseline}.tsv",
        setup = rules.setup.output.install_file
    output:
        up="PipelineData/Tables/Enrichment/KEGG_set{mode}_design_{design}_cond{condition}_vs_{baseline}_ud_upregulated.tsv",
        down="PipelineData/Tables/Enrichment/KEGG_set{mode}_design_{design}_cond{condition}_vs_{baseline}_ud_downregulated.tsv"
    script: "keggenrichment.R"

rule GOEnrichment:
    input:
        setup = rules.GOSetup.output.finished_file,
        defile = "PipelineData/Tables/DESeq/set{mode}_design_{design}_cond{condition}_vs_base{baseline}.tsv",
    output:
        up = "PipelineData/Tables/Enrichment/GO_set{mode}_design_{design}_cond{condition}_vs_base{baseline}_ud_upregulated.tsv",
        down = "PipelineData/Tables/Enrichment/GO_set{mode}_design_{design}_cond{condition}_vs_base{baseline}_ud_downregulated.tsv"
    script: "GOEnrichment.R"


rule PlotEnrichment:
    input:
        table = "PipelineData/Tables/Enrichment/{enrichment}_set{mode}_design_{design}_cond{condition}_vs_base{baseline}_ud_{updown}.tsv"
    output:
        plot =  "PipelineData/Plots/Enrichment/{enrichment}_set{mode}_design_{design}_cond{condition}_vs_base{baseline}_ud_{updown}.html"
    run:
        import plotly.graph_objects as go
        import plotly.express as px
        import pandas as pd
        enriched = pd.read_csv(input.table, sep="\t")
        if len(enriched) == 0:
            raise IndexError("Table is empty: Nothing to plot")
        fig = px.bar(enriched, x="Count", y="Description", color="p.adjust")
        ptitle = f"{wildcards.updown.upper()} Set:{wildcards.mode} Condition:{CONDIMAP[wildcards.condition]} vs Base:{CONDIMAP[wildcards.baseline]}"

        fig.update_layout(
            title=ptitle,
            legend_title="Interesting or not",
            font=dict(
                size=18,
            ),
            yaxis = {'categoryorder': 'total ascending'}
        )
        #fig.update_coloraxes(reversescale=True)
        fig.write_html(output.plot)
