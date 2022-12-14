
GOI = ["psbA2", "rnpB", "pilA1", "psaA", "cmpC", "pixG", "pilB", "pixG;  pisG;  taxP1;  rer1"]

DESIGNS = {
    "sg": ("~ sample_group", "sample_group"),
    "pf": ("~ puromycin + fraction", "fraction"),
    "fp": ("~ fraction + puromycin", "puromycin"),
    "ldpf": ("~ light_dark + puromycin + fraction", "fraction"),
    "pfld": ("~ puromycin + fraction + light_dark", "light_dark"),
    "ldfp": ("~ light_dark + fraction + puromycin", "puromycin"),
    "ld": ("~ light_dark ", "light_dark"),
}

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
        # d = expand(
        #     "PipelineData/Plots/Light_design_{design}_pvalvsfoldchange_cond{condition}_base_{baseline}.html",
        #     condition=["M", "C"], baseline=["TC"], design=["pf"]
        # ),
        # d2 = expand(
        #     "PipelineData/Plots/Light_design_{design}_pvalvsfoldchange_cond{condition}_base_{baseline}.html",
        #     condition=["FM", "FC"], baseline=["FTC"], design=["sg"]
        # ),
        d3= expand(
            "PipelineData/Plots/Light_design_{design}_pvalvsfoldchange_cond{condition}_base_{baseline}.html",
            condition=["M"],baseline=["TC"],design=["pf"]
        ),
        # fld = expand(
        #     "PipelineData/Plots/LightDark_design_{design}_pvalvsfoldchange_cond{condition}_base_{baseline}.html",
        #     condition=["M"],baseline=["TC"],design=["ldpf"]
        # ),
        ld = expand(
            "PipelineData/Plots/LightDark_design_{design}_pvalvsfoldchange_cond{condition}_base_{baseline}.html", zip,
            condition=["L", "M"],baseline=["D", "TC"],design=["pfld", "ldpf"]
        ),
        ldn= expand(
        "PipelineData/Plots/NormLightDark_design_{design}_pvalvsfoldchange_cond{condition}_base_{baseline}.html",zip,
            condition=["L"],baseline=["D"],design=["ld"]
            ),
        pca = expand(
            "PipelineData/Plots/{mode}_design_{design}_pca_plot.html",
            mode=["LightDark"], design=["ldpf"]
        )




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
        df.to_csv(output.all_counts, sep="\t")
        annotation = pd.DataFrame(data)
        annotation.to_csv(output.annotation, index=False, sep="\t")

rule DESeq2:
    input:
        counts = rules.joinDataFrames.output.all_counts,
        annotation = rules.joinDataFrames.output.annotation
    output:
        pca_data = "PipelineData/Analyzed/Light_design_{design}_normed_pca.csv",
        heatmap = "PipelineData/Plots/Light_design_{design}_correltation.pdf",
        dispest = "PipelineData/Plots/Light_design_{design}_dispersionestimates.png",
        deseq_result = "PipelineData/Analyzed/Light_design_{design}_deseq_res_obj.RData"
    params:
        design = lambda wildcards: DESIGNS[wildcards.design][0]
    script:
        "deseq.R"




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
        for col in df.columns:
            name = "_".join(col.split("_")[1:5])
            name = name.replace("plusP", "plus_P")
            name = name.replace("minP", "minus_P")
            print(name)
            data["name"].append(name)
            mp, p, frac, ld, rep = name.split("_")
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
        df.to_csv(output.table, sep="\t")
        annotation = pd.DataFrame(data)
        annotation.to_csv(output.annotation, index=False, sep="\t")



rule dropUnusedColumns:
    input:
        ld = rules.clean_light_dark.output.table,
        ldann = rules.clean_light_dark.output.annotation,
    output:
        counts = "PipelineData/Data/cleaned_LightDarkLight.csv",
        annotation =  "PipelineData/Data/annotation_LightDarkLight.csv"
    run:
        import pandas as pd
        ld = pd.read_csv(input.ld, sep="\t", index_col=0)
        ldann = pd.read_csv(input.ldann, sep="\t", index_col=0)
        df = ld.loc[:, (ld.columns.str.contains("minus_P_TC_L_")) | (ld.columns.str.contains("minus_P_TC_D_"))]
        ann = ldann[(ldann.index.str.contains("minus_P_TC_L_")) | (ldann.index.str.contains("minus_P_TC_D_"))]
        df.to_csv(output.counts, sep="\t")
        ann.to_csv(output.annotation, sep="\t")




rule extract_result_Light:
    input:
        deseq_result = "PipelineData/Analyzed/Light_design_{design}_deseq_res_obj.RData"
    params:
        design = lambda wildcards: DESIGNS[wildcards.design][1]
    output:
        result_table = "PipelineData/Tables/Light_design_{design}_deseq_{condition}_vs_{baseline}.tsv"
    script: "extractDESeqResult.R"



rule LightDarkNormDESeq:
    input:
        counts = rules.dropUnusedColumns.output.counts,
        annotation = rules.dropUnusedColumns.output.annotation,
        norm_genes = "Data/normGene.tsv"
    output:
        pca_data = "PipelineData/Analyzed/NormLightDark_design_{design}_normed_pca.csv",
        heatmap = "PipelineData/Plots/NormLightDark_design_{design}_correltation.pdf",
        dispest = "PipelineData/Plots/NormLightDark_design_{design}_dispersionestimates.png",
        deseq_result = "PipelineData/Analyzed/NormLightDark_design_{design}_deseq_res_obj.RData"
    params:
        design = lambda wildcards: DESIGNS[wildcards.design][0]
    script:
        "deseqLightDark.R"

rule extract_result_NormLightDark:
    input:
        deseq_result = "PipelineData/Analyzed/NormLightDark_design_{design}_deseq_res_obj.RData"
    params:
        design = lambda wildcards: DESIGNS[wildcards.design][1]
    output:
        result_table = "PipelineData/Tables/NormLightDark_design_{design}_deseq_{condition}_vs_{baseline}.tsv"
    script: "extractDESeqResult.R"


rule least_significant_LightDark:
    input:
        ld = expand(
            "PipelineData/Tables/NormLightDark_design_{design}_deseq_{condition}_vs_{baseline}.tsv", zip,
            design=["ld"], condition=["D"], baseline=["L"]
        ),
        tcm = "PipelineData/Tables/Light_design_pf_deseq_M_vs_TC.tsv"
    output:
        least_significant_list = "PipelineData/Tables/least_sig_NormLightDark.tsv"
    run:
        import pandas as pd
        dfs = [input.ld[0], input.tcm]
        subdfs = []
        for file in dfs:
            df = pd.read_csv(file, sep="\t", index_col=0)
            df["absChange"] = df["log2FoldChange"].abs()
            subdf = df
            subdf = df[df["absChange"] < 1]
            subdf = subdf[subdf["baseMean"] > 10]
            subdf = subdf[subdf["padj"] >= 0.05]
            subdf = subdf.sort_values(by=["padj", "absChange"], ascending=[False, True])
            #subdf = subdf.iloc[0:200]
            subdf.index = subdf.index.astype("string")
            subdfs.append(subdf)
        intersect = set(subdfs[0].index)
        for df in subdfs[1:]:
            intersect = intersect & set(df.index)

        outdf = pd.DataFrame(
            index=list(intersect)
        )
        outdf.to_csv(output.least_significant_list, sep="\t")


rule LightDarkDESeq:
    input:
        counts = rules.clean_light_dark.output.table,
        annotation = rules.clean_light_dark.output.annotation,
        norm_genes = rules.least_significant_LightDark.output.least_significant_list
    output:
        pca_data = "PipelineData/Analyzed/LightDark_design_{design}_normed_pca.csv",
        heatmap = "PipelineData/Plots/LightDark_design_{design}_correltation.pdf",
        dispest = "PipelineData/Plots/LightDark_design_{design}_dispersionestimates.png",
        deseq_result = "PipelineData/Analyzed/LightDark_design_{design}_deseq_res_obj.RData"
    params:
        design = lambda wildcards: DESIGNS[wildcards.design][0]
    script:
        "deseqLightDark.R"

rule extract_result_LightDark:
    input:
        deseq_result = "PipelineData/Analyzed/LightDark_design_{design}_deseq_res_obj.RData"
    params:
        design = lambda wildcards: DESIGNS[wildcards.design][1]
    output:
        result_table = "PipelineData/Tables/LightDark_design_{design}_deseq_{condition}_vs_{baseline}.tsv"
    script: "extractDESeqResult.R"


rule plotFoldPval:
    input:
        table = "PipelineData/Tables/{mode}_design_{design}_deseq_{condition}_vs_{baseline}.tsv",
        tag2name = "locusTagtoGeneName.csv"
    output:
        plot = "PipelineData/Plots/{mode}_design_{design}_pvalvsfoldchange_cond{condition}_base_{baseline}.html"
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
        nsig = df[df["padj"] > 0]
        fig = go.Figure()

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
        ptitle = f"{wildcards.mode} {CONDIMAP[wildcards.condition]} vs {CONDIMAP[wildcards.baseline]}"
        fig.update_layout(
            title=ptitle,
            xaxis_title="Log2FoldChange",
            yaxis_title="-Log10(pval)",
            legend_title="Interesting or not",
            font=dict(
                size=18,
            )
        )
        fig.add_hline(-1* np.log10(0.005))
        fig.write_html(output.plot)

rule PlotPCA:
    input:
        pca_data = "PipelineData/Analyzed/{mode}_design_{design}_normed_pca.csv",
    output:
        pca_plot = "PipelineData/Plots/{mode}_design_{design}_pca_plot.html"
    run:
        import pandas as pd
        import plotly.express as px
        df = pd.read_csv(input.pca_data, index_col=0)
        fig = px.scatter(df, x="PC1", y="PC2", color="sample_group")
        fig.write_html(output.pca_plot)



