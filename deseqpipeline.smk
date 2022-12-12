
GOI = ["psbA2", "rnpB"]

rule all:
    input:
        d = expand(
            "PipelineData/Plots/pvalvsfoldchange_cond{condition}_base_{baseline}.html",
            condition=["M", "C"], baseline=["TC"]
        )


rule joinDataFrames:
    input:
        counts = "PipelineData/Data/All_CDS.csv",
        spike_ins = "PipelineData/Data/All_ERCC.csv"
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
        #dedata = "Analyzed/deseq_data.csv",
        pca_data = "PipelineData/Analyzed/normed_pca.csv",
        heatmap = "PipelineData/Plots/correltation.pdf",
        dispest = "PipelineData/Plots/dispersionestimates.png",
        deseq_result = "PipelineData/Analyzed/deseq_res_obj.RData"
    script:
        "deseq.R"

rule PlotPCA:
    input:
        pca_data = rules.DESeq2.output.pca_data
    output:
        pca_plot = "PipelineData/Plots/pca_plot.html"
    run:
        import pandas as pd
        import plotly.express as px
        df = pd.read_csv(input.pca_data, index_col=0)
        fig = px.scatter(df, x="PC1", y="PC2", color="sample_group")


rule extract_result:
    input:
        deseq_result = rules.DESeq2.output.deseq_result
    output:
        result_table = "PipelineData/Tables/deseq_{condition}_vs_{baseline}.tsv"
    script: "extractDESeqResult.R"


rule plotFoldPval:
    input:
        table = rules.extract_result.output.result_table,
        tag2name = "locusTagtoGeneName.csv"
    output:
        plot = "PipelineData/Plots/pvalvsfoldchange_cond{condition}_base_{baseline}.html"
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
        sig = df[df["padj"] <= 0.05]
        sig = df[df["gene_name"].isin(GOI)]
        nsig = df[df["padj"] > 0]
        fig = go.Figure()

        fig.add_trace(go.Scatter(
            x=nsig["log2FoldChange"],y=nsig["-log10padj"], mode="markers", marker=dict(color="blue"),
            hovertemplate=
            '<i>Y</i>: %{y:.2f}' +
            '<br><b>X</b>: %{x}<br>' +
            '<b>%{text}</b>',
            text=nsig["gene_name"]
        ))
        fig.add_trace(go.Scatter(
            x=sig["log2FoldChange"], y=sig["-log10padj"], mode="markers", marker=dict(color="red"),
            hovertemplate=
            '<i>Y</i>: %{y:.2f}' +
            '<br><b>X</b>: %{x}<br>' +
            '<b>%{text}</b>',
            text=sig["gene_name"]
        ))
        fig.add_hline(-1* np.log10(0.005))
        fig.write_html(output.plot)




