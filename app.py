import matplotlib
matplotlib.use("Agg")
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import os
from pathlib import Path
from upsetplot import UpSet, from_contents
import matplotlib.pyplot as plt
from collections import defaultdict
import networkx as nx
from pyvis.network import Network
import tempfile

# ===============================
# Streamlit config
# ===============================
st.set_page_config(
    page_title="Neurodegeneration Regulon Database",
    layout="wide"
)

st.sidebar.image("images/regulome.png", width=250)
st.sidebar.title("Neurodegeneration Regulon Database")

st.sidebar.markdown(
    """
    ** Neurodegeneration Regulon Database** is a reconstructed transcriptional regulatory network resource actually focused on **Alzheimer’s disease (AD)**.

    It integrates:
    - **Inferred transcription factor (TF) regulons**
    - **Brain-region–specific gene expression**
    - **Differential expression and enrichment analyses**
    - **Regulatory network reconstruction**
    This database is a resource to explore **deregulated transcription factors and regulons**
    in AD pathology.
    """
)

st.sidebar.markdown("---")

st.sidebar.subheader("Data included")

st.sidebar.markdown(
    """
    - TF–target interactions  
    - Region-specific regulons  
    - GO Biological Process, CC, and MF  
    - Cross-region overlap analysis  
    """
)

st.sidebar.markdown("---")

st.sidebar.subheader("Citation")
st.sidebar.markdown(
    """
    *Belém-Souza MV, Barra-Matos G,  
    Santana de Araújo G.*  
    **Regulon reconstruction uncovers novel deregulated transcription
    factors in Alzheimer’s disease**  
    *Manuscript under review*
    """
)

st.sidebar.markdown("---")

st.sidebar.caption(
    "Laboratory of Bioinformatics and Data Science (LBCD) · UFPA"
)

tab1, tab2, tab3, tab4, tab5 = st.tabs(["DGE analysis", "Gene Ontology", "GO Overlap", "Master-regulons", "Supporting TF-Link, AD-PPI data"])

with tab1:

    st.title("Differential Gene Expression Analysis (DEG) - Explorer")
    st.caption("Multi-cohort, multi-region - DEG visualization")

    # ===============================
    # Load files
    # ===============================
    DGE_DATA_DIR = "degs/"

    all_files = [f for f in os.listdir(DGE_DATA_DIR) if f.endswith(".txt")]
    file_names = [os.path.splitext(f)[0] for f in all_files]

    col1, col2, col3 = st.columns(3)
    with col1:
        selected_name = st.selectbox(
            "Select dataset",
            file_names
        )
        selected_file = os.path.join(DGE_DATA_DIR, selected_name + ".txt")

        st.write(f"You selected: {selected_name}")

        fc_thr = st.slider(
            "|logFC| threshold",
            min_value=0.0,
            max_value=10.0,
            value=0.5,
            step=0.1
        )

        fdr_thr = st.slider(
            "FDR threshold",
            min_value=0.0,
            max_value=0.2,
            value=0.05,
            step=0.005
        )
    with col2:
        # ===============================
        # Load DEG table
        # ===============================
        @st.cache_data
        def load_deg(path):
            df = pd.read_csv(path, sep="\t")
            return df

        df = load_deg(os.path.join(selected_file))
        del df['keyvals']
        # ===============================
        # DEG classification
        # ===============================
        df["significant"] = (
            (df["FDR"] <= fdr_thr) &
            (df["logFC"].abs() >= fc_thr)
        )

        df["Direction"] = np.where(
            (df["significant"]) & (df["logFC"] > 0), "Up",
            np.where((df["significant"]) & (df["logFC"] < 0), "Down", "NS")
        )

        # ===============================
        # Summary
        # ===============================
        col1, col2, col3 = st.columns(3)

        col1.metric("Total genes", df.shape[0])
        col2.metric("Up-regulated", (df["Direction"] == "Up").sum())
        col3.metric("Down-regulated", (df["Direction"] == "Down").sum())

        # ===============================
        # Volcano plot
        # ===============================
        st.subheader("Volcano plot")

        volcano = px.scatter(
            df,
            x="logFC",
            y=-np.log10(df["FDR"]),
            color="Direction",
            color_discrete_map={
                "Up": "#d62728",
                "Down": "#1f77b4",
                "NS": "lightgray"
            },
            hover_data=["gene"] if "gene" in df.columns else None,
            opacity=0.7
        )

        volcano.add_vline(x=fc_thr, line_dash="dash", line_color="black")
        volcano.add_vline(x=-fc_thr, line_dash="dash", line_color="black")
        volcano.add_hline(y=-np.log10(fdr_thr), line_dash="dash", line_color="black")

        volcano.update_layout(
            height=550,
            legend_title="Direction",
            xaxis_title="logFC",
            yaxis_title="-log10(FDR)"
        )

        st.plotly_chart(volcano, use_container_width=True)

    with col3:
        pass

    # ===============================
    # DEG table
    # ===============================
    st.subheader("Differentially expressed genes")

    deg_table = df[df["Direction"] != "NS"] \
        .sort_values("FDR") \
        .reset_index(drop=True)

    st.dataframe(
        deg_table,
        use_container_width=True,
        height=400
    )

    # ===============================
    # Download
    # ===============================
    st.download_button(
        label="Download filtered DGE analysis table",
        data=deg_table.to_csv(sep="\t", index=False),
        file_name=f"{selected_file.replace('.txt','')}_DEGs.tsv",
        mime="text/tab-separated-values"
    )

with tab2:
    GO_DATA_DIR = "go/"

    st.set_page_config(page_title="GO Explorer", layout="wide")
    st.title("GO Enrichment Explorer")
    st.caption("Filter by tissue and explore BP, CC, MF GO tables")

    all_files = [f for f in os.listdir(GO_DATA_DIR) if f.endswith(".tsv")]

    tissues = sorted(list(set([f.split("_")[-1].replace(".tsv", "") for f in all_files])))
    categories = ["BP", "CC", "MF"]

    selected_tissue = st.selectbox("Select tissue", tissues)


    @st.cache_data
    def load_go_table(category, tissue):
        fname = f"GO_{category}_{tissue}.tsv"
        path = os.path.join(GO_DATA_DIR, fname)
        if os.path.exists(path):
            df = pd.read_csv(path, sep="\t")
            return df
        else:
            return pd.DataFrame()  # empty if file not found

    for cat in categories:
        st.subheader(f"GO {cat} - {selected_tissue}")
        df_go = load_go_table(cat, selected_tissue)
        if df_go.empty:
            st.write("No data available")
        else:
            st.dataframe(df_go, use_container_width=True)

with tab3:

    st.set_page_config(
        page_title="GO Term overlap by tissue",
        layout="wide"
    )

    st.title("GO Term overlap across brain regions")
    st.markdown(
        """
        Overlap of enriched **Gene Ontology (GO) terms** across brain regions,
        visualized using **UpSet plots**, separated by ontology (BP, CC, MF).
        """
    )

    # ===============================
    # Controls (top of page)
    # ===============================
    col1, col2, col3 = st.columns([2, 1, 1])

    go_dir = "go"

    with col1:
        ontology = st.selectbox(
            "Ontology",
            ["BP", "CC", "MF"]
        )

    with col2:
        min_subset_size = st.number_input(
            "Min GO terms per intersection",
            min_value=1,
            max_value=100,
            value=5
        )

    show_table = st.checkbox(
        "Show table of shared GO terms",
        value=True
    )

    st.divider()

    # ===============================
    # Load GO files on cache
    # ===============================
    @st.cache_data
    def load_go_files(go_dir):
        go_terms = defaultdict(dict)

        for f in go_dir.glob("GO_*.tsv"):
            parts = f.stem.split("_")
            ontology = parts[1]
            tissue = "_".join(parts[2:])

            df = pd.read_csv(f, sep="\t")

            if "ID" in df.columns:
                go_ids = set(df["ID"])
            elif "GO.ID" in df.columns:
                go_ids = set(df["GO.ID"])
            else:
                raise ValueError(f"No GO ID column found in {f.name}")

            go_terms[ontology][tissue] = go_ids

        return go_terms


    go_dir = Path(go_dir)

    if not go_dir.exists():
        st.error("GO directory not found.")
        st.stop()

    go_terms = load_go_files(go_dir)

    if ontology not in go_terms:
        st.error(f"No files found for ontology {ontology}")
        st.stop()

    tissues = sorted(go_terms[ontology].keys())

    st.subheader("Tissue selection")

    selected_tissues = st.multiselect(
        "Select tissues",
        tissues,
        default=tissues
    )

    if len(selected_tissues) < 2:
        st.warning("Please select at least two tissues.")
        st.stop()

    st.divider()

    contents = {
        tissue: go_terms[ontology][tissue]
        for tissue in selected_tissues
    }

    upset_data = from_contents(contents)

    st.subheader(f"{ontology} GO term overlap")

    # fig = plt.figure(figsize=(10, 6))
    # UpSet(
    #     upset_data,
    #     min_subset_size=min_subset_size,
    #     show_counts=True,
    #     sort_by="cardinality"
    # ).plot(fig=fig)
    # st.pyplot(fig)
    # plt.close(fig)

    if show_table:
        st.subheader("Shared GO terms")

        go_to_tissues = defaultdict(set)
        go_desc = {}

        for f in go_dir.glob(f"GO_{ontology}_*.tsv"):
            tissue = "_".join(f.stem.split("_")[2:])

            if tissue not in selected_tissues:
                continue

            df = pd.read_csv(f, sep="\t")

            # handle both column conventions
            if "GO.ID" in df.columns:
                go_col = "GO.ID"
            elif "ID" in df.columns:
                go_col = "ID"
            else:
                continue

            desc_col = "Description" if "Description" in df.columns else None

            for _, row in df.iterrows():
                go = row[go_col]
                go_to_tissues[go].add(tissue)

                if desc_col:
                    go_desc[go] = row[desc_col]

        shared_df = (
            pd.DataFrame({
                "GO_ID": list(go_to_tissues.keys()),
                "GO_description": [go_desc.get(go, "NA") for go in go_to_tissues.keys()],
                "Tissues": [", ".join(sorted(t)) for t in go_to_tissues.values()],
                "Tissue_count": [len(t) for t in go_to_tissues.values()]
            })
            .sort_values("Tissue_count", ascending=False)
            .reset_index(drop=True)
        )

        if shared_df.empty:
            st.warning("No shared GO terms found for the selected tissues.")
        else:
            st.dataframe(shared_df, use_container_width=True, height=450)

        st.caption(
            "Shared GO terms across selected tissues with GO descriptions."
        )

with tab4:
    # Load regulon data
    df_regulons = pd.read_table("analysis/regulons_TFs_masters.txt")
    df_regulons["mode"] = pd.to_numeric(df_regulons["mode"], errors="coerce")

    st.subheader("TF-centered regulon networks grouped by regulatory mode")

    # Selectors
    tissue = st.selectbox(
        "Select tissue",
        sorted(df_regulons["tissue"].dropna().unique())
    )

    tf = st.selectbox(
        "Select transcription factor",
        sorted(
            df_regulons.loc[df_regulons["tissue"] == tissue, "tf"]
            .dropna()
            .unique()
        )
    )

    df_filt = df_regulons[
        (df_regulons["tissue"] == tissue) &
        (df_regulons["tf"] == tf)
        ]

    if df_filt.empty:
        st.warning("No regulon edges found for this TF in this tissue.")
        st.stop()

    #Graph
    G = nx.DiGraph()

    G.add_node(tf, color="grey", size=35, fixed=True, x=0, y=0)

    for _, row in df_filt.iterrows():
        if row["mode"] > 0:
            color = "red"
            group = "positive"
        else:
            color = "blue"
            group = "negative"

        G.add_node(
            row["target"],
            color=color,
            size=20,
            group=group
        )

        G.add_edge(
            tf,
            row["target"],
            color=color,
            title=f"mode = {row['mode']:.3f}",
            width=2 + abs(row["mode"]) * 3
        )

    net = Network(
        height="600px",
        width="100%",
        directed=True,
        bgcolor="white",
        font_color="black"
    )

    net.from_nx(G)
    net.repulsion(node_distance=220, spring_length=200)

    net.set_options("""
    var options = {
      "groups": {
        "positive": {
          "color": { "background": "red", "border": "darkred" }
        },
        "negative": {
          "color": { "background": "blue", "border": "darkblue" }
        }
      },
      "physics": {
        "barnesHut": {
          "gravitationalConstant": -30000,
          "springLength": 180,
          "springConstant": 0.04
        },
        "minVelocity": 0.75
      }
    }
    """)

    # Render
    with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmp:
        net.save_graph(tmp.name)
        html_path = tmp.name

    with open(html_path, "r", encoding="utf-8") as f:
        st.components.v1.html(f.read(), height=600)

    os.unlink(html_path)

    st.caption("🔴 Positive mode (activation-like) | 🔵 Negative mode (repression-like)")

    st.dataframe(
        df_filt,
        use_container_width=True,
        height=400
    )

with tab5:
    supporting_ppi = pd.read_csv("supporting/summary_mrtf_ppi.csv")
    supporting_tflink = pd.read_csv("supporting/summary_mrtf_tflink.csv")

    st.markdown("### Each row corresponds to an MR–TF interaction supported by summarized AD-PPI (BioGRID) evidence.")
    st.dataframe(supporting_ppi)
    st.markdown("### Each row corresponds to an MR–TF interaction supported by summarized TFLink evidence.")
    st.dataframe(supporting_tflink)