import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from matplotlib import cm


def parse_detailed_log_multi_ds(file_path):
    with open(file_path, 'r') as f:
        content = f.read()

    runs = content.split("================== RUN ==================")
    # Structure: data[dataset][r][s] = {'cores': [], 'comps': [], 'nodes': []}
    data = {}

    for run in runs:
        if "dataset" not in run: continue
        ds_match = re.search(r"dataset\s+:\s+.*/(.*?)\.edges", run)
        if not ds_match: continue
        ds_name = ds_match.group(1)

        r_match = re.search(r"r\s+:\s+(\d+)", run)
        s_match = re.search(r"s\s+:\s+(\d+)", run)
        if not (r_match and s_match): continue
        r_val, s_val = int(r_match.group(1)), int(s_match.group(1))

        stats_match = re.search(r"Core Component Stats=+\n(.*?)\n=+", run, re.DOTALL)
        if stats_match:
            lines = re.findall(r"Core:\s+([\d.e+-]+)\s+Components:\s+([\d.e+-]+)\s+LivingNodes:\s+([\d.e+-]+)",
                               stats_match.group(1))
            if lines:
                cores = [float(x[0]) for x in lines]
                comps = [float(x[1]) for x in lines]
                nodes = [float(x[2]) for x in lines]

                user_time_sec = None
                max_rss_kb = None
                t_match = re.search(r"User time \(seconds\):\s+([\d.]+)", run)
                if t_match:
                    user_time_sec = float(t_match.group(1))
                m_match = re.search(r"Maximum resident set size \(kbytes\):\s+(\d+)", run)
                if m_match:
                    max_rss_kb = int(m_match.group(1))

                if ds_name not in data:
                    data[ds_name] = {}
                if r_val not in data[ds_name]:
                    data[ds_name][r_val] = {}

                data[ds_name][r_val][s_val] = {
                    'cores': cores, 'comps': comps, 'nodes': nodes,
                    'user_time_sec': user_time_sec, 'max_rss_kb': max_rss_kb
                }
    return data


def flatten_runs_to_df(data):
    rows = []
    for ds_name, r_dict in data.items():
        for r_val, s_dict in r_dict.items():
            for s_val, metrics in s_dict.items():
                cores = metrics.get('cores') or []
                comps = metrics.get('comps') or []
                nodes = metrics.get('nodes') or []
                user_time_sec = metrics.get('user_time_sec')
                max_rss_kb = metrics.get('max_rss_kb')

                # Prefer per-core rows; fall back to a single aggregate row if cores are missing
                if cores and comps and nodes:
                    for core, comp, node in zip(cores, comps, nodes):
                        rows.append({
                            'dataset': ds_name,
                            'r': r_val,
                            's': s_val,
                            'core': core,
                            'components': comp,
                            'living_nodes': node,
                            'user_time_sec': user_time_sec,
                            'max_rss_kb': max_rss_kb,
                        })
                else:
                    rows.append({
                        'dataset': ds_name,
                        'r': r_val,
                        's': s_val,
                        'core': None,
                        'components': None,
                        'living_nodes': None,
                        'user_time_sec': user_time_sec,
                        'max_rss_kb': max_rss_kb,
                    })
    return pd.DataFrame(rows)


def plot_fixed_s_varying_r(data):
    # Set academic plotting style
    plt.rcParams.update({
        'font.family': 'serif', 'font.size': 14,  # Slightly larger font for single big plots
        'axes.labelsize': 14,
        'axes.titlesize': 16,
        'legend.fontsize': 12,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'savefig.dpi': 300,
        'mathtext.fontset': 'stix', 'pdf.fonttype': 42, 'ps.fonttype': 42,
    })

    base_out = "case_study_fixed_s"
    os.makedirs(base_out, exist_ok=True)

    # 1. Pivot data to be accessed by [dataset][s][r]
    pivoted_data = {}
    for ds_name, r_dict in data.items():
        if ds_name not in pivoted_data: pivoted_data[ds_name] = {}
        for r_val, s_dict in r_dict.items():
            for s_val, metrics in s_dict.items():
                if s_val not in pivoted_data[ds_name]: pivoted_data[ds_name][s_val] = {}
                pivoted_data[ds_name][s_val][r_val] = metrics

    # 2. Plotting loop
    for ds_name, s_dict in pivoted_data.items():
        ds_out = os.path.join(base_out, ds_name)
        os.makedirs(ds_out, exist_ok=True)

        # Sort s values to generate ordered plots
        sorted_s_vals = sorted(s_dict.keys())

        for s_val in sorted_s_vals:
            r_data_map = s_dict[s_val]
            # Need at least one line to plot
            if not r_data_map: continue

            sorted_rs = sorted(r_data_map.keys())

            # Use distinct colors for different r values (usually few, e.g., 1,2,3,4)
            # Plasma or Autumn allows clear distinction between low and high r
            colors = cm.plasma(np.linspace(0, 0.8, len(sorted_rs)))

            # --- Plot 1: Component Evolution (Fixed s, Varying r) ---
            fig_comp, ax_comp = plt.subplots(figsize=(8, 6))
            min_core_val = float('inf')

            for idx, r_val in enumerate(sorted_rs):
                d = r_data_map[r_val]
                # Plot Line
                ax_comp.plot(d['cores'], d['comps'], color=colors[idx],
                             marker='o', markersize=4, linewidth=2, alpha=0.8,
                             label=f'r={r_val}')  # Add label for legend as backup

                # Direct Labeling at START
                start_x = d['cores'][0]
                start_y = d['comps'][0]
                min_core_val = min(min_core_val, start_x)

                ax_comp.text(start_x * 0.9, start_y, f'r={r_val} ',
                             color=colors[idx], ha='right', va='center',
                             fontsize=12, fontweight='bold')

            ax_comp.set_xscale('log')
            ax_comp.set_yscale('log')
            ax_comp.set_xlabel('Core Value ($k$)')
            ax_comp.set_ylabel('Component Count')
            ax_comp.set_title(f'{ds_name}: Components (s={s_val})')
            ax_comp.set_xlim(left=min_core_val * 0.5)
            ax_comp.grid(True, which="both", ls="-", alpha=0.1)

            comp_path = os.path.join(ds_out, f'{ds_name}_s{s_val}_components.png')
            fig_comp.tight_layout()
            fig_comp.savefig(comp_path)
            plt.close(fig_comp)

            # --- Plot 2: Living Node Evolution (Fixed s, Varying r) ---
            fig_node, ax_node = plt.subplots(figsize=(8, 6))
            min_core_val = float('inf')

            for idx, r_val in enumerate(sorted_rs):
                d = r_data_map[r_val]
                ax_node.plot(d['cores'], d['nodes'], color=colors[idx],
                             marker='^', markersize=4, linewidth=2, alpha=0.8,
                             label=f'r={r_val}')

                # Direct Labeling at START
                start_x = d['cores'][0]
                start_y = d['nodes'][0]
                min_core_val = min(min_core_val, start_x)

                ax_node.text(start_x * 0.9, start_y, f'r={r_val} ',
                             color=colors[idx], ha='right', va='center',
                             fontsize=12, fontweight='bold')

            ax_node.set_xscale('log')
            ax_node.set_yscale('log')
            ax_node.set_xlabel('Core Value ($k$)')
            ax_node.set_ylabel('Living Node Count')
            ax_node.set_title(f'{ds_name}: Living Nodes (s={s_val})')
            ax_node.set_xlim(left=min_core_val * 0.5)
            ax_node.grid(True, which="both", ls="-", alpha=0.1)

            node_path = os.path.join(ds_out, f'{ds_name}_s{s_val}_nodes.png')
            fig_node.tight_layout()
            fig_node.savefig(node_path)
            plt.close(fig_node)

            print(f"Generated plots for {ds_name} s={s_val} with r={sorted_rs}")

def plot_fixed_r_varying_s(data):
    # Set academic plotting style
    plt.rcParams.update({
        'font.family': 'serif', 'font.size': 11, 'axes.labelsize': 12,
        'axes.titlesize': 13, 'savefig.dpi': 300,
        'mathtext.fontset': 'stix', 'pdf.fonttype': 42, 'ps.fonttype': 42,
    })

    base_out = "case_study_fixed_r"
    os.makedirs(base_out, exist_ok=True)

    for ds_name, r_dict in data.items():
        ds_out = os.path.join(base_out, ds_name)
        os.makedirs(ds_out, exist_ok=True)

        for r_val, s_dict in r_dict.items():
            sorted_ss = sorted(s_dict.keys())
            # Use a distinctive colormap
            colors = cm.viridis(np.linspace(0, 0.9, len(sorted_ss)))

            # --- Plot 1: Component Evolution ---
            fig_comp, ax_comp = plt.subplots(figsize=(8, 5))
            min_core_val = float('inf')

            for idx, s_val in enumerate(sorted_ss):
                d = s_dict[s_val]
                ax_comp.plot(d['cores'], d['comps'], color=colors[idx],
                             marker='o', markersize=2, linewidth=1.5, alpha=0.8)

                # Direct Labeling at START
                start_x = d['cores'][0]
                start_y = d['comps'][0]
                min_core_val = min(min_core_val, start_x)

                # Label text with s=...
                # Offset slightly to the left
                ax_comp.text(start_x * 0.95, start_y, f's={s_val} ',
                             color=colors[idx], ha='right', va='center',
                             fontsize=9, fontweight='bold')

            ax_comp.set_xscale('log')
            ax_comp.set_yscale('log')
            ax_comp.set_xlabel('Core Value ($k$)')
            ax_comp.set_ylabel('Component Count')
            ax_comp.set_title(f'{ds_name} Components (r={r_val})')

            # Adjust x-limit to make room for labels on the left
            ax_comp.set_xlim(left=min_core_val * 0.5)
            ax_comp.grid(True, which="both", ls="-", alpha=0.1)

            comp_path = os.path.join(ds_out, f'{ds_name}_r{r_val}_components_startlabel.png')
            fig_comp.savefig(comp_path, bbox_inches='tight')
            fig_comp.savefig(comp_path.replace('.png', '.eps'), bbox_inches='tight')
            plt.close(fig_comp)

            # --- Plot 2: Living Node Evolution ---
            fig_node, ax_node = plt.subplots(figsize=(8, 5))
            min_core_val = float('inf')

            for idx, s_val in enumerate(sorted_ss):
                d = s_dict[s_val]
                ax_node.plot(d['cores'], d['nodes'], color=colors[idx],
                             marker='^', markersize=2, linewidth=1.5, alpha=0.8)

                # Direct Labeling at START
                start_x = d['cores'][0]
                start_y = d['nodes'][0]
                min_core_val = min(min_core_val, start_x)

                ax_node.text(start_x * 0.95, start_y, f's={s_val} ',
                             color=colors[idx], ha='right', va='center',
                             fontsize=9, fontweight='bold')

            ax_node.set_xscale('log')
            ax_node.set_yscale('log')
            ax_node.set_xlabel('Core Value ($k$)')
            ax_node.set_ylabel('Living Node Count')
            ax_node.set_title(f'{ds_name} Living Nodes (r={r_val})')

            # Adjust x-limit
            ax_node.set_xlim(left=min_core_val * 0.5)
            ax_node.grid(True, which="both", ls="-", alpha=0.1)

            node_path = os.path.join(ds_out, f'{ds_name}_r{r_val}_nodes_startlabel.png')
            fig_node.savefig(node_path, bbox_inches='tight')
            fig_node.savefig(node_path.replace('.png', '.eps'), bbox_inches='tight')
            plt.close(fig_node)

            print(f"Generated start-labeled plots for {ds_name} r={r_val}")
            print(f"Paths saved in {ds_out}")


# Main execution
log_file = "data/dblp_casestudy_server.log"
data_ready = parse_detailed_log_multi_ds(log_file)
# data_ready.to_csv('data/parsed_dblp_casestudy_server.csv', index=False)
# data_ready.to_csv('/Users/zhangwenqian/Library/CloudStorage/Dropbox/应用/Overleaf/Nuclear CD/figure/pokec_casestudy_data.csv', index=False)
df = flatten_runs_to_df(data_ready)
df.to_csv('/Users/zhangwenqian/Library/CloudStorage/Dropbox/应用/Overleaf/Nuclear CD/figure/pokec_casestudy_data.csv', index=False)
plot_fixed_r_varying_s(data_ready)
plot_fixed_s_varying_r(data_ready)