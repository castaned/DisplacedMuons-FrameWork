#!/usr/bin/env python3

import argparse
import glob
import math
import os
from array import array

import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


def collect_files(patterns):
    files = []
    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        if matches:
            files.extend(matches)
        elif os.path.isfile(pattern):
            files.append(pattern)
    return list(dict.fromkeys(files))


def require_branches(chain, names):
    missing = [name for name in names if not chain.GetBranch(name)]
    if missing:
        raise RuntimeError("Missing global/outer study branches: " + ", ".join(missing))


def print_zero_selection_diagnostics(
    path,
    tree_name,
    candidate,
    upper_candidate,
    lower_candidate,
    raw_selection,
    quality_candidate,
    resolution_candidate,
):
    diagnostic = ROOT.TChain(tree_name)
    diagnostic.Add(path)
    checks = (
        ("events in first file", ""),
        ("events with at least two DGL muons", "Sum$(dmu_isDGL)>=2"),
        ("events with at least two usable DGL outer tracks", f"Sum$({candidate})>=2"),
        (
            "events with at least one upper and one lower candidate",
            f"Sum$({upper_candidate})>=1 && Sum$({lower_candidate})>=1",
        ),
        (
            "events with exactly one upper and one lower candidate",
            f"Sum$({candidate})==2 && Sum$({upper_candidate})==1 && "
            f"Sum$({lower_candidate})==1",
        ),
        ("events passing the opposite-direction pair requirement", raw_selection),
        (
            "events where both candidates pass DGL quality cuts",
            f"({raw_selection}) && Sum$({quality_candidate})==2",
        ),
        (
            "events where both candidates also pass pT-error cuts",
            f"({raw_selection}) && Sum$({resolution_candidate})==2",
        ),
    )
    print(f"Zero-selection diagnostic using first file: {path}")
    for label, cut in checks:
        count = diagnostic.GetEntries(cut) if cut else diagnostic.GetEntries()
        print(f"  {label}: {count}")


def make_hist(chain, name, expression, selection, bins, xmin, xmax):
    hist = ROOT.TH1F(name, "", bins, xmin, xmax)
    chain.Draw(f"{expression}>>{name}", selection, "goff")
    hist.SetDirectory(0)
    return hist


def make_hist2d(
    chain, name, expression, selection, xbins, xmin, xmax, ybins, ymin, ymax
):
    hist = ROOT.TH2F(name, "", xbins, xmin, xmax, ybins, ymin, ymax)
    chain.Draw(f"{expression}>>{name}", selection, "goff")
    hist.SetDirectory(0)
    return hist


def normalize(hist):
    integral = hist.Integral(0, hist.GetNbinsX() + 1)
    if integral > 0:
        hist.Scale(1.0 / integral)


def draw_overlay(
    hist_global,
    hist_outer,
    xtitle,
    output_path,
    first_label="global-track fit",
    second_label="outer-track fit",
):
    normalize(hist_global)
    normalize(hist_outer)
    hist_global.SetLineColor(ROOT.kRed + 1)
    hist_outer.SetLineColor(ROOT.kBlue + 1)
    hist_global.SetLineWidth(2)
    hist_outer.SetLineWidth(2)
    ymax = max(hist_global.GetMaximum(), hist_outer.GetMaximum()) * 1.25
    hist_global.SetMaximum(ymax if ymax > 0 else 1.0)
    hist_global.SetMinimum(0.0)
    hist_global.GetXaxis().SetTitle(xtitle)
    hist_global.GetYaxis().SetTitle("Normalized events")

    canvas = ROOT.TCanvas(f"c_{hist_global.GetName()}", "", 900, 700)
    hist_global.Draw("HIST")
    hist_outer.Draw("HIST SAME")
    legend = ROOT.TLegend(0.62, 0.76, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(hist_global, first_label, "l")
    legend.AddEntry(hist_outer, second_label, "l")
    legend.Draw()
    canvas.SaveAs(output_path)


def draw_distribution(hist, xtitle, output_path):
    normalize(hist)
    hist.SetLineColor(ROOT.kBlue + 1)
    hist.SetLineWidth(2)
    hist.SetMinimum(0.0)
    hist.GetXaxis().SetTitle(xtitle)
    hist.GetYaxis().SetTitle("Normalized events")

    canvas = ROOT.TCanvas(f"c_{hist.GetName()}", "", 900, 700)
    hist.Draw("HIST")
    zero = ROOT.TLine(0.0, 0.0, 0.0, 1.05 * hist.GetMaximum())
    zero.SetLineColor(ROOT.kGray + 2)
    zero.SetLineStyle(2)
    zero.Draw()
    canvas.SaveAs(output_path)


def draw_two_pt_correlations(hist_global, hist_outer, output_path):
    hist_global.GetXaxis().SetTitle("upper p_{T}^{global} [GeV]")
    hist_global.GetYaxis().SetTitle("lower p_{T}^{global} [GeV]")
    hist_outer.GetXaxis().SetTitle("upper p_{T}^{outer} [GeV]")
    hist_outer.GetYaxis().SetTitle("lower p_{T}^{outer} [GeV]")
    canvas = ROOT.TCanvas("c_upper_lower_pt_global_outer", "", 1400, 650)
    canvas.Divide(2, 1)
    diagonals = []
    for pad_index, hist in enumerate((hist_global, hist_outer), start=1):
        canvas.cd(pad_index)
        ROOT.gPad.SetRightMargin(0.14)
        ROOT.gPad.SetLogz()
        hist.Draw("COLZ")
        diagonal = ROOT.TLine(0.0, 0.0, 500.0, 500.0)
        diagonal.SetLineColor(ROOT.kRed + 1)
        diagonal.SetLineStyle(2)
        diagonal.SetLineWidth(2)
        diagonal.Draw()
        diagonals.append(diagonal)
    canvas.SaveAs(output_path)


def draw_global_outer_side_correlations(hist_upper, hist_lower, output_path):
    hist_upper.GetXaxis().SetTitle("upper p_{T}^{global} [GeV]")
    hist_upper.GetYaxis().SetTitle("upper p_{T}^{outer} [GeV]")
    hist_lower.GetXaxis().SetTitle("lower p_{T}^{global} [GeV]")
    hist_lower.GetYaxis().SetTitle("lower p_{T}^{outer} [GeV]")
    canvas = ROOT.TCanvas("c_global_outer_pt_upper_lower", "", 1400, 650)
    canvas.Divide(2, 1)
    diagonals = []
    for pad_index, hist in enumerate((hist_upper, hist_lower), start=1):
        canvas.cd(pad_index)
        ROOT.gPad.SetRightMargin(0.14)
        ROOT.gPad.SetLogz()
        hist.Draw("COLZ")
        diagonal = ROOT.TLine(0.0, 0.0, 500.0, 500.0)
        diagonal.SetLineColor(ROOT.kRed + 1)
        diagonal.SetLineStyle(2)
        diagonal.SetLineWidth(2)
        diagonal.Draw()
        diagonals.append(diagonal)
    canvas.SaveAs(output_path)


def draw_asymmetry_correlation(hist, output_path):
    hist.GetXaxis().SetTitle("global-track upper/lower p_{T} asymmetry")
    hist.GetYaxis().SetTitle("outer-track upper/lower p_{T} asymmetry")
    canvas = ROOT.TCanvas("c_outer_vs_global_asymmetry", "", 900, 750)
    canvas.SetRightMargin(0.14)
    canvas.SetLogz()
    hist.Draw("COLZ")
    diagonal = ROOT.TLine(-2.0, -2.0, 2.0, 2.0)
    diagonal.SetLineColor(ROOT.kRed + 1)
    diagonal.SetLineStyle(2)
    diagonal.SetLineWidth(2)
    diagonal.Draw()
    zero_x = ROOT.TLine(0.0, -2.0, 0.0, 2.0)
    zero_y = ROOT.TLine(-2.0, 0.0, 2.0, 0.0)
    for line in (zero_x, zero_y):
        line.SetLineColor(ROOT.kGray + 2)
        line.SetLineStyle(3)
        line.Draw()
    canvas.SaveAs(output_path)


def fit_central_gaussian(hist, name):
    if hist.GetEntries() < 30 or hist.Integral() <= 0:
        return None

    probabilities = array("d", [0.16, 0.50, 0.84])
    quantiles = array("d", [0.0, 0.0, 0.0])
    hist.GetQuantiles(3, quantiles, probabilities)
    center = quantiles[1]
    width = 0.5 * (quantiles[2] - quantiles[0])
    if width <= 0:
        width = hist.GetRMS()
    if width <= 0:
        return None

    fit_low = max(-2.0, center - 2.0 * width)
    fit_high = min(2.0, center + 2.0 * width)
    function = ROOT.TF1(name, "gaus", fit_low, fit_high)
    function.SetParameters(hist.GetMaximum(), center, width)
    result = hist.Fit(function, "QRSN")
    if int(result) != 0 or function.GetParameter(2) <= 0:
        return None
    return {
        "mean": function.GetParameter(1),
        "mean_error": function.GetParError(1),
        "sigma": abs(function.GetParameter(2)),
        "sigma_error": function.GetParError(2),
        "function": function,
    }


def make_trend_graphs(chain, selection, bin_expression, bin_label, global_asym, outer_asym):
    pt_bins = [20.0, 30.0, 40.0, 50.0, 65.0, 85.0, 120.0, 200.0, 400.0, 1000.0]
    graphs = {
        "mean_global": ROOT.TGraphErrors(),
        "mean_outer": ROOT.TGraphErrors(),
        "sigma_global": ROOT.TGraphErrors(),
        "sigma_outer": ROOT.TGraphErrors(),
    }
    histograms = []
    functions = []
    point_index = 0
    for low, high in zip(pt_bins[:-1], pt_bins[1:]):
        bin_selection = f"({selection}) && ({bin_expression})>={low} && ({bin_expression})<{high}"
        h_global = make_hist(
            chain,
            f"h_global_asymmetry_pt_{low:g}_{high:g}",
            global_asym,
            bin_selection,
            160,
            -2.0,
            2.0,
        )
        h_outer = make_hist(
            chain,
            f"h_outer_asymmetry_pt_{low:g}_{high:g}",
            outer_asym,
            bin_selection,
            160,
            -2.0,
            2.0,
        )
        histograms.extend([h_global, h_outer])
        fit_global = fit_central_gaussian(h_global, f"f_global_asymmetry_pt_{low:g}_{high:g}")
        fit_outer = fit_central_gaussian(h_outer, f"f_outer_asymmetry_pt_{low:g}_{high:g}")
        if not fit_global or not fit_outer:
            print(f"Skipping trend bin {low:g}-{high:g} GeV: insufficient or failed fit")
            continue

        x = math.sqrt(low * high)
        x_error = 0.5 * (high - low)
        for key, fit in (("global", fit_global), ("outer", fit_outer)):
            graphs[f"mean_{key}"].SetPoint(point_index, x, fit["mean"])
            graphs[f"mean_{key}"].SetPointError(point_index, x_error, fit["mean_error"])
            graphs[f"sigma_{key}"].SetPoint(point_index, x, fit["sigma"])
            graphs[f"sigma_{key}"].SetPointError(point_index, x_error, fit["sigma_error"])
            functions.append(fit["function"])
        point_index += 1

    for key, graph in graphs.items():
        graph.SetName(f"g_{key}_vs_pt")
    return graphs, histograms, functions, bin_label


def graph_values(graph, include_errors=False):
    values = []
    for index in range(graph.GetN()):
        value = graph.GetPointY(index)
        if include_errors:
            error = graph.GetErrorY(index)
            values.extend([value - error, value + error])
        else:
            values.append(value)
    return values


def draw_trend(graph_global, graph_outer, ylabel, bin_label, output_path, is_sigma=False):
    graph_global.SetMarkerStyle(20)
    graph_outer.SetMarkerStyle(21)
    graph_global.SetMarkerColor(ROOT.kRed + 1)
    graph_outer.SetMarkerColor(ROOT.kBlue + 1)
    graph_global.SetLineColor(ROOT.kRed + 1)
    graph_outer.SetLineColor(ROOT.kBlue + 1)

    values = graph_values(graph_global, True) + graph_values(graph_outer, True)
    if not values:
        print(f"No valid points for {output_path}; skipping plot")
        return
    if is_sigma:
        ymin = 0.0
        ymax = max(values) * 1.25
    else:
        span = max(max(values) - min(values), 0.05)
        ymin = min(min(values), 0.0) - 0.20 * span
        ymax = max(max(values), 0.0) + 0.20 * span

    canvas = ROOT.TCanvas(f"c_{graph_global.GetName()}", "", 900, 700)
    canvas.SetLogx()
    frame = ROOT.TH1F(f"frame_{graph_global.GetName()}", "", 100, 20.0, 1000.0)
    frame.SetMinimum(ymin)
    frame.SetMaximum(ymax if ymax > ymin else ymin + 1.0)
    frame.GetXaxis().SetTitle(bin_label)
    frame.GetYaxis().SetTitle(ylabel)
    frame.Draw()
    graph_global.Draw("P E1 SAME")
    graph_outer.Draw("P E1 SAME")
    zero = None
    if not is_sigma:
        zero = ROOT.TLine(20.0, 0.0, 1000.0, 0.0)
        zero.SetLineColor(ROOT.kGray + 2)
        zero.SetLineStyle(2)
        zero.Draw()
    legend = ROOT.TLegend(0.62, 0.76, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(graph_global, "global-track fit", "lp")
    legend.AddEntry(graph_outer, "outer-track fit", "lp")
    legend.Draw()
    canvas.SaveAs(output_path)


def main():
    parser = argparse.ArgumentParser(
        description="Compare upper/lower DGL momentum using global and associated outer tracks"
    )
    parser.add_argument("--input", nargs="+", required=True, help="ROOT files or quoted glob patterns")
    parser.add_argument("--tree", default="Events", help="Input tree name")
    parser.add_argument("--outdir", default="dgl_global_outer_study", help="Output directory")
    parser.add_argument(
        "--selection",
        choices=("raw", "quality", "resolution"),
        default="resolution",
        help="Symmetric DGL pair selection used for both momentum estimators",
    )
    args = parser.parse_args()

    files = collect_files(args.input)
    if not files:
        raise RuntimeError("No ROOT files found from --input")
    os.makedirs(args.outdir, exist_ok=True)
    chain = ROOT.TChain(args.tree)
    for path in files:
        chain.Add(path)
    if chain.GetEntries() == 0:
        raise RuntimeError("The input chain contains no events")

    required = [
        "ndmu",
        "dmu_isDGL",
        "dmu_dgl_pt",
        "dmu_dgl_eta",
        "dmu_dgl_phi",
        "dmu_dgl_ptError",
        "dmu_dgl_nMuonHits",
        "dmu_dgl_nValidStripHits",
    ]
    has_explicit_outer = bool(
        chain.GetBranch("dmu_dgl_hasOuterTrack")
        and chain.GetBranch("dmu_dgl_outer_pt")
        and chain.GetBranch("dmu_dgl_outer_side")
    )
    if has_explicit_outer:
        required.extend(
            ["dmu_dgl_hasOuterTrack", "dmu_dgl_outer_pt", "dmu_dgl_outer_side"]
        )
    else:
        required.extend(["dmu_isDSA", "dmu_dsa_pt", "dmu_dsa_side"])
    require_branches(chain, required)

    if has_explicit_outer:
        outer_branch = "dmu_dgl_outer_pt"
        side_branch = "dmu_dgl_outer_side"
        candidate = (
            "dmu_isDGL && dmu_dgl_hasOuterTrack && "
            "dmu_dgl_pt>0 && dmu_dgl_outer_pt>0"
        )
        print("Using explicit DGL-associated outer-track branches")
    else:
        outer_branch = "dmu_dsa_pt"
        side_branch = "dmu_dsa_side"
        candidate = (
            "dmu_isDGL && dmu_isDSA && "
            "dmu_dgl_pt>0 && dmu_dsa_pt>0"
        )
        print("Using legacy same-index dmu_dsa_pt as the DGL outer-track momentum")

    # MiniAOD may not retain the TrackExtra positions needed by the geometric side
    # branch. For collision-like DGL tracks, the y direction is given by sin(phi).
    side_expression = (
        f"({side_branch}) + (({side_branch})==0)*sin(dmu_dgl_phi)"
    )
    upper_candidate = f"({candidate}) && ({side_expression})>0"
    lower_candidate = f"({candidate}) && ({side_expression})<0"
    global_upper = f"Sum$(dmu_dgl_pt*({upper_candidate}))"
    global_lower = f"Sum$(dmu_dgl_pt*({lower_candidate}))"
    outer_upper = f"Sum$({outer_branch}*({upper_candidate}))"
    outer_lower = f"Sum$({outer_branch}*({lower_candidate}))"
    eta_upper = f"Sum$(dmu_dgl_eta*({upper_candidate}))"
    eta_lower = f"Sum$(dmu_dgl_eta*({lower_candidate}))"
    phi_upper = f"Sum$(dmu_dgl_phi*({upper_candidate}))"
    phi_lower = f"Sum$(dmu_dgl_phi*({lower_candidate}))"
    cos_alpha = (
        f"tanh({eta_upper})*tanh({eta_lower}) + "
        f"cos(({phi_upper})-({phi_lower}))/(cosh({eta_upper})*cosh({eta_lower}))"
    )

    raw_selection = (
        f"Sum$({candidate})==2 && "
        f"Sum$({upper_candidate})==1 && Sum$({lower_candidate})==1 && "
        f"({cos_alpha})<{math.cos(2.8):.9f}"
    )
    quality_candidate = (
        f"({candidate}) && dmu_dgl_pt>20 && abs(dmu_dgl_eta)<0.9 && "
        "dmu_dgl_nMuonHits>12 && dmu_dgl_nValidStripHits>5"
    )
    resolution_candidate = (
        f"({quality_candidate}) && dmu_dgl_ptError/dmu_dgl_pt<0.3"
    )
    selection = {
        "raw": raw_selection,
        "quality": f"({raw_selection}) && Sum$({quality_candidate})==2",
        "resolution": f"({raw_selection}) && Sum$({resolution_candidate})==2",
    }[args.selection]

    global_asymmetry = f"2.0*({global_lower}-{global_upper})/({global_lower}+{global_upper})"
    outer_asymmetry = f"2.0*({outer_lower}-{outer_upper})/({outer_lower}+{outer_upper})"

    print(f"Input files: {len(files)}")
    print(f"Input events: {chain.GetEntries()}")

    # Read remote files once, then make all plots from a compact local selected tree.
    selected_path = os.path.join(args.outdir, "selected_dgl_pairs.root")
    chain.SetBranchStatus("*", 0)
    for branch_name in required:
        chain.SetBranchStatus(branch_name, 1)
    selected_file = ROOT.TFile(selected_path, "RECREATE")
    selected_tree = chain.CopyTree(selection)
    selected_count = selected_tree.GetEntries()
    selected_tree.Write()
    selected_file.Close()
    print(f"Selected upper/lower DGL pairs ({args.selection}): {selected_count}")

    if selected_count == 0:
        print_zero_selection_diagnostics(
            files[0],
            args.tree,
            candidate,
            upper_candidate,
            lower_candidate,
            raw_selection,
            quality_candidate,
            resolution_candidate,
        )
        raise RuntimeError(
            "No upper/lower DGL pairs passed the selection; no PNG files were produced"
        )

    chain = ROOT.TChain(args.tree)
    chain.Add(selected_path)
    selection = "1"

    objects = []
    h_global_correlation = make_hist2d(
        chain,
        "h_upper_lower_pt_global",
        f"{global_lower}:{global_upper}",
        selection,
        100,
        0.0,
        500.0,
        100,
        0.0,
        500.0,
    )
    h_outer_correlation = make_hist2d(
        chain,
        "h_upper_lower_pt_outer",
        f"{outer_lower}:{outer_upper}",
        selection,
        100,
        0.0,
        500.0,
        100,
        0.0,
        500.0,
    )
    objects.extend([h_global_correlation, h_outer_correlation])
    draw_two_pt_correlations(
        h_global_correlation,
        h_outer_correlation,
        os.path.join(args.outdir, "upper_lower_pt_global_outer.png"),
    )

    h_upper_global_outer = make_hist2d(
        chain,
        "h_upper_global_outer_pt",
        f"{outer_upper}:{global_upper}",
        selection,
        100,
        0.0,
        500.0,
        100,
        0.0,
        500.0,
    )
    h_lower_global_outer = make_hist2d(
        chain,
        "h_lower_global_outer_pt",
        f"{outer_lower}:{global_lower}",
        selection,
        100,
        0.0,
        500.0,
        100,
        0.0,
        500.0,
    )
    objects.extend([h_upper_global_outer, h_lower_global_outer])
    draw_global_outer_side_correlations(
        h_upper_global_outer,
        h_lower_global_outer,
        os.path.join(args.outdir, "global_vs_outer_pt_upper_lower.png"),
    )

    h_upper_global_pt = make_hist(
        chain, "h_upper_global_pt", global_upper, selection, 100, 0.0, 500.0
    )
    h_upper_outer_pt = make_hist(
        chain, "h_upper_outer_pt", outer_upper, selection, 100, 0.0, 500.0
    )
    h_lower_global_pt = make_hist(
        chain, "h_lower_global_pt", global_lower, selection, 100, 0.0, 500.0
    )
    h_lower_outer_pt = make_hist(
        chain, "h_lower_outer_pt", outer_lower, selection, 100, 0.0, 500.0
    )
    objects.extend(
        [h_upper_global_pt, h_upper_outer_pt, h_lower_global_pt, h_lower_outer_pt]
    )
    draw_overlay(
        h_upper_global_pt,
        h_upper_outer_pt,
        "upper DGL p_{T} [GeV]",
        os.path.join(args.outdir, "pt_upper_global_vs_outer.png"),
    )
    draw_overlay(
        h_lower_global_pt,
        h_lower_outer_pt,
        "lower DGL p_{T} [GeV]",
        os.path.join(args.outdir, "pt_lower_global_vs_outer.png"),
    )
    draw_overlay(
        h_upper_global_pt,
        h_lower_global_pt,
        "DGL global-track p_{T} [GeV]",
        os.path.join(args.outdir, "pt_upper_vs_lower_global.png"),
        "upper DGL",
        "lower DGL",
    )
    draw_overlay(
        h_upper_outer_pt,
        h_lower_outer_pt,
        "DGL outer-track p_{T} [GeV]",
        os.path.join(args.outdir, "pt_upper_vs_lower_outer.png"),
        "upper DGL",
        "lower DGL",
    )

    h_global_asymmetry = make_hist(
        chain, "h_global_upper_lower_pt_asymmetry", global_asymmetry, selection, 160, -2.0, 2.0
    )
    h_outer_asymmetry = make_hist(
        chain, "h_outer_upper_lower_pt_asymmetry", outer_asymmetry, selection, 160, -2.0, 2.0
    )
    objects.extend([h_global_asymmetry, h_outer_asymmetry])
    draw_overlay(
        h_global_asymmetry,
        h_outer_asymmetry,
        "2(p_{T}^{lower}-p_{T}^{upper})/(p_{T}^{lower}+p_{T}^{upper})",
        os.path.join(args.outdir, "upper_lower_pt_asymmetry_global_outer.png"),
    )

    outer_inverse_pt_residual = (
        f"((1.0/({outer_upper}))-(1.0/({outer_lower})))/(1.0/({outer_lower}))"
    )
    h_outer_inverse_pt_residual = make_hist(
        chain,
        "h_outer_inverse_pt_relative_residual",
        outer_inverse_pt_residual,
        selection,
        200,
        -5.0,
        5.0,
    )
    objects.append(h_outer_inverse_pt_residual)
    draw_distribution(
        h_outer_inverse_pt_residual,
        "[(1/p_{T}^{upper})-(1/p_{T}^{lower})]/(1/p_{T}^{lower})",
        os.path.join(args.outdir, "outer_inverse_pt_relative_residual.png"),
    )

    global_inverse_pt_residual = (
        f"((1.0/({global_upper}))-(1.0/({global_lower})))/(1.0/({global_lower}))"
    )
    h_global_inverse_pt_residual = make_hist(
        chain,
        "h_global_inverse_pt_relative_residual",
        global_inverse_pt_residual,
        selection,
        200,
        -5.0,
        5.0,
    )
    objects.append(h_global_inverse_pt_residual)
    draw_distribution(
        h_global_inverse_pt_residual,
        "[(1/p_{T}^{upper})-(1/p_{T}^{lower})]/(1/p_{T}^{lower})",
        os.path.join(args.outdir, "global_inverse_pt_relative_residual.png"),
    )

    h_asymmetry_correlation = make_hist2d(
        chain,
        "h_outer_vs_global_pt_asymmetry",
        f"({outer_asymmetry}):({global_asymmetry})",
        selection,
        120,
        -2.0,
        2.0,
        120,
        -2.0,
        2.0,
    )
    objects.append(h_asymmetry_correlation)
    draw_asymmetry_correlation(
        h_asymmetry_correlation,
        os.path.join(args.outdir, "outer_vs_global_pt_asymmetry.png"),
    )

    h_upper_outer_over_global = make_hist(
        chain,
        "h_upper_outer_over_global_pt",
        f"{outer_upper}/{global_upper}",
        selection,
        120,
        0.0,
        3.0,
    )
    h_lower_outer_over_global = make_hist(
        chain,
        "h_lower_outer_over_global_pt",
        f"{outer_lower}/{global_lower}",
        selection,
        120,
        0.0,
        3.0,
    )
    objects.extend([h_upper_outer_over_global, h_lower_outer_over_global])
    draw_overlay(
        h_upper_outer_over_global,
        h_lower_outer_over_global,
        "p_{T}^{outer}/p_{T}^{global}",
        os.path.join(args.outdir, "outer_over_global_pt_upper_lower.png"),
    )

    common_average_pt = (
        f"0.25*({global_upper}+{global_lower}+{outer_upper}+{outer_lower})"
    )
    trend_selection = selection
    bin_expression = common_average_pt
    bin_label = "common average reconstructed p_{T} [GeV]"

    graphs, trend_hists, functions, bin_label = make_trend_graphs(
        chain,
        trend_selection,
        bin_expression,
        bin_label,
        global_asymmetry,
        outer_asymmetry,
    )
    objects.extend(trend_hists)
    objects.extend(functions)
    objects.extend(graphs.values())
    draw_trend(
        graphs["mean_global"],
        graphs["mean_outer"],
        "Gaussian mean of upper/lower p_{T} asymmetry",
        bin_label,
        os.path.join(args.outdir, "asymmetry_mean_global_outer_vs_pt.png"),
    )
    draw_trend(
        graphs["sigma_global"],
        graphs["sigma_outer"],
        "Gaussian #sigma of upper/lower p_{T} asymmetry",
        bin_label,
        os.path.join(args.outdir, "asymmetry_sigma_global_outer_vs_pt.png"),
        is_sigma=True,
    )

    output = ROOT.TFile(os.path.join(args.outdir, "dgl_global_outer_study.root"), "RECREATE")
    for obj in objects:
        if obj:
            obj.Write()
    output.Close()
    print(f"Wrote DGL global-vs-outer plots and ROOT objects to {args.outdir}")


if __name__ == "__main__":
    main()
