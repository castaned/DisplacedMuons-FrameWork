#!/usr/bin/env python3

import argparse
import glob
import os

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
        raise RuntimeError("Missing DSA study branches: " + ", ".join(missing))


def make_hist(chain, name, expression, selection, bins, xmin, xmax):
    hist = ROOT.TH1F(name, "", bins, xmin, xmax)
    chain.Draw(f"{expression}>>{name}", selection, "goff")
    hist.SetDirectory(0)
    return hist


def normalize(hist):
    integral = hist.Integral(0, hist.GetNbinsX() + 1)
    if integral > 0:
        hist.Scale(1.0 / integral)


def draw_overlay(hist_a, hist_b, label_a, label_b, xtitle, output_path):
    normalize(hist_a)
    normalize(hist_b)
    hist_a.SetLineColor(ROOT.kRed + 1)
    hist_b.SetLineColor(ROOT.kBlue + 1)
    hist_a.SetLineWidth(2)
    hist_b.SetLineWidth(2)
    ymax = max(hist_a.GetMaximum(), hist_b.GetMaximum()) * 1.25
    hist_a.SetMaximum(ymax if ymax > 0 else 1.0)
    hist_a.SetMinimum(0.0)
    hist_a.GetXaxis().SetTitle(xtitle)
    hist_a.GetYaxis().SetTitle("Normalized events")

    canvas = ROOT.TCanvas(f"c_{hist_a.GetName()}", "", 900, 700)
    hist_a.Draw("HIST")
    hist_b.Draw("HIST SAME")
    legend = ROOT.TLegend(0.66, 0.76, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(hist_a, label_a, "l")
    legend.AddEntry(hist_b, label_b, "l")
    legend.Draw()
    canvas.SaveAs(output_path)


def draw_single(hist, xtitle, output_path):
    hist.SetLineColor(ROOT.kBlue + 1)
    hist.SetLineWidth(2)
    hist.GetXaxis().SetTitle(xtitle)
    hist.GetYaxis().SetTitle("Events")
    canvas = ROOT.TCanvas(f"c_{hist.GetName()}", "", 900, 700)
    hist.Draw("HIST")
    canvas.SaveAs(output_path)


def draw_2d(chain, expression, selection, output_path):
    hist = ROOT.TH2F("h_pt_upper_vs_lower", "", 80, 0.0, 400.0, 80, 0.0, 400.0)
    chain.Draw(f"{expression}>>h_pt_upper_vs_lower", selection, "goff")
    hist.SetDirectory(0)
    hist.GetXaxis().SetTitle("upper DSA p_{T} [GeV]")
    hist.GetYaxis().SetTitle("lower DSA p_{T} [GeV]")
    canvas = ROOT.TCanvas("c_pt_upper_vs_lower", "", 900, 750)
    canvas.SetRightMargin(0.14)
    hist.Draw("COLZ")
    diagonal = ROOT.TLine(0.0, 0.0, 400.0, 400.0)
    diagonal.SetLineColor(ROOT.kRed + 1)
    diagonal.SetLineStyle(2)
    diagonal.Draw()
    canvas.SaveAs(output_path)
    return hist


def threshold_scan(chain, thresholds, base_selection, outdir):
    mean_graph = ROOT.TGraphErrors(len(thresholds))
    count_graph = ROOT.TGraph(len(thresholds))
    mean_graph.SetName("g_dsa_residual_mean_vs_min_pt")
    count_graph.SetName("g_dsa_event_count_vs_min_pt")
    scan_hists = []

    for index, threshold in enumerate(thresholds):
        selection = (
            f"({base_selection}) && evt_dsa_pt_upper>{threshold} && "
            f"evt_dsa_pt_lower>{threshold} && abs(evt_dsa_residual_lower_upper)<5"
        )
        hist = make_hist(
            chain,
            f"h_residual_minpt_{threshold:g}",
            "evt_dsa_residual_lower_upper",
            selection,
            200,
            -5.0,
            5.0,
        )
        scan_hists.append(hist)
        entries = hist.GetEntries()
        mean = hist.GetMean() if entries else 0.0
        mean_error = hist.GetMeanError() if entries else 0.0
        mean_graph.SetPoint(index, threshold, mean)
        mean_graph.SetPointError(index, 0.0, mean_error)
        count_graph.SetPoint(index, threshold, entries)

    canvas = ROOT.TCanvas("c_threshold_mean", "", 900, 700)
    mean_graph.SetMarkerStyle(20)
    mean_graph.SetMarkerColor(ROOT.kBlue + 1)
    mean_graph.SetLineColor(ROOT.kBlue + 1)
    mean_graph.SetTitle(";common minimum p_{T} [GeV];mean lower/upper q/p_{T} residual")
    mean_graph.Draw("APL")
    canvas.Update()
    zero = ROOT.TLine(
        thresholds[0],
        0.0,
        thresholds[-1],
        0.0,
    )
    zero.SetLineStyle(2)
    zero.SetLineColor(ROOT.kGray + 2)
    zero.Draw()
    canvas.SaveAs(os.path.join(outdir, "threshold_scan_residual_mean.png"))

    canvas_count = ROOT.TCanvas("c_threshold_count", "", 900, 700)
    count_graph.SetMarkerStyle(20)
    count_graph.SetMarkerColor(ROOT.kBlue + 1)
    count_graph.SetLineColor(ROOT.kBlue + 1)
    count_graph.SetTitle(";common minimum p_{T} [GeV];selected events")
    count_graph.Draw("APL")
    canvas_count.SaveAs(os.path.join(outdir, "threshold_scan_event_count.png"))
    return mean_graph, count_graph, scan_hists


def main():
    parser = argparse.ArgumentParser(description="Plot the exactly-two-DSA ntuplizer bias study")
    parser.add_argument("--input", nargs="+", required=True, help="ROOT files or quoted glob patterns")
    parser.add_argument("--tree", default="Events", help="Input tree name")
    parser.add_argument("--outdir", default="dsa_bias_study", help="Output directory")
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
        "evt_dsa_nReco",
        "event",
        "evt_dsa_pt_1",
        "evt_dsa_pt_2",
        "evt_dsa_pt_asymmetry",
        "evt_dsa_absqoverpt_1",
        "evt_dsa_absqoverpt_2",
        "evt_dsa_qoverpt_1",
        "evt_dsa_qoverpt_2",
        "evt_dsa_residual_12",
        "evt_dsa_residual_21",
        "evt_dsa_cosAlpha_12",
        "evt_dsa_oppositeSides",
        "evt_dsa_index_upper",
        "evt_dsa_index_lower",
        "evt_dsa_pt_upper",
        "evt_dsa_pt_lower",
        "evt_dsa_residual_lower_upper",
        "dmu_dsa_pt",
        "dmu_dsa_eta",
        "dmu_dsa_ptError",
        "dmu_dsa_normalizedChi2",
        "dmu_dsa_nValidMuonDTHits",
    ]
    require_branches(chain, required)

    exactly_two_selection = "evt_dsa_nReco==2"
    selection_flags = [
        "evt_dsa_passRawPair",
        "evt_dsa_passQualityPair",
        "evt_dsa_passResolutionPair",
    ]
    if all(chain.GetBranch(name) for name in selection_flags):
        raw_selection = "evt_dsa_passRawPair"
        quality_selection = "evt_dsa_passQualityPair"
        resolution_selection = "evt_dsa_passResolutionPair"
        print("Using stored DSA pair-selection flags")
    else:
        upper = "evt_dsa_index_upper"
        lower = "evt_dsa_index_lower"
        raw_selection = (
            "evt_dsa_nReco==2 && evt_dsa_oppositeSides && "
            "evt_dsa_cosAlpha_12<-0.5048461"
        )
        quality_selection = (
            f"({raw_selection}) && "
            f"dmu_dsa_pt[{upper}]>20 && dmu_dsa_pt[{lower}]>20 && "
            f"abs(dmu_dsa_eta[{upper}])<0.7 && abs(dmu_dsa_eta[{lower}])<0.7 && "
            f"dmu_dsa_nValidMuonDTHits[{upper}]>=31 && "
            f"dmu_dsa_nValidMuonDTHits[{lower}]>=31 && "
            f"dmu_dsa_normalizedChi2[{upper}]<5 && dmu_dsa_normalizedChi2[{lower}]<5"
        )
        resolution_selection = (
            f"({quality_selection}) && "
            f"dmu_dsa_ptError[{upper}]/dmu_dsa_pt[{upper}]<0.5 && "
            f"dmu_dsa_ptError[{lower}]/dmu_dsa_pt[{lower}]<0.5"
        )
        print("Stored pair-selection flags not found; rebuilding equivalent selections from branches")
    exactly_two_count = chain.GetEntries(exactly_two_selection)
    raw_count = chain.GetEntries(raw_selection)
    quality_count = chain.GetEntries(quality_selection)
    resolution_count = chain.GetEntries(resolution_selection)
    print(f"Exactly-two-DSA events: {exactly_two_count}")
    print(f"Raw opposite-side angular pairs: {raw_count}")
    print(f"Quality-matched pairs: {quality_count}")
    print(f"Resolution pairs with common ptError/pt cut: {resolution_count}")
    if raw_count == 0:
        print("WARNING: no upper/lower pairs found; check TrackExtra availability and dmu_dsa_side")
    objects = []

    h_pt_1 = make_hist(chain, "h_pt_1", "evt_dsa_pt_1", raw_selection, 100, 0.0, 400.0)
    h_pt_2 = make_hist(chain, "h_pt_2", "evt_dsa_pt_2", raw_selection, 100, 0.0, 400.0)
    objects.extend([h_pt_1, h_pt_2])
    draw_overlay(
        h_pt_1,
        h_pt_2,
        "collection index 1",
        "collection index 2",
        "DSA p_{T} [GeV]",
        os.path.join(args.outdir, "pt_collection_order_comparison.png"),
    )

    h_pt_upper = make_hist(chain, "h_pt_upper", "evt_dsa_pt_upper", quality_selection, 100, 0.0, 400.0)
    h_pt_lower = make_hist(chain, "h_pt_lower", "evt_dsa_pt_lower", quality_selection, 100, 0.0, 400.0)
    objects.extend([h_pt_upper, h_pt_lower])
    draw_overlay(
        h_pt_upper,
        h_pt_lower,
        "upper detector side",
        "lower detector side",
        "DSA p_{T} [GeV]",
        os.path.join(args.outdir, "pt_upper_lower_comparison.png"),
    )

    quality_plots = [
        (
            "pt_error_over_pt",
            "dmu_dsa_ptError[evt_dsa_index_upper]/dmu_dsa_pt[evt_dsa_index_upper]",
            "dmu_dsa_ptError[evt_dsa_index_lower]/dmu_dsa_pt[evt_dsa_index_lower]",
            100,
            0.0,
            1.0,
            "p_{T}^{error}/p_{T}",
        ),
        (
            "normalized_chi2",
            "dmu_dsa_normalizedChi2[evt_dsa_index_upper]",
            "dmu_dsa_normalizedChi2[evt_dsa_index_lower]",
            100,
            0.0,
            10.0,
            "normalized #chi^{2}",
        ),
        (
            "valid_dt_hits",
            "dmu_dsa_nValidMuonDTHits[evt_dsa_index_upper]",
            "dmu_dsa_nValidMuonDTHits[evt_dsa_index_lower]",
            81,
            -0.5,
            80.5,
            "valid DT hits",
        ),
    ]
    for name, upper_expr, lower_expr, bins, xmin, xmax, xtitle in quality_plots:
        h_upper = make_hist(chain, f"h_upper_{name}", upper_expr, quality_selection, bins, xmin, xmax)
        h_lower = make_hist(chain, f"h_lower_{name}", lower_expr, quality_selection, bins, xmin, xmax)
        objects.extend([h_upper, h_lower])
        draw_overlay(
            h_upper,
            h_lower,
            "upper detector side",
            "lower detector side",
            xtitle,
            os.path.join(args.outdir, f"{name}_upper_lower_comparison.png"),
        )

    qoverpt_selection = (
        f"({quality_selection}) && evt_dsa_absqoverpt_1>0 && evt_dsa_absqoverpt_2>0"
    )
    h_qoverpt_consistency_1 = make_hist(
        chain,
        "h_qoverpt_consistency_1",
        "(abs(evt_dsa_qoverpt_1)-evt_dsa_absqoverpt_1)/evt_dsa_absqoverpt_1",
        qoverpt_selection,
        100,
        -0.01,
        0.01,
    )
    h_qoverpt_consistency_2 = make_hist(
        chain,
        "h_qoverpt_consistency_2",
        "(abs(evt_dsa_qoverpt_2)-evt_dsa_absqoverpt_2)/evt_dsa_absqoverpt_2",
        qoverpt_selection,
        100,
        -0.01,
        0.01,
    )
    objects.extend([h_qoverpt_consistency_1, h_qoverpt_consistency_2])
    draw_overlay(
        h_qoverpt_consistency_1,
        h_qoverpt_consistency_2,
        "collection index 1",
        "collection index 2",
        "relative difference: |q/p| cosh(#eta) vs |q|/p_{T}",
        os.path.join(args.outdir, "qoverpt_parameterization_consistency.png"),
    )

    h_residual_12 = make_hist(
        chain, "h_residual_12", "evt_dsa_residual_12", resolution_selection, 200, -5.0, 5.0
    )
    h_residual_21 = make_hist(
        chain, "h_residual_21", "evt_dsa_residual_21", resolution_selection, 200, -5.0, 5.0
    )
    objects.extend([h_residual_12, h_residual_21])
    draw_overlay(
        h_residual_12,
        h_residual_21,
        "index 1 as reference",
        "index 2 as reference",
        "signed q/p_{T} relative residual",
        os.path.join(args.outdir, "residual_order_swap_comparison.png"),
    )

    h_residual_random = make_hist(
        chain,
        "h_residual_random_order",
        "event%2==0 ? evt_dsa_residual_12 : evt_dsa_residual_21",
        resolution_selection,
        200,
        -5.0,
        5.0,
    )
    objects.append(h_residual_random)
    draw_single(
        h_residual_random,
        "signed q/p_{T} relative residual with deterministic random ordering",
        os.path.join(args.outdir, "residual_randomized_order.png"),
    )

    h_asymmetry = make_hist(
        chain, "h_abs_pt_asymmetry", "abs(evt_dsa_pt_asymmetry)", resolution_selection, 100, 0.0, 1.0
    )
    h_side_residual = make_hist(
        chain,
        "h_residual_lower_upper",
        "evt_dsa_residual_lower_upper",
        resolution_selection,
        200,
        -5.0,
        5.0,
    )
    objects.extend([h_asymmetry, h_side_residual])
    draw_single(h_asymmetry, "|p_{T,2}-p_{T,1}|/(p_{T,2}+p_{T,1})", os.path.join(args.outdir, "pt_asymmetry.png"))
    draw_single(
        h_side_residual,
        "(q/p_{T})_{lower}-(q/p_{T})_{upper} over (q/p_{T})_{upper}",
        os.path.join(args.outdir, "residual_lower_upper.png"),
    )

    h_pt_2d = draw_2d(
        chain,
        "evt_dsa_pt_lower:evt_dsa_pt_upper",
        quality_selection,
        os.path.join(args.outdir, "pt_upper_vs_lower.png"),
    )
    objects.append(h_pt_2d)

    mean_graph, count_graph, scan_hists = threshold_scan(
        chain, [20.0, 30.0, 40.0, 50.0], resolution_selection, args.outdir
    )
    objects.extend(scan_hists)
    objects.extend([mean_graph, count_graph])

    output_root = ROOT.TFile(os.path.join(args.outdir, "dsa_bias_study.root"), "RECREATE")
    for obj in objects:
        obj.Write()
    output_root.Close()
    print(f"Wrote DSA bias-study plots and ROOT objects to {args.outdir}")


if __name__ == "__main__":
    main()
