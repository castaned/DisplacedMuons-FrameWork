#!/usr/bin/env python3

import argparse
import glob
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


def draw_pair_momentum_diagnostics(chain, selection, outdir):
    objects = []

    zoom = ROOT.TH2F("h_pt_upper_vs_lower_zoom", "", 65, 20.0, 150.0, 65, 20.0, 150.0)
    chain.Draw("evt_dsa_pt_lower:evt_dsa_pt_upper>>h_pt_upper_vs_lower_zoom", selection, "goff")
    zoom.SetDirectory(0)
    zoom.GetXaxis().SetTitle("upper DSA p_{T} [GeV]")
    zoom.GetYaxis().SetTitle("lower DSA p_{T} [GeV]")
    canvas_zoom = ROOT.TCanvas("c_pt_upper_vs_lower_zoom", "", 900, 750)
    canvas_zoom.SetRightMargin(0.14)
    canvas_zoom.SetLogz()
    zoom.Draw("COLZ")
    diagonal_zoom = ROOT.TLine(20.0, 20.0, 150.0, 150.0)
    diagonal_zoom.SetLineColor(ROOT.kRed + 1)
    diagonal_zoom.SetLineStyle(2)
    diagonal_zoom.SetLineWidth(2)
    diagonal_zoom.Draw()
    canvas_zoom.SaveAs(os.path.join(outdir, "pt_upper_vs_lower_zoom_logz.png"))
    objects.append(zoom)

    profile = ROOT.TProfile("p_mean_lower_pt_vs_upper_pt", "", 65, 20.0, 150.0)
    chain.Draw("evt_dsa_pt_lower:evt_dsa_pt_upper>>p_mean_lower_pt_vs_upper_pt", selection, "goff")
    profile.SetDirectory(0)
    profile.SetMarkerStyle(20)
    profile.SetMarkerColor(ROOT.kBlue + 1)
    profile.SetLineColor(ROOT.kBlue + 1)
    profile.GetXaxis().SetTitle("upper DSA p_{T} [GeV]")
    profile.GetYaxis().SetTitle("mean lower DSA p_{T} [GeV]")
    profile.SetMinimum(20.0)
    profile.SetMaximum(150.0)
    canvas_profile = ROOT.TCanvas("c_mean_lower_pt_vs_upper_pt", "", 900, 700)
    profile.Draw("E1")
    diagonal_profile = ROOT.TLine(20.0, 20.0, 150.0, 150.0)
    diagonal_profile.SetLineColor(ROOT.kRed + 1)
    diagonal_profile.SetLineStyle(2)
    diagonal_profile.SetLineWidth(2)
    diagonal_profile.Draw()
    canvas_profile.SaveAs(os.path.join(outdir, "mean_lower_pt_vs_upper_pt_profile.png"))
    objects.append(profile)

    ratio = ROOT.TH2F("h_lower_over_upper_pt_vs_upper_pt", "", 80, 20.0, 400.0, 90, 0.0, 3.0)
    chain.Draw(
        "(evt_dsa_pt_lower/evt_dsa_pt_upper):evt_dsa_pt_upper"
        ">>h_lower_over_upper_pt_vs_upper_pt",
        selection,
        "goff",
    )
    ratio.SetDirectory(0)
    ratio.GetXaxis().SetTitle("upper DSA p_{T} [GeV]")
    ratio.GetYaxis().SetTitle("lower p_{T} / upper p_{T}")
    ratio_profile = ratio.ProfileX("p_lower_over_upper_pt_vs_upper_pt")
    ratio_profile.SetDirectory(0)
    ratio_profile.SetMarkerStyle(20)
    ratio_profile.SetMarkerColor(ROOT.kRed + 1)
    ratio_profile.SetLineColor(ROOT.kRed + 1)
    canvas_ratio = ROOT.TCanvas("c_lower_over_upper_pt_vs_upper_pt", "", 900, 750)
    canvas_ratio.SetRightMargin(0.14)
    canvas_ratio.SetLogz()
    ratio.Draw("COLZ")
    unity = ROOT.TLine(20.0, 1.0, 400.0, 1.0)
    unity.SetLineColor(ROOT.kBlack)
    unity.SetLineStyle(2)
    unity.Draw()
    ratio_profile.Draw("E1 SAME")
    canvas_ratio.SaveAs(os.path.join(outdir, "lower_over_upper_pt_vs_upper_pt.png"))
    objects.extend([ratio, ratio_profile])

    symmetric = ROOT.TH2F(
        "h_symmetric_residual_vs_average_pt", "", 80, 20.0, 400.0, 100, -2.0, 2.0
    )
    average_pt = "0.5*(evt_dsa_pt_upper+evt_dsa_pt_lower)"
    symmetric_residual = (
        "2.0*((1.0/evt_dsa_pt_lower)-(1.0/evt_dsa_pt_upper))"
        "/((1.0/evt_dsa_pt_lower)+(1.0/evt_dsa_pt_upper))"
    )
    chain.Draw(
        f"({symmetric_residual}):({average_pt})>>h_symmetric_residual_vs_average_pt",
        selection,
        "goff",
    )
    symmetric.SetDirectory(0)
    symmetric.GetXaxis().SetTitle("average upper/lower p_{T} [GeV]")
    symmetric.GetYaxis().SetTitle("symmetric |q|/p_{T} residual")
    symmetric_profile = symmetric.ProfileX("p_symmetric_residual_vs_average_pt")
    symmetric_profile.SetDirectory(0)
    symmetric_profile.SetMarkerStyle(20)
    symmetric_profile.SetMarkerColor(ROOT.kRed + 1)
    symmetric_profile.SetLineColor(ROOT.kRed + 1)
    canvas_symmetric = ROOT.TCanvas("c_symmetric_residual_vs_average_pt", "", 900, 750)
    canvas_symmetric.SetRightMargin(0.14)
    canvas_symmetric.SetLogz()
    symmetric.Draw("COLZ")
    zero = ROOT.TLine(20.0, 0.0, 400.0, 0.0)
    zero.SetLineColor(ROOT.kBlack)
    zero.SetLineStyle(2)
    zero.Draw()
    symmetric_profile.Draw("E1 SAME")
    canvas_symmetric.SaveAs(os.path.join(outdir, "symmetric_residual_vs_average_pt.png"))
    objects.extend([symmetric, symmetric_profile])
    return objects


def fit_single_gaussian(hist, name):
    """Fit one Gaussian model in a robust central range."""
    if hist.GetEntries() < 20 or hist.Integral() <= 0:
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

    fit_low = max(hist.GetXaxis().GetXmin(), center - 2.0 * width)
    fit_high = min(hist.GetXaxis().GetXmax(), center + 2.0 * width)
    gaussian = ROOT.TF1(name, "gaus", fit_low, fit_high)
    gaussian.SetParameters(hist.GetMaximum(), center, width)
    gaussian.SetLineColor(ROOT.kRed + 1)
    gaussian.SetLineWidth(2)
    fit_result = hist.Fit(gaussian, "RQ0S")
    if int(fit_result) != 0 or gaussian.GetParameter(2) <= 0:
        return None
    return gaussian


def draw_residual_grid(histograms, fits, pt_bins, output_path):
    canvas = ROOT.TCanvas("c_upper_lower_residual_by_upper_pt", "", 1350, 1050)
    canvas.Divide(3, 3, 0.002, 0.002)
    labels = []
    for index, hist in enumerate(histograms):
        canvas.cd(index + 1)
        ROOT.gPad.SetLeftMargin(0.14)
        ROOT.gPad.SetBottomMargin(0.13)
        hist.SetLineColor(ROOT.kBlue + 1)
        hist.SetLineWidth(2)
        hist.SetTitle(f"{pt_bins[index]:g} < upper p_{{T}} < {pt_bins[index + 1]:g} GeV")
        hist.GetXaxis().SetTitle("upper/lower |q|/p_{T} relative residual")
        hist.GetYaxis().SetTitle("Events")
        hist.Draw("HIST")
        if fits[index]:
            fits[index].Draw("SAME")
            label = ROOT.TLatex()
            label.SetNDC()
            label.SetTextSize(0.045)
            label.DrawLatex(
                0.53,
                0.78,
                f"#mu = {fits[index].GetParameter(1):.3f} #pm {fits[index].GetParError(1):.3f}",
            )
            label.DrawLatex(
                0.53,
                0.70,
                f"#sigma = {fits[index].GetParameter(2):.3f} #pm {fits[index].GetParError(2):.3f}",
            )
            labels.append(label)
        else:
            label = ROOT.TLatex()
            label.SetNDC()
            label.SetTextColor(ROOT.kRed + 1)
            label.SetTextSize(0.05)
            label.DrawLatex(0.48, 0.76, "fit unavailable")
            labels.append(label)
    canvas.SaveAs(output_path)


def draw_summary_graph(graph, ytitle, output_path, zero_line=False):
    canvas = ROOT.TCanvas(f"c_{graph.GetName()}", "", 900, 700)
    canvas.SetLogx()
    graph.SetMarkerStyle(20)
    graph.SetMarkerColor(ROOT.kBlue + 1)
    graph.SetLineColor(ROOT.kBlue + 1)
    graph.SetTitle(f";upper-track p_{{T}} [GeV];{ytitle}")
    graph.Draw("AP")
    if zero_line:
        canvas.Update()
        line = ROOT.TLine(
            graph.GetXaxis().GetXmin(),
            0.0,
            graph.GetXaxis().GetXmax(),
            0.0,
        )
        line.SetLineStyle(2)
        line.SetLineColor(ROOT.kGray + 2)
        line.Draw()
    canvas.SaveAs(output_path)


def make_pair_resolution_plots(chain, selection, outdir):
    pt_bins = [20.0, 30.0, 40.0, 50.0, 65.0, 85.0, 120.0, 200.0, 1000.0]
    residual = (
        "((1.0/evt_dsa_pt_lower)-(1.0/evt_dsa_pt_upper))"
        "/(1.0/evt_dsa_pt_upper)"
    )

    inclusive = make_hist(
        chain,
        "h_upper_lower_resolution_residual_inclusive",
        residual,
        selection,
        200,
        -5.0,
        5.0,
    )
    inclusive_fit = fit_single_gaussian(inclusive, "fit_upper_lower_resolution_residual_inclusive")
    canvas = ROOT.TCanvas("c_upper_lower_resolution_residual_inclusive", "", 900, 700)
    inclusive.SetLineColor(ROOT.kBlue + 1)
    inclusive.SetLineWidth(2)
    inclusive.GetXaxis().SetTitle("upper/lower |q|/p_{T} relative residual")
    inclusive.GetYaxis().SetTitle("Events")
    inclusive.Draw("HIST")
    if inclusive_fit:
        inclusive_fit.Draw("SAME")
    canvas.SaveAs(os.path.join(outdir, "upper_lower_resolution_residual_inclusive.png"))

    histograms = []
    fits = []
    mean_graph = ROOT.TGraphErrors()
    sigma_graph = ROOT.TGraphErrors()
    mean_graph.SetName("g_upper_lower_residual_mean_vs_upper_pt")
    sigma_graph.SetName("g_upper_lower_residual_sigma_vs_upper_pt")
    graph_point = 0
    for index, (low, high) in enumerate(zip(pt_bins[:-1], pt_bins[1:])):
        bin_selection = f"({selection}) && evt_dsa_pt_upper>={low} && evt_dsa_pt_upper<{high}"
        hist = make_hist(
            chain,
            f"h_upper_lower_resolution_residual_pt_{index}",
            residual,
            bin_selection,
            200,
            -5.0,
            5.0,
        )
        fit = fit_single_gaussian(hist, f"fit_upper_lower_resolution_residual_pt_{index}")
        histograms.append(hist)
        fits.append(fit)
        if fit:
            center = 0.5 * (low + high)
            half_width = 0.5 * (high - low)
            mean_graph.SetPoint(graph_point, center, fit.GetParameter(1))
            mean_graph.SetPointError(graph_point, half_width, fit.GetParError(1))
            sigma_graph.SetPoint(graph_point, center, fit.GetParameter(2))
            sigma_graph.SetPointError(graph_point, half_width, fit.GetParError(2))
            graph_point += 1

    draw_residual_grid(
        histograms,
        fits,
        pt_bins,
        os.path.join(outdir, "upper_lower_resolution_residual_by_upper_pt.png"),
    )
    draw_summary_graph(
        mean_graph,
        "Gaussian mean of relative residual",
        os.path.join(outdir, "upper_lower_residual_mean_vs_upper_pt.png"),
        zero_line=True,
    )
    draw_summary_graph(
        sigma_graph,
        "Gaussian #sigma of relative residual",
        os.path.join(outdir, "upper_lower_residual_sigma_vs_upper_pt.png"),
    )
    print(f"Quality-matched geometric-pair residual entries: {inclusive.GetEntries():.0f}")
    return [inclusive, inclusive_fit, *histograms, *fits, mean_graph, sigma_graph]


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
    objects.extend(draw_pair_momentum_diagnostics(chain, resolution_selection, args.outdir))

    mean_graph, count_graph, scan_hists = threshold_scan(
        chain, [20.0, 30.0, 40.0, 50.0], resolution_selection, args.outdir
    )
    objects.extend(scan_hists)
    objects.extend([mean_graph, count_graph])

    objects.extend(make_pair_resolution_plots(chain, resolution_selection, args.outdir))

    output_root = ROOT.TFile(os.path.join(args.outdir, "dsa_bias_study.root"), "RECREATE")
    for obj in objects:
        if obj:
            obj.Write()
    output_root.Close()
    print(f"Wrote DSA bias-study plots and ROOT objects to {args.outdir}")


if __name__ == "__main__":
    main()
