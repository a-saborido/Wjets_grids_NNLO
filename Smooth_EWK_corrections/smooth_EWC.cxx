// smooth_EWC.cxx  ---------------------------------------------------
#include <cmath>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraphSmooth.h"

/*
// ------------------------------------------------------------------
// 1) manual 5-bin weighted smoother
// ------------------------------------------------------------------
TH1D* SmoothHistogramWithErrors(const TH1D* h)
{
   const int nBins = h->GetNbinsX();
   auto* h_smoothed = (TH1D*)h->Clone("h_manual");
   h_smoothed->Reset();

   for (int i = 1; i <= nBins; ++i) {
      double sumWeightedContent = 0., sumWeights = 0., sumErr2 = 0.;
      for (int j = i - 2; j <= i + 2; ++j) {
         if (j < 1 || j > nBins) continue;
         const double c = h->GetBinContent(j);
         const double e = h->GetBinError  (j);
         const double w = h->GetBinWidth  (j);
         if (e <= 0. || w <= 0.) continue;
         const double weight = 1.0/(e*e*w);
         sumWeightedContent += weight * c;
         sumWeights         += weight;
         sumErr2            += weight * weight * e * e;
      }
      if (sumWeights > 0.) {
         h_smoothed->SetBinContent(i, sumWeightedContent / sumWeights);
         h_smoothed->SetBinError  (i, std::sqrt(sumErr2) / sumWeights);
      }
   }
   return h_smoothed;
}
*/

// ------------------------------------------------------------------
// 1) kernel–weighted 5-bin smoother (built-in kernel x 1/sigma^2 weights)
// ------------------------------------------------------------------
TH1D* SmoothHistogramWithErrors(const TH1D* h)
{
   const int nBins = h->GetNbinsX();

   // kernel (1 2 3 2 1)/9  – the one TH1::Smooth() uses
   constexpr int    KSIZE = 5;
   constexpr int    HALF  = KSIZE / 2;             // = 2
   const     double K[KSIZE] = {1., 2., 3., 2., 1.};
   constexpr double KNORM = 9.0;

   auto* h_smoothed = (TH1D*)h->Clone("h_manual");
   h_smoothed->Reset("ICE");

   for (int i = 1; i <= nBins; ++i) {

      double num   = 0.0;   // sum w · y
      double denom = 0.0;   // sum w
      double varW  = 0.0;   // sum w^2 · sigma^2

      for (int j = i - HALF; j <= i + HALF; ++j) {
         if (j < 1 || j > nBins) continue;         // skip off-hist bins

         const int    kIdx = j - i + HALF;         // 0…4
         const double kVal = K[kIdx] / KNORM;

         const double y  = h->GetBinContent(j);
         const double sigma  = h->GetBinError  (j);
         if (sigma <= 0.0) continue;

         const double w = kVal / (sigma*sigma);            // kernel x 1/sigma^2

         num   += w * y;
         denom += w;
         varW  += w * w * sigma * sigma;                   // variance propagation
      }

      if (denom > 0.0) {
         const double ySmooth = num / denom;
         const double errSmooth = std::sqrt(varW) / denom;

         h_smoothed->SetBinContent(i, ySmooth);
         h_smoothed->SetBinError  (i, errSmooth);
      }
   }

   return h_smoothed;
}

// ------------------------------------------------------------------
// 2) statistically-weighted Super-Smoother (ROOT’s TGraphSmooth)
// ------------------------------------------------------------------
TH1D* SmoothHistogramSuper(const TH1D* h,
                           double bass = 0.0,
                           double span = 0.0)
{
   const int nBins = h->GetNbinsX();

   auto* g = new TGraphErrors(nBins);
   std::vector<double> weights(nBins);

   for (int i = 1; i <= nBins; ++i) {
      const double x  = h->GetBinCenter (i);
      const double y  = h->GetBinContent(i);
      const double dy = h->GetBinError  (i);

      g->SetPoint      (i-1, x, y);
      g->SetPointError (i-1, 0.0, dy);
      weights[i-1] = (dy > 0.0) ? 1.0 / (dy*dy) : 0.0;
   }

   TGraphSmooth gs;
   TGraph* gSmooth = gs.SmoothSuper(g, "", bass, span, kFALSE, weights.data());

   auto* hSmooth = (TH1D*)h->Clone("h_super");
   hSmooth->Reset("ICE");

   for (int i = 0; i < nBins; ++i) {
      double xx, yy;
      gSmooth->GetPoint(i, xx, yy);
      hSmooth->SetBinContent(i+1, yy);
      hSmooth->SetBinError  (i+1, 0.0);
   }

   delete g;    // input graph
   return hSmooth;
}

// ------------------------------------------------------------------
// 3) main macro
// ------------------------------------------------------------------
void smooth_EWC()
{
   gStyle->SetOptStat(0);

   // -------- user parameters --------------------------------------
   const std::string histName = "jet_y1_1j_Wp";          // choose histogram: jet_y1_1j_Wp, jet_pt1_1j_Wp, HT_1j_Wp, W_pt_1j_Wp
   const std::string histPath = "EWC/" + histName;       // directory in file
   const double YMIN = 0.2, YMAX = 1.17;                 // display range
   // ----------------------------------------------------------------

   TFile* f = TFile::Open("EWC.root");
   if (!f || f->IsZombie()) { printf("Cannot open EWC.root\n"); return; }

   TH1D* h_raw = nullptr;
   f->GetObject(histPath.c_str(), h_raw);
   if (!h_raw) { printf("Histogram %s not found.\n", histPath.c_str()); return; }

   // 2) smoothers ---------------------------------------------------
   TH1D* h_manual = SmoothHistogramWithErrors(h_raw);
   TH1D* h_super  = SmoothHistogramSuper    (h_raw);

   // 2b) average ----------------------------------------------------
   TH1D* h_avg = (TH1D*)h_manual->Clone("h_avg");
   h_avg->Reset("ICE");

   const int nBins = h_avg->GetNbinsX();
   for (int i = 1; i <= nBins; ++i) {
      const double yM = h_manual->GetBinContent(i);
      const double yS = h_super ->GetBinContent(i);

      const double yA = 0.5 * (yM + yS);
      const double eA = std::fabs(yM - yS);
      h_avg->SetBinContent(i, yA);
      h_avg->SetBinError  (i, eA);
   }

   // 3) styling -----------------------------------------------------
   h_raw   ->SetMarkerStyle(20); h_raw   ->SetMarkerSize(0.8);
   h_manual->SetMarkerStyle(24); h_manual->SetMarkerSize(0.8);
   h_super ->SetMarkerStyle(22); h_super ->SetMarkerSize(0.8);
   h_avg   ->SetMarkerStyle(21); h_avg   ->SetMarkerSize (0.8);

   h_raw   ->SetLineColor(kBlack);   h_raw   ->SetMarkerColor(kBlack);
   h_manual->SetLineColor(kBlue+1);  h_manual->SetMarkerColor(kBlue+1);
   h_super ->SetLineColor(kGreen+2); h_super ->SetMarkerColor(kGreen+2);
   h_avg   ->SetLineColor(kMagenta+2); h_avg->SetMarkerColor(kMagenta+2);

   // 4) canvas & axes ----------------------------------------------
   auto* c = new TCanvas("c","RAW vs smoothed",1200,850);
   c->SetGrid();
   c->SetLeftMargin(0.12);
   c->SetRightMargin(0.05);
   c->SetBottomMargin(0.13);
   c->SetTopMargin(0.05);

   h_raw->SetTitle("");
   h_raw->GetXaxis()->SetTitle(histName.c_str());
   h_raw->GetXaxis()->SetTitleOffset(1.1);
   h_raw->GetYaxis()->SetTitleOffset(1.4);
   h_raw->GetXaxis()->SetTitleSize(0.045);
   h_raw->GetYaxis()->SetTitleSize(0.045);
   h_raw->GetXaxis()->SetLabelSize(0.040);
   h_raw->GetYaxis()->SetLabelSize(0.040);

   h_raw->SetMinimum(YMIN);
   h_raw->SetMaximum(YMAX);

   // 5) draw --------------------------------------------------------
   h_raw   ->Draw("E1");
   h_manual->Draw("E1 SAME");
   h_super ->Draw("E1 SAME");
   h_avg   ->Draw("E1 SAME");

   // legend
   auto* leg = new TLegend();
   leg->SetTextSize(0.032);
   leg->AddEntry(h_raw,    "raw data"                  ,"lep");
   leg->AddEntry(h_manual, "kernel 5-bin weighted"     ,"lep");
   leg->AddEntry(h_super,  "Super-Smoother (1/sigma^{2})"  ,"lep");
   leg->AddEntry(h_avg,    "average (manual+super)/2"  ,"lep");
   leg->SetBorderSize(0);
   leg->Draw();

   c->SaveAs((histName + "_smoothed.pdf").c_str());

   // 6) terminal dumps ---------------------------------------------
   auto dump = [](const char* tag, const TH1D* h)
   {
      const int n = h->GetNbinsX();
      printf("\n%s\n%-14s %-14s %-14s\n",
             tag,"low","high","content");
      for (int i = 1; i <= n; ++i) {
         const double low  = h->GetBinLowEdge(i);
         const double high = low + h->GetBinWidth(i);
         const double val  = h->GetBinContent(i);
         printf("%14.6g %14.6g %14.6g\n", low, high, val);
      }
   };

   dump("Kernel-weighted smoother",        h_manual);
   dump("Super-Smoother (1/sigma²)",           h_super );
   dump("Average (manual + super)/2",      h_avg   );

   // 7) difference manual − super ---------------------------------
   printf("\nDifference  (manual − super)\n%-14s %-14s %-14s\n",
          "low","high","Δcontent");

   for (int i = 1; i <= nBins; ++i) {
      const double low  = h_manual->GetBinLowEdge(i);
      const double high = low + h_manual->GetBinWidth(i);
      const double delta = abs(h_manual->GetBinContent(i) -
                           h_super ->GetBinContent(i))*100;
      //printf("%14.6g %14.6g %14.6g\n", low, high, delta);
      printf("%14.6g\n", delta);
   }
}
