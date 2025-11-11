// --------------------------------------------------------------
//  plot_pdf_thesis.C
//  - LHAPDF comparison (replica / (a)symm-hessian handled)
//  - OTHERPARAM-augmented bands for baseline graphs
//  - Ratio (middle) + Relative-uncertainty (bottom) panels
// --------------------------------------------------------------

//root -l -b -q 'plot_pdf_thesis.C("uv",10,"output_HERA_wzProduction/plots.root","uv","output_HERA_wzProduction_WpJ+WmJ_ptw/plots.root","ATLASepWZVjet20-EIG","output_HERA_wzProduction_NegativeGluon/plots.root","output_HERA_wzProduction_WpJ+WmJ_ptw_NegativeGluon/plots.root")'
//root -l -b -q 'plot_pdf_thesis.C("uv",10,"output_HERA_wzProduction/plots.root","uv","output_HERA_wzProduction_WpJ+WmJ_ptw/plots.root","ATLASepWZVjet20-EIG")'
//root -l -b -q 'plot_pdf_thesis.C("uv",10,"output_HERA_wzProduction/plots.root","uv","output_HERA_wzProduction_WpJ+WmJ_ptw/plots.root","output_HERA_wzProduction_NegativeGluon/plots.root","output_HERA_wzProduction_WpJ+WmJ_ptw_NegativeGluon/plots.root")'

// to test the reference:
// root -l -b -q 'plot_pdf_thesis.C("sbar",10,"output_HERA_wzProduction_NegativeGluon/plots.root","sbar"," ","ATLAS-epWZ16-EIG","output_HERA_WZ/plots.root"," ")'

//Q2VAL = 1.9, 3.0, 4.0, 5., 10., 100., 6464, 8317 


#include "Rtypes.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TLatex.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <map>

// ---------- LHAPDF helper --------------------------------------
#include "LHAPDF/LHAPDF.h"

// ====================== FLAVOUR EXPRESSIONS =======================================================
enum class Op { Linear, Ratio };           // Linear = sum_i c_i * pdf(flavor_i) ; Ratio = (Linear)/(Linear)
struct LinTerm { int pdg; double coeff; }; // PDG: 21=g, +-1..+-6 quarks (no photon by default)
struct PidCombo { std::vector<LinTerm> terms; };
struct PidExpr  { Op op; PidCombo num, den; const char* note; bool lh_ok; };

// Evaluate x*f(x,Q2) for a linear combo in one member 'm'
inline double eval_combo(const std::vector<LHAPDF::PDF*>& pdfs, size_t m,
                         const PidCombo& pc, double x, double Q2){
  double s=0.0; for(const auto& t: pc.terms) s += t.coeff * pdfs[m]->xfxQ2(t.pdg, x, Q2);
  return s;
}

// Evaluate expression (linear or ratio) in member 'm'
inline double eval_expr(const std::vector<LHAPDF::PDF*>& pdfs, size_t m,
                        const PidExpr& px, double x, double Q2){
  if(px.op==Op::Linear) return eval_combo(pdfs, m, px.num, x, Q2);
  const double N = eval_combo(pdfs, m, px.num, x, Q2);
  const double D = eval_combo(pdfs, m, px.den, x, Q2);
  return (D!=0.0 ? N/D : 0.0);
}

// Catalogue mirroring the ROOT keys
inline const std::map<std::string,PidExpr>& pid_catalogue(){
  static std::map<std::string,PidExpr> M;

  auto L = [](std::initializer_list<LinTerm> t, const char* note)->PidExpr{
    return {Op::Linear, PidCombo{std::vector<LinTerm>(t)}, {}, note, true};
  };
  auto R = [](std::initializer_list<LinTerm> n,
              std::initializer_list<LinTerm> d,
              const char* note)->PidExpr{
    return {Op::Ratio, PidCombo{std::vector<LinTerm>(n)}, PidCombo{std::vector<LinTerm>(d)}, note, true};
  };
  auto BASELINE_ONLY = [](const char* note)->PidExpr{
    return {Op::Linear, PidCombo{{}}, {}, note, false};
  };

  // Elementary and valence (present in file)
  M["uv"]   = L({{+2,+1},{-2,-1}},             "u - ubar");
  M["dv"]   = L({{+1,+1},{-1,-1}},             "d - dbar");
  M["g"]    = L({{21,+1}},                     "gluon");
  M["ubar"] = L({{-2,+1}},                     "anti-up");
  M["dbar"] = L({{-1,+1}},                     "anti-down");
  M["s"]    = L({{+3,+1}},                     "strange");
  M["sbar"] = L({{-3,+1}},                     "anti-strange");
  M["c"]    = L({{+4,+1}},                     "charm (commonly c+cbar in plots)");
  M["b"]    = L({{+5,+1}},                     "bottom (commonly b+bbar in plots)");

  // Type sums (HERA/xFitter naming in your files)
  M["U"]    = L({{+2,+1},{+4,+1}},             "u + c");
  M["D"]    = L({{+1,+1},{+3,+1},{+5,+1}},             "d + s + b");
  M["UBar"] = L({{-2,+1},{-4,+1}},             "ubar + cbar");
  M["DBar"] = L({{-1,+1},{-3,+1},{-5,+1}},             "dbar + sbar + bbar");
  M["Sea"]  = L({{-2,+1},{-1,+1},{-3,+1},{-4,+1}}, "ubar + dbar + sbar + cbar");

  // Linear combos present explicitly
  M["uv+dv"]        = L({{+2,+1},{-2,-1},{+1,+1},{-1,-1}}, "uv + dv");
  M["uv+dv+2Sea"]   = L({{+2,+1},{+1,+1},                   // uv+dv
                         {-2,+1},{-1,+1},{-3,+2},{-4,+2}},  // + 2*Sea
                        "uv + dv + 2*Sea");
  M["23uv+13dv"]    = L({{+2, 2.0/3},{-2,-2.0/3},{+1,1.0/3},{-1,-1.0/3}},
                        "(2/3)uv + (1/3)dv");
  M["uv-dv"]        = L({{+2,+1},{-2,-1},{+1,-1},{-1,+1}}, "uv - dv");
  M["dbar-ubar"]    = L({{-1,+1},{-2,-1}},                 "dbar - ubar");

  // Ratios present in file
  M["doveru"]        = R({{+1,+1}}, {{+2,+1}},                  "d/u");
  M["dbaroverubar"]  = R({{-1,+1}}, {{-2,+1}},                  "dbar/ubar");
  M["dvoveruv"]      = R({{+1,+1},{-1,-1}}, {{+2,+1},{-2,-1}}, "dv/uv");
  M["goversea"]      = R({{21,+1}}, {{-2,+1},{-1,+1},{-3,+1},{-4,+1}}, "g/Sea");
  M["soversbar"]     = R({{+3,+1}}, {{-3,+1}},                  "s/sbar");

  //Rs
  M["Rs"] = R({{+3,+1},{-3,+1}},   // s + sbar
            {{-2,+1},{-1,+1}},   // ubar + dbar
            "R_{s}=(s+#bar{s})/(#bar{u}+ #bar{d})");

  // Ambiguous / analysis-specific keys → allow baseline, disable overlays
  M["ph"]  = BASELINE_ONLY("photon or analysis-specific; overlay disabled");
  M["sg"]  = BASELINE_ONLY("analysis-specific; overlay disabled");
  M["gg"]  = BASELINE_ONLY("analysis-specific; overlay disabled");
  M["rs"]  = BASELINE_ONLY("strangeness fraction; overlay disabled");

  // Standard PDF singlet over {u,d,s,c,b}: Sigma = sum_q (q + qbar)
  M["Sigma"] = L({{+2,+1},{-2,+1},{+1,+1},{-1,+1},{+3,+1},{-3,+1},{+4,+1},{-4,+1},{+5,+1},{-5,+1}},
                 "quark singlet (u..b): Σ=∑(q+qbar)");

  return M;
}

// Parse user key
inline PidExpr decode(const std::string& key){
  const auto& M = pid_catalogue();
  auto it = M.find(key);
  if(it==M.end()) throw std::runtime_error("Unknown flavour '"+key+"'");
  return it->second;
}

// ========================================================================

struct MuSig { double mu, sigP, sigM; };

// Uncertainty evaluator (handles linear/ratio and all error types)
MuSig mu_sigmas_Q2(const LHAPDF::PDFSet& set,
                   const std::vector<LHAPDF::PDF*>& pdfs,
                   const PidExpr& px, double x, double Q2)
{
  const std::string et = set.errorType();
  const bool isRep = et.find("replica")!=std::string::npos;
  const bool isSym = et.find("symmhessian")!=std::string::npos;

  // 68% CL rescale if needed
  const double fac = isRep ? 1.0 : set.errorInfo().conflevel/68.268949;

  const size_t first = isRep ? 1 : 0, last = pdfs.size();
  const double mu = eval_expr(pdfs, 0, px, x, Q2);

  if(isRep){
    double s2=0.0;
    for(size_t m=first; m<last; ++m){
      const double y = eval_expr(pdfs, m, px, x, Q2);
      const double d = y - mu; s2 += d*d;
    }
    const double sig = std::sqrt(s2 / std::max<size_t>(1, (last-first-1))) * fac;
    return {mu, sig, sig};
  }

  if(isSym){
    double s2=0.0;
    for(size_t p=1; p+1<pdfs.size(); p+=2){
      const double up = eval_expr(pdfs, p,   px, x, Q2);
      const double dn = eval_expr(pdfs, p+1, px, x, Q2);
      const double d = 0.5*(up-dn); s2 += d*d;
    }
    const double sig = std::sqrt(s2) * fac;
    return {mu, sig, sig};
  }

  // Asymmetric Hessian (PDG)
  double sP=0.0, sM=0.0;
  for(size_t p=1; p+1<pdfs.size(); p+=2){
    const double up = eval_expr(pdfs, p,   px, x, Q2);
    const double dn = eval_expr(pdfs, p+1, px, x, Q2);
    const double dP = std::max({0.0, up-mu, dn-mu});
    const double dM = std::max({0.0, mu-up, mu-dn});
    sP += dP*dP; sM += dM*dM;
  }
  return {mu, std::sqrt(sP)*fac, std::sqrt(sM)*fac};
}

// ---- Helpers to read graphs and (if needed) build Sigma from components ----
inline TGraphAsymmErrors* get_graph(TFile& f, const TString& key){
  TDirectory* g = (TDirectory*)f.Get("Graphs");
  if(!g) return nullptr;
  TGraphAsymmErrors* gr=nullptr; g->GetObject(key,gr); return gr;
}

inline TGraphAsymmErrors* linCombGraphs(const std::vector<const TGraphAsymmErrors*>& gs,
                                        const std::vector<double>& coeffs){
  if(gs.empty() || gs.size()!=coeffs.size()) return nullptr;
  for(const auto* g: gs) if(!g) return nullptr;
  const int N = gs[0]->GetN();
  for(const auto* g: gs) if(g->GetN()!=N) return nullptr;

  auto* out = (TGraphAsymmErrors*)gs[0]->Clone();
  for(int i=0;i<N;++i){
    const double x = gs[0]->GetX()[i];
    double y=0.0, eL2=0.0, eH2=0.0;
    for(size_t j=0;j<gs.size();++j){
      const double c  = coeffs[j];
      const double yj = gs[j]->GetY()[i];
      const double l  = gs[j]->GetErrorYlow(i);
      const double h  = gs[j]->GetErrorYhigh(i);
      y   += c * yj;
      eL2 += c*c * l*l;
      eH2 += c*c * h*h;
    }
    out->SetPoint(i, x, y);
    out->SetPointError(i, 0, 0, std::sqrt(eL2), std::sqrt(eH2));
  }
  return out;
}

// Build Sigma = (uv + 2*ubar) + (dv + 2*dbar) + (s + sbar) + (c + cbar) + (b + bbar)
// Falls back to 2*c and/or 2*b if cbar/bbar graphs are missing.
inline TGraphAsymmErrors* buildSigmaIfMissing(TFile& f, double q2){
  auto g_uv = get_graph(f, Form("dir1_q2_%g_pdf_uv",   q2));
  auto g_dv = get_graph(f, Form("dir1_q2_%g_pdf_dv",   q2));
  auto g_ub = get_graph(f, Form("dir1_q2_%g_pdf_ubar", q2));
  auto g_db = get_graph(f, Form("dir1_q2_%g_pdf_dbar", q2));
  auto g_s  = get_graph(f, Form("dir1_q2_%g_pdf_s",    q2));
  auto g_sb = get_graph(f, Form("dir1_q2_%g_pdf_sbar", q2));
  auto g_c  = get_graph(f, Form("dir1_q2_%g_pdf_c",    q2));
  auto g_cb = get_graph(f, Form("dir1_q2_%g_pdf_cbar", q2)); // may be null
  auto g_b  = get_graph(f, Form("dir1_q2_%g_pdf_b",    q2));
  auto g_bb = get_graph(f, Form("dir1_q2_%g_pdf_bbar", q2)); // may be null

  // Require the light pieces to exist
  if(!(g_uv && g_dv && g_ub && g_db && g_s && g_sb && g_c && g_b)){
    std::cout << "[Sigma] Missing one or more required light/heavy graphs at Q2="<< q2 <<"\n";
    return nullptr;
  }

  std::vector<const TGraphAsymmErrors*> gs;
  std::vector<double> coeffs;

  // uv + 2*ubar
  gs.push_back(g_uv); coeffs.push_back(1.0);
  gs.push_back(g_ub); coeffs.push_back(1.0);
  gs.push_back(g_ub); coeffs.push_back(1.0);

  // dv + 2*dbar
  gs.push_back(g_dv); coeffs.push_back(1.0);
  gs.push_back(g_db); coeffs.push_back(1.0);
  gs.push_back(g_db); coeffs.push_back(1.0);

  // s + sbar
  gs.push_back(g_s ); coeffs.push_back(1.0);
  gs.push_back(g_sb); coeffs.push_back(1.0);

  // c + cbar  (fallback: 2*c)
  if(g_cb){
    gs.push_back(g_c ); coeffs.push_back(1.0);
    gs.push_back(g_cb); coeffs.push_back(1.0);
  }else{
    std::cout << "[Sigma] No cbar graph found → using 2*c at Q2="<< q2 <<"\n";
    gs.push_back(g_c ); coeffs.push_back(2.0);
  }

  // b + bbar  (fallback: 2*b)
  if(g_bb){
    gs.push_back(g_b ); coeffs.push_back(1.0);
    gs.push_back(g_bb); coeffs.push_back(1.0);
  }else{
    std::cout << "[Sigma] No bbar graph found → using 2*b at Q2="<< q2 <<"\n";
    gs.push_back(g_b ); coeffs.push_back(2.0);
  }

  return linCombGraphs(gs, coeffs);
}

// Pretty TLatex labels for axis titles (concise but readable)
inline TString latex_pdf_label(const TString& key){
  static const std::map<std::string,std::string> L = {
    {"uv","u_{v}"}, {"dv","d_{v}"},
    {"g","g"},
    {"ubar","#bar{u}"}, {"dbar","#bar{d}"},
    {"s","s"}, {"sbar","#bar{s}"},
    {"c","c"}, {"b","b"},
    {"U","u + c"}, {"D","d + s"},
    {"UBar","#(bar{u} + #bar{c})"}, {"DBar","#bar{D}"},
    {"Sea","Sea"},
    {"uv+dv","u_{v} + d_{v}"},
    {"uv+dv+2Sea","u_{v} + d_{v} + 2\\,Sea"},
    {"23uv+13dv","#frac{2}{3}u_{v} + #frac{1}{3}d_{v}"},
    {"uv-dv","u_{v} - d_{v}"},
    {"dbar-ubar","#bar{d} - #bar{u}"},
    {"doveru","#frac{d}{u}"},
    {"dbaroverubar","#frac{#bar{d}}{#bar{u}}"},
    {"dvoveruv","#frac{d_{v}}{u_{v}}"},
    {"goversea","#frac{g}{Sea}"},
    {"soversbar","#frac{s}{#bar{s}}"},
    {"ph","#gamma"},
    {"sg","sg"}, {"gg","gg"},
    {"Rs","R_{s}"}, {"rs","r_{s}"},
    {"Sigma","#Sigma"}
  };
  auto it = L.find(key.Data());
  return (it!=L.end()? TString(it->second.c_str()) : key);
}

// ----------------------------------------------------------------------
void plot_pdf_thesis( TString flav="dv", double q2=10.,
                      TString file1="output/plots.root",
                      TString flav2=" ", TString file2=" ",
                      // Default: disable LHAPDF unless explicitly passed
                      TString csvSets=" ",
                      // Default: no OTHERPARAM unless explicitly passed
                      TString file1_other=" ", TString file2_other=" " )
{
  // style zero borders (helps remove hairline gaps)
  gStyle->SetOptTitle(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);

  // --- argument disambiguation ---------------------------------------
  auto looksLikeFile = [](const TString& s){
    return s.EndsWith(".root") || s.Contains("/") || s.Contains("\\");
  };
  if( looksLikeFile(csvSets) && looksLikeFile(file1_other) && file2_other==" " ){
    file2_other = file1_other;
    file1_other = csvSets;
    csvSets = " ";
  }

  const bool wantOther  = (file1_other!=" " || file2_other!=" ");
  const bool wantLHAPDF = (csvSets!=" " && csvSets!="" );

  const PidExpr pid = decode(flav.Data());

  // ---------- parse LHAPDF set list (optional) -----------------------
  std::vector<std::string> sets;
  if(wantLHAPDF){
    std::stringstream ss(csvSets.Data()); std::string s;
    while(std::getline(ss,s,',')) if(!s.empty()) sets.push_back(s);
  }

  auto nicelabel = [](TString s){
    s.ReplaceAll("/plots.root",""); s.ReplaceAll("output_","");
    return s;
  };

  // ---------- open baseline ------------------------------------------
  TFile f1(file1);  TDirectory* g1=(TDirectory*)f1.Get("Graphs");
  if(!g1){ std::cout<<"No Graphs in "<<file1<<"\n"; return; }
  TString key1; key1.Form("dir1_q2_%g_pdf_%s",q2,flav.Data());
  TGraphAsymmErrors* pref=(TGraphAsymmErrors*)g1->Get(key1);
  // Fallback build for Sigma if missing
  if(!pref && flav=="Sigma"){
    pref = buildSigmaIfMissing(f1, q2);
    if(pref){
      pref->SetName("Sigma_built");
      pref->SetTitle("");
    }
  }
  if(!pref){ std::cout<<"Cannot find "<<key1<<"\n"; return; }
  pref->SetTitle("");
  pref->SetFillColor(kRed-9); pref->SetFillStyle(1001); pref->SetFillColorAlpha(kRed-9,0.65);
  pref->SetLineColor(kRed);   pref->SetLineWidth(1);

  // ---------- OTHERPARAM for baseline (optional) ---------------------
  TGraphAsymmErrors* pref_other=nullptr;  TFile *f1o_ptr=nullptr;   // (explicit type to appease Cling)
  if(wantOther && file1_other!=" "){
    f1o_ptr = new TFile(file1_other);
    if(!f1o_ptr->Get("Graphs")){ std::cout<<"No Graphs in "<<file1_other<<"\n"; return; }
    pref_other = get_graph(*f1o_ptr, key1);
    if(!pref_other && flav=="Sigma"){ // build Sigma in OTHERPARAM file too
      pref_other = buildSigmaIfMissing(*f1o_ptr, q2);
      if(pref_other){ pref_other->SetName("Sigma_built_other"); pref_other->SetTitle(""); }
    }
    if(!pref_other){ std::cout<<"Cannot find "<<key1<<" in "<<file1_other<<"\n"; return; }
  }
  TGraphAsymmErrors* pref_aug = nullptr;
  if(wantOther){
    pref_aug = (TGraphAsymmErrors*)pref->Clone("pref_aug");
    pref_aug->SetFillColor(kRed-3); pref_aug->SetFillStyle(3004); pref_aug->SetLineColor(kRed+1);
  }

  // ---------- optional second baseline + OTHERPARAM -------------------
  TGraphAsymmErrors* p2=nullptr;  TFile *f2_ptr=nullptr;
  TGraphAsymmErrors* p2_other=nullptr; TFile *f2o_ptr=nullptr;
  if(flav2!=" "&&file2!=" "){
     f2_ptr = new TFile(file2);
     if(!f2_ptr->Get("Graphs")){ std::cout<<"No Graphs in "<<file2<<"\n"; return; }
     TString key2; key2.Form("dir1_q2_%g_pdf_%s",q2,flav2.Data());
     p2=get_graph(*f2_ptr,key2);
     if(!p2 && flav2=="Sigma"){
       p2 = buildSigmaIfMissing(*f2_ptr, q2);
       if(p2){ p2->SetName("Sigma_built_2"); p2->SetTitle(""); }
     }
     if(!p2){ std::cout<<"Cannot find "<<key2<<"\n"; return; }
     p2->SetTitle("");
     p2->SetFillColor(kBlue-9); p2->SetFillStyle(3001); p2->SetFillColorAlpha(kBlue-9,0.50);
     p2->SetLineColor(kBlue);   p2->SetLineWidth(1);
     if(wantOther && file2_other!=" "){
       f2o_ptr = new TFile(file2_other);
       if(!f2o_ptr->Get("Graphs")){ std::cout<<"No Graphs in "<<file2_other<<"\n"; return; }
       p2_other = get_graph(*f2o_ptr,key2);
       if(!p2_other && flav2=="Sigma"){
         p2_other = buildSigmaIfMissing(*f2o_ptr, q2);
         if(p2_other){ p2_other->SetName("Sigma_built_2_other"); p2_other->SetTitle(""); }
       }
       if(!p2_other){ std::cout<<"Cannot find "<<key2<<" in "<<file2_other<<"\n"; return; }
     }
  }
  TGraphAsymmErrors* p2_aug = nullptr;
  if(wantOther && p2){ p2_aug=(TGraphAsymmErrors*)p2->Clone("p2_aug");
                       p2_aug->SetFillColor(kBlue-4); p2_aug->SetFillStyle(3005); p2_aug->SetLineColor(kBlue+1); }

  const int N=pref->GetN();  // grid
  std::vector<int> lhcol={kGreen+2,kMagenta+2,kOrange+7,kAzure+7,kViolet+1};

  // ---------- LHAPDF overlays (skip if flavour flagged baseline-only) ------
  std::vector<TGraphAsymmErrors*> glha; glha.reserve(sets.size());
  if(wantLHAPDF && !pid.lh_ok){
    std::cout << "LHAPDF overlay disabled for flavour '" << flav << "' (baseline-only key)\n";
  }
  if(wantLHAPDF && pid.lh_ok){
    for(size_t s=0;s<sets.size();++s){
       LHAPDF::PDFSet pdfset(sets[s]);
       const auto& pdfs=pdfset.mkPDFs();
       auto* g=new TGraphAsymmErrors(N);
       for(int i=0;i<N;++i){
          double x=pref->GetX()[i];
          const MuSig ms = mu_sigmas_Q2(pdfset,pdfs,pid,x,q2);
          g->SetPoint(i,x,ms.mu);
          g->SetPointError(i,0,0,ms.sigM,ms.sigP);
       }
       g->SetFillColorAlpha(lhcol[s%lhcol.size()],0.35);
       g->SetLineColor(lhcol[s%lhcol.size()]); g->SetLineWidth(2);
       glha.push_back(g);
    }
  }

  // ---------- ratio & rel-unc graphs ---------------------------------
  TGraphAsymmErrors* r  = new TGraphAsymmErrors(N);
  TGraphAsymmErrors* r2 = p2? new TGraphAsymmErrors(N):nullptr;
  std::vector<TGraphAsymmErrors*> rlha(glha.size());
  for(size_t s=0;s<glha.size();++s) rlha[s]=new TGraphAsymmErrors(N);

  TGraphAsymmErrors* r_aug  = wantOther ? new TGraphAsymmErrors(N) : nullptr;
  TGraphAsymmErrors* r2_aug = (wantOther && p2)? new TGraphAsymmErrors(N):nullptr;

  TGraphAsymmErrors* u    = new TGraphAsymmErrors(N);
  TGraphAsymmErrors* u2   = p2? new TGraphAsymmErrors(N):nullptr;
  TGraphAsymmErrors* u_aug  = wantOther ? new TGraphAsymmErrors(N) : nullptr;
  TGraphAsymmErrors* u2_aug = (wantOther && p2)? new TGraphAsymmErrors(N):nullptr;

  std::vector<TGraphAsymmErrors*> ulha(glha.size());
  for(size_t s=0;s<glha.size();++s) ulha[s]=new TGraphAsymmErrors(N);

  for(int i=0;i<N;++i){
     const double x   = pref->GetX()[i];
     const double ref = pref->GetY()[i];
     const double eLr = pref->GetErrorYlow(i), eHr = pref->GetErrorYhigh(i);
     const double ref_other = (pref_other ? pref_other->GetY()[i] : ref);
     const double d_ref     = std::abs(ref_other - ref);
     const double eLr_tot   = std::sqrt(eLr*eLr + d_ref*d_ref);
     const double eHr_tot   = std::sqrt(eHr*eHr + d_ref*d_ref);

     if(pref_aug) pref_aug->SetPointError(i,0,0,eLr_tot,eHr_tot);
     if(p2 && p2_aug){
       const double n2  = p2->GetY()[i];
       const double eL2 = p2->GetErrorYlow(i), eH2 = p2->GetErrorYhigh(i);
       const double n2_other = (p2_other ? p2_other->GetY()[i] : n2);
       const double d2 = std::abs(n2_other - n2);
       const double eL2_tot = std::sqrt(eL2*eL2 + d2*d2);
       const double eH2_tot = std::sqrt(eH2*eH2 + d2*d2);
       p2_aug->SetPointError(i,0,0,eL2_tot,eH2_tot);
     }

     r->SetPoint(i,x,1.0);
     r->SetPointError(i,0,0, (ref? eLr/ref:0.), (ref? eHr/ref:0.));
     if(r_aug){
       r_aug->SetPoint(i,x,1.0);
       r_aug->SetPointError(i,0,0, (ref? eLr_tot/ref:0.), (ref? eHr_tot/ref:0.));
     }

     if(p2){
        const double num   = p2->GetY()[i];
        const double eL2   = p2->GetErrorYlow(i), eH2 = p2->GetErrorYhigh(i);
        const double ratio = (ref? num/ref : 0.0);

        r2->SetPoint(i,x,ratio);
        r2->SetPointError(i,0,0, (ref? eL2/ref:0.), (ref? eH2/ref:0.));

        if(r2_aug){
           const double n2_other = (p2_other ? p2_other->GetY()[i] : num);
           const double d2 = std::abs(n2_other - num);
           const double eL2_tot = std::sqrt(eL2*eL2 + d2*d2);
           const double eH2_tot = std::sqrt(eH2*eH2 + d2*d2);
           r2_aug->SetPoint(i,x,ratio);
           r2_aug->SetPointError(i,0,0, (ref? eL2_tot/ref:0.), (ref? eH2_tot/ref:0.));
        }
     }

     for(size_t s=0;s<glha.size();++s){
        const auto* g = glha[s];
        const double num = g->GetY()[i];
        const double eL  = g->GetErrorYlow(i), eH=g->GetErrorYhigh(i);
        rlha[s]->SetPoint(i,x, ref? num/ref : 0.0);
        rlha[s]->SetPointError(i,0,0, ref? eL/ref : 0.0, ref? eH/ref : 0.0);

        ulha[s]->SetPoint(i,x,1.0);
        ulha[s]->SetPointError(i,0,0, num? eL/num : 0.0, num? eH/num : 0.0);
     }

     u->SetPoint(i,x,1.0);
     u->SetPointError(i,0,0, (ref? eLr/ref:0.), (ref? eHr/ref:0.));
     if(u_aug){
       u_aug->SetPoint(i,x,1.0);
       u_aug->SetPointError(i,0,0, (ref? eLr_tot/ref:0.), (ref? eHr_tot/ref:0.));
     }

     if(p2){
        const double num = p2->GetY()[i];
        const double eL2 = p2->GetErrorYlow(i), eH2 = p2->GetErrorYhigh(i);
        u2->SetPoint(i,x,1.0);
        u2->SetPointError(i,0,0, (num? eL2/num:0.), (num? eH2/num:0.));
        if(u2_aug){
          const double n2_other = (p2_other ? p2_other->GetY()[i] : num);
          const double d2 = std::abs(n2_other - num);
          const double eL2_tot = std::sqrt(eL2*eL2 + d2*d2);
          const double eH2_tot = std::sqrt(eH2*eH2 + d2*d2);
          u2_aug->SetPoint(i,x,1.0);
          u2_aug->SetPointError(i,0,0, (num? eL2_tot/num:0.), (num? eH2_tot/num:0.));
        }
     }
  }

  // ---------- style ---------------------------------------------------
  r->SetFillColorAlpha(kRed-9,0.40); r->SetLineColor(kRed);
  if(r_aug){ r_aug->SetFillColor(kRed-3);  r_aug->SetFillStyle(3004);  r_aug->SetLineColor(kRed+1); }
  if(r2){ r2->SetFillColor(kBlue-9); r2->SetFillStyle(3001); r2->SetLineColor(kBlue); }
  if(r2_aug){ r2_aug->SetFillColor(kBlue-4); r2_aug->SetFillStyle(3005); r2_aug->SetLineColor(kBlue+1); }
  for(size_t s=0;s<rlha.size();++s){ rlha[s]->SetFillColorAlpha(lhcol[s%lhcol.size()],0.35);
                                     rlha[s]->SetLineColor(lhcol[s%lhcol.size()]); }
  for(size_t s=0;s<ulha.size();++s){ ulha[s]->SetFillColorAlpha(lhcol[s%lhcol.size()],0.35);
                                     ulha[s]->SetLineColor(lhcol[s%lhcol.size()]); }

  // ---------- canvas and pads -----------------------------------------
  const double relsize = 0.25;    // bottom
  const double ratsize = 0.25;    // middle
  const double top_frac = 1.0 - ratsize - relsize;

  const double SCALE = 1.45;
  const double label_base  = 0.10 * SCALE;
  const double title_base  = 0.12 * SCALE;

  const double top_label = label_base * (ratsize / top_frac);
  const double top_title = title_base * (ratsize / top_frac);
  const double bot_label = label_base * (ratsize / relsize);
  const double bot_title = title_base * (ratsize / relsize);

  const double yoff_mid = 0.36;
  const double yoff_top = yoff_mid * (title_base / top_title);
  const double yoff_bot = yoff_mid * (title_base / bot_title);

  TCanvas* c=new TCanvas("PDF","pdf",550,690);

  const double left = 0.12, right = 0.02;

  // ---------- TOP PAD -------------------------------------------------
  TPad* p1=new TPad("p1","p1",0.,ratsize+relsize,1.,1.);
  p1->SetLogx(); p1->SetLeftMargin(left); p1->SetRightMargin(right);
  p1->SetTopMargin(0.02);  p1->SetBottomMargin(0.025);
  p1->SetTickx(); p1->SetTicky();
  p1->Draw(); p1->cd();

  pref->GetXaxis()->Set(101,1e-4,1.);
  {
    TString ytitle = "#it{x} " + latex_pdf_label(flav);
    pref->GetYaxis()->SetTitle(ytitle);
  }
  pref->GetYaxis()->SetLabelSize(top_label);
  pref->GetYaxis()->SetTitleSize(top_title);
  pref->GetYaxis()->SetTitleOffset(yoff_top);
  pref->GetXaxis()->SetLabelSize(0);

  pref->Draw("ACE3");
  if(pref_aug) pref_aug->Draw("CE3 SAME");
  if(p2){ p2->Draw("CE3 SAME"); if(p2_aug) p2_aug->Draw("CE3 SAME"); }
  for(auto* g:glha) g->Draw("CE3 SAME");

  // Legend
  //TLegend* leg = new TLegend(0.14, 0.62, 0.54, 0.92);
  //TLegend* leg = new TLegend(0.14, 0.14, 0.54, 0.44); // bottom-left box
  //TLegend* leg = new TLegend(0.30, 0.14, 0.70, 0.44); // bottom-center box
  //TLegend* leg = new TLegend(0.38, 0.14, 0.78, 0.44); // bottom-center, a bit further right
  TLegend* leg = new TLegend(0.40, 0.55, 0.80, 0.85); // upper-right, a bit left & lower 


  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->SetTextSize(0.06);
  leg->SetMargin(0.1);
  leg->SetNColumns(1);
  TString lab1=nicelabel(file1);
  leg->AddEntry(pref,lab1.Data(),"f");
  if(pref_aug) leg->AddEntry(pref_aug,(lab1+" total unc.").Data(),"f");
  if(p2){
     TString lab2=nicelabel(file2);
     leg->AddEntry(p2,lab2.Data(),"f");
     if(p2_aug) leg->AddEntry(p2_aug,(lab2+" total unc.").Data(),"f");
  }
  for(size_t s=0;s<glha.size();++s) leg->AddEntry(glha[s],sets[s].c_str(),"f");
  leg->Draw();

  TLatex tex; tex.SetNDC(); tex.SetTextAlign(34); tex.SetTextSize(top_label);
  tex.DrawLatex(0.92,0.9,Form("Q^{2} = %.2f GeV^{2}",q2));

  // ---------- MIDDLE PAD (ratio) -------------------------------------
  c->cd();
  TPad* p2pad=new TPad("p2","p2",0.,relsize,1.,ratsize+relsize);
  p2pad->SetLogx(); p2pad->SetLeftMargin(left); p2pad->SetRightMargin(right);
  p2pad->SetTopMargin(0.00); p2pad->SetBottomMargin(0.00);
  p2pad->SetTickx(); p2pad->SetTicky();
  p2pad->Draw(); p2pad->cd();

  r->SetMaximum(1.35); r->SetMinimum(0.65); r->GetYaxis()->SetNdivisions(504);
  r->GetYaxis()->SetTitle("Ratio");
  r->GetYaxis()->SetLabelSize(label_base);
  r->GetYaxis()->SetTitleSize(title_base);
  r->GetYaxis()->SetTitleOffset(yoff_mid);
  r->GetXaxis()->Set(101,1e-4,1.);
  r->GetXaxis()->SetLabelSize(0); r->GetXaxis()->SetTitleSize(0);

  r->Draw("ACE3");
  if(r_aug)   r_aug->Draw("CE3 SAME");
  if(r2)      r2->Draw("CE3 SAME");
  if(r2_aug)  r2_aug->Draw("CE3 SAME");
  for(auto* g:rlha) g->Draw("CE3 SAME");

  // ---------- BOTTOM PAD (relative uncertainties) --------------------
  c->cd();
  TPad* p3=new TPad("p3","p3",0.,0.,1.,relsize);
  p3->SetLogx(); p3->SetLeftMargin(left); p3->SetRightMargin(right);
  p3->SetTopMargin(0.00);
  p3->SetBottomMargin(0.40);
  p3->SetTickx(); p3->SetTicky();
  p3->Draw(); p3->cd();

  u->SetMaximum(1.35); u->SetMinimum(0.65); u->GetYaxis()->SetNdivisions(504);
  u->GetYaxis()->SetTitle("Rel. unc.");
  u->GetYaxis()->SetLabelSize(bot_label);
  u->GetYaxis()->SetTitleSize(bot_title);
  u->GetYaxis()->SetTitleOffset(yoff_bot);

  u->GetXaxis()->Set(101,1e-4,1.);
  u->GetXaxis()->SetLabelSize(bot_label);
  u->GetXaxis()->SetTitle("#it{x}");
  u->GetXaxis()->SetTitleSize(bot_title);
  u->GetXaxis()->SetTitleOffset(1.00);

  u->SetFillColor (kRed-9);  u->SetFillStyle (1001); u->SetLineColor (kRed);
  if(u_aug){ u_aug->SetFillColor (kRed-3); u_aug->SetFillStyle (3004); u_aug->SetLineColor (kRed+1); }
  if(u2){ u2->SetFillColor (kBlue-9); u2->SetFillStyle (3001); u2->SetLineColor (kBlue); }
  if(u2_aug){ u2_aug->SetFillColor (kBlue-4); u2_aug->SetFillStyle (3005); u2_aug->SetLineColor (kBlue+1); }

  u->Draw("ACE3");
  if(u_aug)   u_aug->Draw("CE3 SAME");
  if(u2)      u2->Draw("CE3 SAME");
  if(u2_aug)  u2_aug->Draw("CE3 SAME");
  for(auto* g:ulha) g->Draw("CE3 SAME");

  // ---------- save and cleanup -----------------------------------------
  c->SaveAs(Form("pdf_%s_Q2_%g_withLHAPDF_relunc.pdf",flav.Data(),q2));

  if(f1o_ptr){ f1o_ptr->Close(); delete f1o_ptr; }
  if(f2_ptr){  f2_ptr->Close();  delete f2_ptr;  }
  if(f2o_ptr){ f2o_ptr->Close(); delete f2o_ptr; }
}
