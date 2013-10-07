#ifndef LATINO_STYLE_H
#define LATINO_STYLE_H
#include "TStyle.h"
#include "TColor.h"

TStyle *LatinoStyle() {

  TStyle* GloStyle;
  GloStyle = gStyle;

  TStyle* ls = new TStyle("LatinosStyle", "LatinosStyle");
  gStyle = ls;


  //----------------------------------------------------------------------------
  // Canvas
  //----------------------------------------------------------------------------
  ls->SetCanvasBorderMode(  0);
  ls->SetCanvasBorderSize( 10);
  ls->SetCanvasColor     (  0);
  ls->SetCanvasDefH      (600);
  ls->SetCanvasDefW      (550);
  ls->SetCanvasDefX      ( 10);
  ls->SetCanvasDefY      ( 10);


  //----------------------------------------------------------------------------
  // Pad
  //----------------------------------------------------------------------------
  ls->SetPadBorderMode  (   0);
  ls->SetPadBorderSize  (  10);
  ls->SetPadColor       (   0);
  ls->SetPadBottomMargin(0.20);
  ls->SetPadTopMargin   (0.08);
  ls->SetPadLeftMargin  (0.18);
  ls->SetPadRightMargin (0.05);


  //----------------------------------------------------------------------------
  // Frame
  //----------------------------------------------------------------------------
  ls->SetFrameFillStyle ( 0);
  ls->SetFrameFillColor ( 0);
  ls->SetFrameLineColor ( 1);
  ls->SetFrameLineStyle ( 0);
  ls->SetFrameLineWidth ( 2);
  ls->SetFrameBorderMode( 0);
  ls->SetFrameBorderSize(10);


  //----------------------------------------------------------------------------
  // Hist
  //----------------------------------------------------------------------------
  ls->SetHistFillColor(0);
  ls->SetHistFillStyle(1);
  ls->SetHistLineColor(1);
  ls->SetHistLineStyle(0);
  ls->SetHistLineWidth(1);


  //----------------------------------------------------------------------------
  // Axis
  //----------------------------------------------------------------------------
  ls->SetLabelFont  (   42, "xyz");
  ls->SetLabelOffset(0.015, "xyz");
  ls->SetLabelSize  (0.050, "xyz");
  ls->SetNdivisions (  505, "xyz");
  ls->SetTitleFont  (   42, "xyz");
  ls->SetTitleSize  (0.050, "xyz");

  //  ls->SetNdivisions ( -503, "y");

  ls->SetTitleOffset(  1.4,   "x");
  ls->SetTitleOffset(  1.2,   "y");
  ls->SetPadTickX   (           1);  // Tick marks on the opposite side of the frame
  ls->SetPadTickY   (           1);  // Tick marks on the opposite side of the frame


  //----------------------------------------------------------------------------
  // Title
  //----------------------------------------------------------------------------
  ls->SetTitleBorderSize(    0);
  ls->SetTitleFillColor (   10);
  ls->SetTitleAlign     (   12);
  ls->SetTitleFontSize  (0.040);
  ls->SetTitleX         (0.590);
  ls->SetTitleY         (0.860);

  ls->SetTitleFont(42, "");


  //----------------------------------------------------------------------------
  // Stat
  //----------------------------------------------------------------------------
  ls->SetOptStat       (1110);
  ls->SetStatBorderSize(   0);
  ls->SetStatColor     (  10);
  ls->SetStatFont      (  42);
  ls->SetStatX         (0.94);
  ls->SetStatY         (0.91);

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  ls->SetNumberContours(NCont);

  return ls;
}
#endif
