#!/bin/bash

#cat log/* | grep HOGE2 | grep -v nan | awk '{if( $1 != 0 && $2 != 0 && $1/$2 < 2 ) print $1" "$2" "$1/$2" "(1.01811*$1-$2)/sqrt(1.01868*1.01868+1)}' > tmp_2nd.log
#cat log/* | grep HOGE3 | grep -v nan | awk '{if( $1 != 0 && $2 != 0 && $1/$2 < 2 ) print $1" "$2" "$1/$2" "(1.01811*$1-$2)/sqrt(1.00340*1.00340+1)}' > tmp_3rd.log


root -l <<EOF
TGraph* g = new TGraph("tmp_3rd.log","%*lg %*lg %lg %lg");
TH1D* h = new TH1D("h","h",100, 0.1, 0.1 );
for( int i=0; i<g->GetN(); i++ ) h->Fill( g->GetY()[i] );
h->Draw();
h->Fit("gaus");
c1->Print("pic/beta_3rd_dev.eps");
EOF

root -l <<EOF
TGraph* g = new TGraph("tmp_3rd.log");
g->SetMarkerStyle(20);
g->SetMarkerSize(0.3);
g->Draw("AP");
g->Fit("pol1");
TGraph* g_SM = new TGraph();
g_SM->SetMarkerStyle(20);
g_SM->SetMarkerSize(0.6);
g_SM->SetMarkerColor(2);
//g_SM->SetPoint(0, 0.110736, 0.131560); # 2nd bin
g_SM->SetPoint(0, 0.319769, 0.321939); # 3rd bin
g_SM->Draw("Psame");
c1->Print("pic/beta_3rd.eps");
EOF

exit

for A7  in 100
do
for A9  in `seq -200 20 200`
#for A9  in 100 50 -50
do
for A10 in `seq -200 20 200`
#for A10 in 100 50 -50
do
   bsub -q e ./theory_curve_dev    $A7  $A9  $A10
done
done
done
